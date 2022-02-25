%% RAnalyzer performs the Roughness calculation of an audio signal. 
%
% Calculation of Roughness using the ERB scale  described in:
% Estimation of Perceived Roughness
% V.J.P Jourdes - 29 July 2004
%
% This model is based on the optimized Aures described in: 
% P. Daniel and R. Weber, “Psychoacoustic Roughness: Implementation of an 
% Optimized Model,” Acustica 83, 113~123 (1997).
%
% Usage:
% [R] = RAnalyzer(audioSignal)
%
% Input parameters:
%   audioSignal  -  column array of 16000 audio samples
%
% Output parameters:
%   R  - roughness in aspers.
%
%********************************************************************************
% Copyright (c) 2018-2019 Department of Computer Science                        *
%                         ETSE, Universitat de València                         *
%                         46100, Burjassot, Valencia, Spain                     *
%                                                                               *
% This file is part of the URPAA: URBAURAMON Psycho-acoustic Annoyance Analizer *
%                                                                               *
% URPAA is free software:  you can redistribute it and/or modify it  under      *
% the terms of the  GNU  General  Public  License  as published by the  Free    *
% Software Foundation, either version 3 of the License,  or (at your option)    *
% any later version.                                                            *
%                                                                               *
% URPAA  is distributed in the hope that it will be useful, but WITHOUT ANY     *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS    *
% FOR A PARTICULAR PURPOSE.                                                     *
% See the GNU General Public License for more details.                          *
%                                                                               *
% You should  have received a copy  of the GNU General Public License  along    *
% with this program.  If not, see <http://www.gnu.org/licenses/>.               *
%                                                                               *
% https://github.com/jausegar/urbauramon/tree/master/URPAA  jaume.segura@uv.es  *
%********************************************************************************

function R  = RAnalyzer(xr,a,b,c,d,Nch,fcoefs,gzi,fso) %#codegen

%load('Rpar.mat', 'a', 'b', 'c', 'd',  'Nch', 'fcoefs', 'filt', 'gzi', 'fs');

%coder.extrinsic('AdaptLevel','ERBFilterBank');

 fs = 16000.;
 %if(size(xr,1)==1)
 %    xr = xr'; 
 %end
     
 if(fso~=fs)
     xs = resample1(xr,fs,fso);
 else
     xs = xr(:);
 end
 
 
xs = AdaptLevel(xs,60);

winlen = size(xs,2);
wlen   = winlen*fs;
win_= blackman(wlen);
win = win_;
xf = enframe(xs,win);
nn=conv(b, d);
dd=conv(a, c);
signal = filter(nn, dd, xf');

exc = zeros(76,16000);
exc = ERBFilterBank(signal, fcoefs);
exc = max(exc,0);

etmp=zeros(76,16000);
etmp = abs(exc);

h0=zeros(76,16000);
h0 = kron(ones(1,length(exc)), mean(etmp,2)); %<-- Valor DC

exdc=zeros(76,16000);
excd = etmp - h0;

hBPi=zeros(76,16000);
numH7=1e-06 *[-0.921451216121396   0.982888560318938   0.295260688176911  -0.286857741098071  -0.069840689764262];
denH7=[1.000000000000000  -5.899617218896303  14.513555995920429 -19.057265213007142  14.086668700939789  -5.557679124645693   0.914336863008540];
hBPi(1:11,:) = filter(numH7,denH7,excd(1:11,:), [], 2);
%filt('numH7'), filt('denH7'), excd(1:11,:), [], 2);
numH14=[0.000306736050582  -0.000939253934583   0.001106419814007  -0.000593687072990  0.000119785307609];
denH14=[1.000000000000000  -5.298993260459181  11.670499756774747 -13.672869088042340   8.986920269391314  -3.142229897723366   0.456672622395994];
hBPi(12:28,:) = filter(numH14,denH14,excd(12:28,:), [], 2);
%filt('numH14'), filt('denH14'), excd(12:28,:), [], 2);
numH30 = 1e-04 *[0.129727267205170  -0.147102163717761   0.011692413564215  -0.038565888535768   0.044284366915656];
denH30 = [1.000000000000000  -5.801221661924524  14.053701947571739 -18.198219047957338  13.284982096699689  -5.183997259147179   0.844754071932504];
hBPi(29:36,:) = filter(numH30,denH30,excd(29:36,:), [], 2);
%filt('numH30'), filt('denH30'), excd(29:36,:), [], 2);
numH36=1e-04 * [0.305545894165424  -0.482784630733298   0.177214106193544   0.000053943523087   0.000000815705521];
denH36=[1.000000000000000  -5.700385282122440  13.570085247681225 -17.268667566010468  12.389932821492810  -4.752259715846983   0.761294748827704];
hBPi(37:65,:) = filter(numH36,denH36,excd(37:65,:), [], 2);
%filt('numH36'), filt('denH36'), excd(37:65,:), [], 2);
numH66=1e-04 *[0.149958281527685  -0.042586768984454  -0.110078886192767   0.001578826028454   0.001285865584613];
denH66=[1.000000000000000  -5.400367006426416  12.113760050230690 -14.435484703808234   9.628677305908159  -3.404071762415771   0.497486522472457];
hBPi(66:76,:) = filter(numH66,denH66,excd(66:76,:), [], 2);
%filt('numH66'), filt('denH66'), excd(66:76,:), [], 2);
hBPrms = rms(hBPi,2); 

rexc=zeros(76,1);
rexc = rms(exc,2);
maxi=0;
maxi = max(rexc);
if maxi > 0
    calib = rexc / maxi;
else
    calib = 0;
end

mdepth = zeros(1,Nch);

for n = 1:Nch
    if h0(n) > 0
        mdepth(n) = hBPrms(n) / h0(n);
        mdepth(n) = mdepth(n) * calib(n);
    else
        mdepth(n) = 0;
    end
end
ki = zeros(1,Nch);

for n = 1:Nch
    if (n < Nch - 1)
        amount = 0.003 * fs;
        ki(n) = shiftcov(hBPi(n,:), hBPi(n+2,:), amount);
    end
end
                
Cr = 1.1998;  
ri = zeros(1,Nch);
ri(1:7) = (gzi(1:7)' .* mdepth(1:7) .* ki(1:7)).^2;
ri(8:Nch-3) = (gzi(8:Nch-3)'.* mdepth(8:Nch-3).* ki(6:Nch-5).* ki(8:Nch-3)).^2;
ri(Nch-2) = (gzi(Nch-2) * mdepth(Nch-2) * ki(Nch-4))^2;
ri(Nch-1) = (gzi(Nch-1) * mdepth(Nch-1) * ki(Nch-3))^2;
ri(Nch) = (gzi(Nch) * mdepth(Nch) * ki(Nch-2))^2;
R = Cr*sum(ri);

end

