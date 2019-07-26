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

function [R] = RAnalyzer(xr)

load('Rpar.mat', 'a', 'b', 'c', 'd',  'Nch', 'fcoefs', 'filt', 'gzi', 'fs');

xr = xr.'; 
xr = AdaptLevel(xr,60);

winlen = size(xr,1);
wlen   = winlen*fs;
xf = enframe(xr,blackman(wlen));
signal = filter(conv(b, d), conv(a, c), xf);
exc = ERBFilterBank(signal, fcoefs);
exc = max(exc,0);

etmp = abs(exc);
h0 = kron(ones(1,length(exc)), mean(etmp,2)); %<-- Valor DC
excd = etmp - h0;

hBPi(1:11,:) = filter(filt('numH7'), filt('denH7'), excd(1:11,:), [], 2);
hBPi(12:28,:) = filter(filt('numH14'), filt('denH14'), excd(12:28,:), [], 2);
hBPi(29:36,:) = filter(filt('numH30'), filt('denH30'), excd(29:36,:), [], 2);
hBPi(37:65,:) = filter(filt('numH36'), filt('denH36'), excd(37:65,:), [], 2);
hBPi(66:76,:) = filter(filt('numH66'), filt('denH66'), excd(66:76,:), [], 2);
hBPrms = rms(hBPi,2); 

rexc = rms(exc,2);
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
ri(1:7) = (gzi(1:7) .* mdepth(1:7) .* ki(1:7)).^2;
ri(8:Nch-3) = (gzi(8:Nch-3) .* mdepth(8:Nch-3) .* ki(6:Nch-5) .* ki(8:Nch-3)).^2;
ri(Nch-2) = (gzi(Nch-2) * mdepth(Nch-2) * ki(Nch-4))^2;
ri(Nch-1) = (gzi(Nch-1) * mdepth(Nch-1) * ki(Nch-3))^2;
ri(Nch) = (gzi(Nch) * mdepth(Nch) * ki(Nch-2))^2;
R = Cr*sum(ri);

end

