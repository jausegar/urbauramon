%% FAnalyzer performs the Roughness calculation of an audio signal. 
%
% Calculation of Fluctuaton Strength using the ERB scale  described in:
% Roughness - Fluctuation strength
% Salford Innovation Research Centre (SIRC)
%
% This model is based on the optimized Aures described in: 
% Psychoacoustical Roughness: Implementation of an Optimized Model
% P. Daniel and R. Weber, “Psychoacoustic Roughness: Implementation of an 
% Optimized Model,” Acustica 83, 113~123 (1997).
%
% This model is based on the Zwicker Model:
% Zwicker E., Fastl H. ‘Psychoacoustics: Facts and Models’(1990).
%
% Usage:
% [F] = FAnalyzer(audioSignal)
%
% Input parameters:
%   audioSignal  -  column array of 16000 audio samples
%
% Output parameters:
%   F  - Fluctuation Strength in vacils.
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
function [F] = FAnalyzer(xr)

load('Fpar.mat','b', 'a', 'd', 'c', 'Nch', 'fcoefs', 'filt', 'gzi','fs');

xr = xr.';
xr = AdaptLevel(xr,60); 
signal = filter(conv(b, d), conv(a, c), xr);
exc = ERBFilterBank(signal, fcoefs);
exc = max(exc,0);

etmp = abs(exc);
h0 = kron(ones(1,length(exc)), mean(etmp,2)); 
excd = etmp - h0;

hBPi = zeros(Nch,size(xr,2));
for n = 1:Nch
    hBPi(n,:) = filter(filt, excd(n,:));
end
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

mdepth(mdepth > 1) = 1;

ki = zeros(1,Nch);
for n = 1:Nch
    if (n < Nch - 1)
        amount = 0.003 * fs;
        ki(n) = shiftcov(hBPi(n,:), hBPi(n+2,:), amount); %<--- Le cuesta
    end
end

Cf = 0.3423;  
fi = zeros(1,Nch);
fi(1:2) = (gzi(1:2) .* mdepth(1:2) .* ki(1:2) .* 1).^2;
fi(3:Nch-2) = (gzi(3:Nch-2) .* mdepth(3:Nch-2) .* ki(3:Nch-2) .* ki(1:Nch-4)).^2;
fi(Nch-1:Nch) = (gzi(Nch-1:Nch) .* mdepth(Nch-1:Nch) .* 1 .* ki(Nch-1:Nch)).^2;
F = Cf*sum(fi);

end
