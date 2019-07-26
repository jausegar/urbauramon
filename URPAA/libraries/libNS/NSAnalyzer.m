%% NSAnalyzer performs the Loudness and Sharpness calculation of an audio signal. 
%
% Usage:
% [N,S] = NSAnalyzer(audioSignal)
%
% Input parameters:
%   audioSignal  -  column array of 16000 audio samples
%
% Output parameters:
%   N  - loudness in sones.
%   S  - sharpness in accums
%
% Loudness calculation is based in DIN 45631 
% Copyright 2003 Aaron Hastings, Ray W. Herrick Laboratories and Purdue University 
% 
% Sharpness calculation is based in function  
% Claire Churchill(Sep 2004) obtained by fitting an 
% equation to data given in 'Psychoacoustics: Facts and Models'
% using MATLAB basic fitting function
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
% https://github.com/jausegar/urbauramon/URPAA            jaume.segura@uv.es    *
%********************************************************************************
function [N, S] = NSAnalyzer(xr)
load('H.mat','H');

wlen = size(xr,1);
Xf = fft(enframe(xr,hanning(wlen)).');
Xf = Xf(1:wlen/2+1,:);
Xf_spl = 88.14 + 20*log10(abs(Xf)/(wlen/4));

mxFlt=28;
YdBZinput=zeros(1,mxFlt);
for n=1:mxFlt
    YdBZinput(n)=10*log10(sum((10.^(Xf_spl'/10)).*(abs(H(n,:).^2))));
end
YdBZinput(YdBZinput<-70) = -70;
[N,lspec,~]=DIN45631(YdBZinput,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nl = length(lspec);
gz = ones(1,140);
z = 141:nl;
gz(z) = 0.00012*(z/10).^4-0.0056*(z/10).^3+0.1*(z/10).^2-0.81*(z/10)+3.5;
z = 0.1:0.1:(nl/10);

S = 0.11 * sum(lspec.*gz.*z.*0.1) / sum(lspec.*0.1);

end

