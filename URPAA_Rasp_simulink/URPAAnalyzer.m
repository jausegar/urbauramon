
%% URPAAnalizer - URBAURAMON Psycho-acoustic Annoyance analizer
%
% URPAAnalyzer is a function to perform psycho-acoustic annoyance analisys of   
% an audio input file in *.wav format. 
%
% The Psycho-acoustic model used is based on the Zwicker Model:
% Zwicker E., Fastl H. ‘Psychoacoustics: Facts and Models’(1990).
%
% Input parameters:
%   audioSignal  -  string with *.wav filename or path to a *.wav file
%
% Output parameters:
%   PA  - Psycho-acoustic Annoyance indicator based in Psycho-acoustic 
%         parameters: Loudness, Sharpness, Roughness and Fluctuation strength
%
% It should be noted that only single-channel audio is accepted and that  
% it is resampled at 16Khz before analysis.
%
% Example of use:
%   PA = URPAAnalyzer('lowPAexample.wav')
%
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
% URPAAnalyzer is a function to perform psycho-acoustic annoyance analisys of   *
% an audio input file in *.wav format.                                          *
% The psycho-acoustic model used to perform the analysis is based on Zwicker's  *
% psycho-acoustic annoyance model.												*
%                                                                               *
% In case of use of the implemented code, please cite us.                       *
%                                                                               *
% https://github.com/jausegar/urbauramon/tree/master/URPAA  jaume.segura@uv.es  *
%********************************************************************************

function [PA,N,S,R,F,T] = URPAAnalyzer(xr,H,H16k,a,b,c,d,Nch,fcoefs,gzi,fsam) %%%%%%%%%#codegen
    %addpath('libraries');
    %addpath('libraries\libNS');
    %addpath('libraries\libR');
    %addpath('libraries\libF');
    %addpath('libraries\libTon');
    xr=xr(:)';
    
    %coder.extrinsic('audioread');
    % Create a Raspberry Pi object
    %r= raspixº();
    
    %audioFile='highPAsound.wav';
    %[xr,fsam] = audioread(audioFile);
    %fsam=44100.;
    %tmp=0.;
    %tmp=fsam;
    fs = 16000.;
    %tmp
    
    if (fsam ~= fs)
        xs = resample1(xr,fs,fsam);
    else
        xs = xr;
    end
    
    tic
    reH=real(H);
    imH=imag(H);
    [N, S] = NSAnalyzer(xs,reH,imH,fs);
    R = RAnalyzer(xs,a,b,c,d,Nch,fcoefs,gzi,fs);
    F = FAnalyzer(xs,b,a,d,c,Nch,fcoefs,gzi,fs);
%
%  G.-Q. Di, X.-W. Chen, K. Song, B. Zhou,Ch.-M. Pei, "Improvement of Zwicker’s psychoacoustic annoyance model aiming at tonal noises".
%  Applied Acoustics 105 (2016) pp.164–170
    reH16k=real(H16k);
    imH16k=imag(H16k);
    [T,ff,LX,WT,Wgr]=TonAnalyzer(xr,fsam,reH16k,imH16k);

    Stmp=0.;
    Ntmp=0.;
    Stmp=S;
    Ntmp=N;
    wS = (Stmp-1.75).* 0.25 .* log10(Ntmp+10);
    wFR = (2.18 ./ (N.^0.4)) .* ((0.4.*F) + (0.6 .* R)); 
    wT = (6.41./(N.^0.52))*T;
    if (T>=0.3)
        PA = N .* (1 + sqrt((wS.^2) + (wFR.^2) + (wT.^2)));
    else
        PA = N .* (1 + sqrt((wS.^2) + (wFR.^2)));
    end
    toc

end

