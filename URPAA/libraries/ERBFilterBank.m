% function output = ERBFilterBank(x, fcoefs)
% Process an input waveform with a gammatone filter bank. This function 
% takes a single sound vector, and returns an array of filter outputs, one 
% channel per row.
%
% The fcoefs parameter, which completely specifies the Gammatone filterbank,
% should be designed with the MakeERBFilters function.  If it is omitted,
% the filter coefficients are computed for you assuming a 22050Hz sampling
% rate and 64 filters regularly spaced on an ERB scale from fs/2 down to 100Hz.
%
% Malcolm Slaney @ Interval, June 11, 1998.
% (c) 1998 Interval Research Corporation  
% Thanks to Alain de Cheveigne' for his suggestions and improvements.
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



function output = ERBFilterBank(x, fcoefs)


if nargin < 1
	error('Syntax: output_array = ERBFilterBank(input_vector[, fcoefs]);');
end

if nargin < 2
	fcoefs = MakeERBFilters(22050,64,100);
end

if size(fcoefs,2) ~= 10
	error('fcoefs parameter passed to ERBFilterBank is the wrong size.');
end

if size(x,2) < size(x,1)
	x = x';
end

A0  = fcoefs(:,1);
A11 = fcoefs(:,2);
A12 = fcoefs(:,3);
A13 = fcoefs(:,4);
A14 = fcoefs(:,5);
A2  = fcoefs(:,6);
B0  = fcoefs(:,7);
B1  = fcoefs(:,8);
B2  = fcoefs(:,9);
gain= fcoefs(:,10);	

output = zeros(size(gain,1), length(x));
for chan = 1: size(gain,1)
	y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
		   A2(chan)/gain(chan)], ...
				[B0(chan) B1(chan) B2(chan)], x);
	y2=filter([A0(chan) A12(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y1);
	y3=filter([A0(chan) A13(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y2);
	y4=filter([A0(chan) A14(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y3);
	output(chan, :) = y4;
end
