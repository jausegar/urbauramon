% Author: V.J.P. Jourdes
% This function calculates the maximum of the cross correlation between two
% signal calculated after shifting one of them of a certain number of
% samples.
%
% y = shifcov(f1, f2, amount);
%
% f1: signal 1
% f2: signal 2
% amount: defines the maximum number of samples that the signal will be
% shifted.
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
function y = shiftcov(f1, f2, amount)
    L = length(f1);
    ff2 = f2;
    x = 1;
    for i = 1:10:amount % The signal is shifted 10 by 10 samples for less
        cfac = cov(f1, f2);
        den = diag(cfac);
        den = sqrt(den * den');
        if den(2,1) > 0
            r(x) = cfac(2,1) / den(2,1);
        else
            r(x) = 0;
        end
        x = x + 1;
        f2(1) = [];
        f2(L) = 0;
    end
    for i = 1:10:amount
        f1 = [zeros(1,i) f1(1:L-i)];
        cfac = cov(f1, ff2);
        den = diag(cfac);
        den = sqrt(den * den');
        if den(2,1) > 0
            r(x) = cfac(2,1) / den(2,1);
        else
            r(x) = 0;
        end
        x = x + 1;
    end
    y = max(r);
end
