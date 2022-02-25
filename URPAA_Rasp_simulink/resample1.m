function yi=resample1(y,Fs2,Fs1)
%Fs1 = 44.1e3;
%Fs2 = 16e3;

% Load your data
%y = load('yourdata'); %y=sin(0:1/Fs1:1);
%coder.extrinsic('length');
tmp = 0;
tmp = length(y);
Ttime = (tmp-1)*(1/Fs1);

% X-Axis
x  = 0:1/Fs1:Ttime;
xi = 0:1/Fs2:Ttime;

% Interpolate
method = 'cubic';
yi = interp1(x,y,xi,method);
yi=yi';
%% or
 % Zero-order hold
%yi = zeros(length(xi),1);
%jj = 1;
%xi(1) = x(1);
%for ii=2:length(y)
%    % Update value
%    if (x(ii)>=xi(jj)),
%        yi(jj) = y(ii-1);
%        jj = jj+1;
%    end
%end
