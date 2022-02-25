function f=enframe(x,win,inc) %#codegen
%ENFRAME split signal up into (overlapping) frames: one per row. F=(X,WIN,INC)

nx=length(x);
nwin=length(win);
if (nwin == 1)
   len_ = win;
   len=len_(1);
else
   len_ = nwin;
   len=len_(1);
end
if (nargin < 3)
   inc = len;
end
nli=nx-len+inc;
nf_ = max(fix(nli/inc),0);   % number of full frames
nf = nf_(1);
%nf = fix((nx-len+inc)/inc);
f=zeros(nf,len);
indf= inc*(0:(nf-1)).';
inds = (1:len);
f(:) = x(indf(:,ones(1,len))+inds(ones(nf,1),:));
if (nwin > 1)
    w = win(:)';
    f = f .* w(ones(nf,1),:);
end