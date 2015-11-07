function l=isodd(v)
% function l=isodd(v)
%
% is the argument even, odd or else
%
% input  :	v		: data array
%
% output :	l		: -1 if v is not a natural number
%           			   0 if v is even
%	    			   1 if v is odd
%
% version 1.0.0		last change 04.09.1995

% Gerd Krahmann, IfM Kiel, Aug 1993

if (length(v)>1)
  dv=v-fix(v);
  bad=find(dv~=0);
  s=size(v);
  l=ones(s(1),s(2));
  l(bad)=-ones(1,length(bad));
  v2=v/2;
  ev=find(fix(v/2)==v/2);
  l(ev)=zeros(1,length(ev));
else
  l=0;
  if fix(v)<v,
    l=-1;
  end
  if fix(v/2)<v/2,
    l=1;
  end
end
