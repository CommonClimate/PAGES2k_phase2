function [y,mi] = nmean(x,dim)
%NMEAN  Average or mean value, ignoring NaN.
%   Same as MEAN, but NaN's are ignored.
%
% uses :	intvers.m

%   Copyright (c) 1997 by Toby Driscoll.
%   Adapted from MEAN.M, written by The MathWorks, Inc.
%   added backwards compatibility	G.Krahmann, LODYC Paris
%   removed bug on undoc. feature	GK 	28.5.1999

if prod(size(x))==1, 
  y = x; 
  if nargout==2
    mi = 1;
  end
  return 
end

if nargin==1, 
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

if intvers>4
  mask = isnan(x);
  if nargout==2
    x = shiftdim(x, dim-1);
    s = size(x);
    s2 = s;
    s2(2:length(s2)) = 1;
    s(1) = 1;
    mi = reshape([1:s2(1)],s2);
    mi = repmat(mi,s);
    mi(mask) = 0;
    mi = shiftdim(mi, ndims(x)-(dim-1));
    x = shiftdim(x, ndims(x)-(dim-1));
  end
  x(mask) = 0;
  s = sum(~mask,dim);
  s(s==0) = NaN;
  if ~isempty(x),
    y = sum(x,dim)./s;
  else
    y = NaN;
  end
  if nargout==2
    if ~isempty(x),
      mi = sum(mi,dim)./s;
    else
      mi = NaN;
    end
  end
else
  [m,n] = size(x);
  y = zeros(m,n);
  valid = ~isnan(x);
  y(valid) = x(valid);
  s = sum(valid);
  if ~isempty(x)
    y = sum(y) ./ s;
  else
    y = NaN;
  end
  if nargout==2
    mi = zeros(m,n);
    mi2 = [1:m]'*ones(1,n);
    mi(valid) = mi2(valid);	
    if ~isempty(x)
      mi = sum(mi) ./ s;
    else
      mi = NaN;
    end
  end
end
