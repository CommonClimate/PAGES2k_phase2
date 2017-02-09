function y = nsum(x,dim)
%NSUM Sum of elements with missing data.
%  Y = NSUM(X) returns the sum of each column of X as a row vector
%  where missing data values are encoded as NaNs. For vectors, NSUM(X)
%  returns the sum of the elements in X.
%  By default NSUM averages the first dimension of non-zero size.
%  Y = NSUM(X,DIM) sums up the dimension DIM
%

%  C. Mertens, IfM Kiel
%  $Revision: 1.2 $ $Date: 1996/01/15 13:24:41 $
%  added compatibility to MATLAB 5

if nargin<2
  dim=min(find(size(x)>1));
end
if isempty(dim)
  dim=1;
end
y = zeros(size(x));
valid = ~isnan(x);
y(valid) = x(valid);
y = sum(y,dim);
