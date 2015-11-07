function [H,P] =  symmetry_test(X,alpha)
% function [H,P] =  symmetry_test(X,alpha)
% Tests for the symmetry of a distribution, as in Lazante et al [1996]. 
% inputs: X : data vector 
%         alpha : test level [default = 0.05]
% outputs: 
%  - H, takes value 1 if test is rejected (i.e. the data are significantly
%         skewed), 0 otherwise
%  - P, test p-value
%
%  Reference: 
%  Lanzante, J. R. (1996), Resistant, robust and non?paramteric 
%  techniques for the analysis of climate data: Theory and examples, including 
% applications to historical radiosonde station data, Int. J. Climatol.,
% 16, 1197?1226, doi:10.1002/(SICI)1097-0088(199611)16:11<1197::AID-JOC89>
% 3.0.CO;2-L.
%
% originally written by Deborah Khider (USC, 2011)
% edited by Julien Emile-Geay, 08/10/2013 
%  (streamlined description, included arbitrary test level & NaN immunity)
% ------------------------------------------------------

if nargin < 2
   alpha = 0.05; % default value
end

% control for NaNs
X   = X(~isnan(X));
% calculate the median
med = median(X);

% partition the X
XL = X(find(X>=med));
XG = X(find(X<med));
% reflect around the median
Y  = med + (med-XL);

%perform a rank test
[P,H] = ranksum(XG,Y,'alpha',alpha);
end