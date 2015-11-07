function Xn = gaussianize(X)
%   function Xn = gaussianize(X)
%   Transform each column of data matrix X to normality using the inverse
%   Rosenblatt transform. Tolerates NaNs.
%
% inspired by split.m in normal.m by Van Albada, S.J., Robinson P.A. (2006)
% Transformation of arbitrary distributions to the normal distribution with application to EEG
% test-retest reliability. J Neurosci Meth, doi:10.1016/j.jneumeth.2006.11.004
%
%  Written 26/06/2015 by Julien Emile-Geay (USC)


[n,p] = size(X);
Xn    = NaN(n,p);
for j = 1:p
    % Sort the data in ascending order and retain permutation indices
    Z = X(:,j); nz = ~isnan(Z);
    [sorted,index]  = sort(Z(nz));
    % Make 'rank' the rank number of each observation
    [x, rank]       = sort(index);
    % The cumulative distribution function
    CDF = rank./n - 1/(2*n);
    % Apply the inverse Rosenblatt transformation
    Xn(nz,j) = sqrt(2)*erfinv(2*CDF - 1);  % Xn is now normally distributed
end

end