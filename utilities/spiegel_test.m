function [H,pval] = spiegel_test(x, alpha)
% FUNCTION [H, pval] = spiegel_test(x,alpha)
% compute p- value under null hypothesis of x normally distributed;
% x should be a vector; alpha is a scalar (test-level, default = 0.05)
%  
% H = 0 indicates a failure to reject the null hypothesis at the alpha significance level
%  (in less crypto-statistical parlance, it means the series is Gaussian as far as we can tell)
% H = 1 means the series is non-Gaussianby this test.
% 
% Ref: D. J. Spiegelhalter, 'Diagnostic tests of distributional shape,' Biometrika, 1983
% =====================================================
if nargin < -1
   alpha = 0.05;
end

xm = mean(x);
xs = std(x);
xz = (x - xm) ./ xs;
xz2 = xz.^2;
N = sum(xz2 .* log(xz2));
n = numel(x);
ts = (N - 0.73 * n) / (0.8969 * sqrt(n)); %under the null, ts ~ N(0,1)
pval = 1 - abs(erf(ts / sqrt(2)));    %2-sided test.

H = logical(pval < 0.05);
