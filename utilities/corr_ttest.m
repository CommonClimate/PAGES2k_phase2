function [r,signif,pval] = corr_ttest(x,y,alpha);
% CORR_TTEST - estimates the significance of correlations between 2
%             time series using the classical T-test with degrees of
%             freedom modified for autocorrelation.
%
%   [r,signif,F] = corr_ttest(x,y,alpha)
%
% This function creates 'nsim' random time series that have the same power
% spectrum as the original time series but with random phases.
%
% Input
% x, y : vector of (real) numbers of identical length
% alpha : [optional] significance level for critical value estimation [0.05]
%
% Output
% r [real] : correlation between x and y
% signif [boolean]: true (1) if significant; false (0) otherwise
% pval [real] : test p-value (the probability of the test statstic exceeding the
% observed one by chance alone) 
% =====================================

% Calculate the Pearson correlation coefficient
r = corr(x,y); 

% Find the effective degrees of freedom
g1  = ar1(x);  g2  = ar1(y);     % sample lag-1 correlation coefficient

N   = length(x);
Nex = N*(1-g1)/(1+g1);
Ney = N*(1-g2)/(1+g2);

Ne  = geomean(Nex + Ney);
if Ne < 10
   error('Too few effective d.o.f. to apply this method')
end

% degrees of freedom
df  = Ne - 2;
% Calculate the t statistic
t   = abs(r).*sqrt(df./(1-r.^2));

% Handle the NANs
t(isnan(t))=0;

% Calculate the p-values
% p = 2*(1-tcdf(abs(t),df));
pval = 2*tcdf(-abs(t),df);

% Determine statistical significance
signif = (pval <= alpha);
return