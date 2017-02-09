function [r,signif,p] = corr_sig(x,y,options)
% CORR_SIG  [r,signif,p] = corr_sig(x,y,options)
% - estimates the significance of correlations between non IID 
%             time series by 3 independent methods
%  1) T-test where d.o.f are corrected for the
%        effect of serial correlation  (options.method = 'ttest')
%   2) options.method = 'isopersistent': AR(1) modeling of x and y. 
%   3) options.method = 'isospectral': phase randomization of original inputs.
%     (3 is the default)
%   
% Inputs
% ======
%      x, y : vector of (real) numbers of identical length
% 
%       options: structure specifying options [default]:
%        - nsim : the number of simulations [1000]
%        - level : significance level[0.05] 
%        - method : methods 1-3 above ['isospectral']
%
% Output
% =======
% r [real] : correlation between x and y
% signif [boolean]: true (1) if significant; false (0) otherwise
% p [real] : Fraction of time series with higher correlation coefficents than
% observed (approximates the p-value). Note that signif = 1 if and only if p <= level. 
%
%
% The T-test is parametric test, hence cheap but usually wrong except in idyllic circumstances.
% The others are non-parametric, but their computational requirements
% scales with nsim. 
% 
%  

x = x(:); y = y(:); % column vectors

% process options
if nargin < 3 || isempty(options)
    fopts   = [];
else
    fopts   = fieldnames(options);
end

if strmatch('method', fopts)
    method = options.method;
else
    method = 'isospectral';
end

if max(strcmp('alpha', fopts)) == 0
    alpha = 0.05;
else
    alpha = options.alpha;
    if ~isscalar(alpha) || alpha <= 0 || alpha >= 1
        error('alpha must be within the range [0,1]');
    end
    
end
if max(strcmp('nsim', fopts)) == 0
   nsim = 2000;
else
   nsim = options.nsim;
end

% x and y myst have equal number of rows
if size(x,1)~=size(y,1)
    error('x and y must have equal number of rows.');
end

% apply 1 of the 3 methods

switch method
   case {'ttest'}
      [r,signif,p] = corr_ttest(x,y,alpha);
   case{'isopersistent'}
      [r,signif,p] = corr_isopersist(x,y,alpha,nsim);
   case{'isospectral'}      
      [r,signif,p] = corr_isospec(x,y,alpha,nsim);
end

end

