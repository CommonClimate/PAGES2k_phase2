function [r,signif,F] = corr_isospec(x,y,level,nsim)

% CORR_ISOSPEC: estimates the significance of correlations between non IID 
%             time series by phase randomization of original inputs.
%
%   [r,signif,F] = corr_isospec(x,y,level,nsim)
%
% This function creates 'nsim' random time series that have the same power
% spectrum as the original time series but random phases.
%
% Input
%  x, y : vector of (real) numbers of identical length
%  level : [optional] significance level for critical value estimation [0.05]
%  nsim : [optional] the number of simulations [1000]
%
% Output
%  r [real] : correlation between x and y
%  signif [boolean]: true (1) if significant; false (0) otherwise
%  F : Fraction of time series with higher correlation coefficents than
% observed (approximates the p-value). Note that signif = 1 if and only if F <= level. 
% 
% requires phaseran.m by Carlos Gias
% http://www.mathworks.nl/matlabcentral/fileexchange/32621-phase-randomization/content/phaseran.m

%  References:
%  ==========
% - Ebisuzaki, W, 1997: A method to estimate the statistical
% significance of a correlation when the data are serially correlated.
% J. of Climate, 10, 2147-2153.
%
% - Prichard, D., Theiler, J. Generating Surrogate Data for Time Series
% with Several Simultaneously Measured Variables (1994)
% Physical Review Letters, Vol 73, Number 7
%
% (Some Rights Reserved)  USC Climate Dynamics Lab, 2012.
% ====================================================


if nargin < 4; nsim = 1000; end
if nargin < 3; level = 0.05; end
if nargin < 2; error('Minimum input is the 2 time series vectors'); end

% rand('seed',-22459); % randn('seed',-22459); % Exact from Ebisuzaki's C code
% rand('seed',sum(100*clock)); randn('seed',sum(100*clock)); % better?
% rand('twister',sum(100*clock)); randn('twister',sum(100*clock)); % best?
rng('default')

x  = x(:); y=y(:);
n  = length(x);

if length(y) ~= n; error('Size of x and y must be the same'); end

r    = corr(x,y);
rSim = NaN(nsim,1);

% generate phase-randomized samples using the Theiler & Prichard method
X = sq(phaseran(x,nsim));
Y = sq(phaseran(y,nsim));

% compute correlations
Xs = standardize(X);
Ys = standardize(Y);
C = Xs'*Ys/(n-1); 
rSim = diag(C); 

% compute fraction of values higher than observed
F = sum(abs(rSim) >= abs(r))/nsim;

% establish significance
signif = logical(F < level); % significant or not?


return




