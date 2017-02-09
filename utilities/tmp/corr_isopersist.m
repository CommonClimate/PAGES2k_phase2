function [r,signif,pval,g,rcrit] = corr_isopersist(x,y,alpha,nsim,iplot);
%   CORR_SIGNIF   
% 
%   [r,signif,pval,g,rcrit]=corr_isopersist(x,y,alpha,[iplot]);
%
%    Computes correlation between two timeseries, and their significance.
%  The latter is gauged via a non-parametric (Monte Carlo) simulation of
%  correlations with nsim AR(1) processes with identical persistence 
%  properties as x and y ; the measure of which is the lag-1 autocorrelation (g).
%  
%    Inputs :  - x, y ; vector timeseries 
%					- alpha ; level of the test (probability of a type I error) (default = 0.05)
%              - nsim = number of simulations (default = 1000);
%					- (optional) iplot = 1 activates plotting of the sampling distribution (default = 0)
%
%    Outputs : - r (Pearson's product moment correlation). 
%              - signif : significance (logical variable)
%					- pval : the p-value of the test *
%					- g : the lag-1 autocorrelation of x and y. (doublet)
%              - rcrit: the 95% quantile of the simulated |r|
%  
% * (the probability of obtaining a test statistic at least
%   as extreme as the one actually observed, assuming that 
%  the null hypothesis is true)
%  
%   The test is 1 tailed on |r|: Ho = { |r| = 0 }, Ha = { |r| > 0 }
% 
%   The test is rejected (signif = 1) if pval <= alpha, otherwise signif=0;
% 
%    (Some Rights Reserved) Hepta Technologies, 2009
%     v1.0 USC, Aug 10 2012, based on corr_signif.m
%
% ======================================================== 
% Argument check
error(nargchk(2, 5, nargin, 'struct'))
if nargin < 3, alpha = 0.05; end
if nargin < 4, nsim = 1000; end
if nargin < 5, iplot = 0; end

% ensure column form
x = x(:); y = y(:);
nx = length(x); ny = length(y);
if (nx ~= ny)
	error('x and y must be of equal length')
end		
n = nx; 
r = corr(x,y); % Correlation coefficient;
ra = abs(r);

% Generate Matrix of isopersistent AR(1) processes
[x_red,g(1)] = isopersistent_rn(x,nsim);
[y_red,g(2)] = isopersistent_rn(y,nsim);
for k = 1:nsim
   rs(k) = corr(x_red(:,k),y_red(:,k));
end
rsa = abs(rs);

% Find p-value
xi   = linspace(0,1.1*max([ra rsa]),200);
prob = ksdensity(rsa,xi);
[min_diff, pos] = min(abs(ra-xi));
pval = trapz(xi(pos:end),prob(pos:end)); % Pr( R > Rs)
%pval = sum(rsa >= ra)/n;
% Plot if asked
if iplot 
   fig('Correlation significance'),clf
   % Estimate the Empirical (Kaplan-Meier) CDF 
	[F,xf] = ecdf(rsa);
	[H,xh] = ecdfhist(F,xf,30);
	figure(1), clf
	bar(xh,H,'w'), hold on
	plot(xi,prob,'r-','linewidth',[2]);
	plot(ra,0.1,'kx','linewidth',[6])
	text(1.1*ra,0.5,'Observed')
end

% Percentile
rcrit  = prctile(rsa,100*(1-alpha));
signif = logical(ra >= rcrit); % significant or not?

return
end

%%    AUXILIARY FUNCTIONs %%%%%%%%%%%%%%%%%%%%%%

 function [red,g] = isopersistent_rn(X,M)
% function [red,g] = isopersistent_rn(X,M)
%    ISOPERSISTENT_RN :
%   Generates M realization of a red noise [ i.e. AR(1)]
% 		process with same persistence properties as X.
%   	(Mean and variance are also preserved).
%     Inputs : - Timeseries X (length = N).  
%					- Number of realizations M
%     Outputs :
%					- matrix of red nosie series (N by M)
%					- lag-1 autocorrelation g
%	
%		Rmk : preserves missing values if any
%
%   (Some Rights Reserved)  	Hepta Technologies, 2008
%   uses ar1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
warning off
%    
randn('state',sum(100*clock)); %set seed for random number generator;
 %
% extract dimensions 
N=length(X);
mu  = nmean(X);  % no NaNs !
sig = nstd(X,0,1); % 0 means \sqrt{N-1} normalization

%  get lag-1 autocorrelation coefficient 
%gam=autocorr(X);
%g1=gam(2); %second index is lag-1

g = ar1(X); % can fail for series too short or with large trends

if isnan(g)
	gam=autocorr(X);
	g=gam(2); %second index is lag-1
end	
red = red_noise(N,M,g);

% Center and normalize 
m = mean(red);  s = std(red,0,1);
red_n=(red-repmat(m,N,1))./repmat(s,N,1);
% Restore X's mean and std 
red=red_n.*repmat(sig,N,M)+ repmat(mu,N,M); 		

% Put NaNs where values were missing
red(isnan(X),:)=NaN; 	
	
return
end

%=======================================================
function [g,a,mu2]=ar1(x)
% AR1 - Allen and Smith AR(1) model estimation.
% Syntax: [g,a,mu2]=ar1(x);
%
% Input:  x - time series (univariate).
%
% Output: g - estimate of the lag-one autocorrelation.
%         a - estimate of the noise variance.
%         mu2 - estimated square on the mean.
%
% AR1 uses the algorithm described by Allen and Smith 1995, except
% that Matlab's 'fzero' is used rather than Newton-Raphson.
%
% Fzero in general can be rather picky - although
% I haven't had any problem with its implementation
% here, I recommend occasionally checking the output
% against the simple estimators in AR1NV.
%
% Alternative AR(1) estimatators: ar1cov, ar1nv, arburg, aryule
%
% Written by Eric Breitenberger.      Version 1/21/96
% Please send comments and suggestions to eric@gi.alaska.edu
%
% Updated,optimized&stabilized by Aslak Grinsted 2003-2005
%


x=x(:);
N=length(x);
m=mean(x);
x=x-m;

% Lag zero and one covariance estimates:
c0=x'*x/N;
c1=x(1:N-1)'*x(2:N)/(N-1);


g0=c1/c0; % Initial estimate for gamma

o=optimset('fzero');
o.TolX=0.0001;
o.Display='off';

% Find g by getting zero of 'gammest':
g=fzero(@gammest,g0,o); %-updated to accomodate newer versions of matlab
if isnan(g) %optimization failed use homemade zero finder... (aslak 2005)
    ginout=[g0*.999;g0];
    ginout(:,2)=[gammest(ginout(1,1)); gammest(ginout(2,1));];
	ii=0;
    while (ii<10)&(abs(ginout(end,2))>1e-6)
		ii=ii+1;
        dx=ginout(end,1)-ginout(end-1,1);
        dy=ginout(end,2)-ginout(end-1,2);
        xguess=ginout(end,1)-ginout(end,2)*dx/dy; %linear extrap
        ginout(end+1,:)=[xguess gammest(xguess)];
    end
    [m,mi]=nanmin(abs(ginout(:,2)));
    g=ginout(mi,1);
    if (m>1e-6)
        g=nan;
    end
end

gk=1:N-1;
gk=g.^gk;
mu2=(1/N)+(1/N^2)*2*sum((N-1:-1:1).*gk);

c0est=c0/(1-mu2);
a=sqrt((1-g^2)*c0est);


    function gout=gammest(gin)
        % GAMMEST - used by AR1 to compute
        % a function for minimization by fzero.
        %
        % Written by Eric Breitenberger.      Version 1/21/96
        % Please send comments and suggestions to eric@gi.alaska.edu
        %
        
        gk=1:N-1;
        gk=gin.^gk;
        mu2=(1/N)+(2/N^2)*sum((N-1:-1:1).*gk);
        gout=(1-g0)*mu2+g0-gin;
        if gout>1, gout=nan; end
    end
end

 function [red] = red_noise(N,M,g);
% function RED_NOISE
%		Produce AR1 process matrix with nv = [N M] and lag-1 autocorrelation g
%  
%      (Some Rights Reserved)  	Hepta Technologies, 2008
%  				J.E.G., GaTech, Oct 20th 2008
% ========================================================================

% Declare matrix
red = zeros(N,M);
% generate N rows by M realizations of red noise data
red(1,:)=randn(1,M);
for jt=2:N  % 
  red(jt,:)=g*red(jt-1,:)+randn(1,M);
end

return
end




