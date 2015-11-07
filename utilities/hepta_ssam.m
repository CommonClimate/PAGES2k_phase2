function [spec,eig_vec,PC,RC,RCp,modes] = hepta_ssam(X,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function hepta_ssa : performs Singular Spectrum Analysis on time series X
%
% syntax : [spec,eig_vec,PC,RC,RCp,modes] = hepta_ssam(X,options)
% input: 
%    X, vector
%   options, structure with optional parameters [default value]:
%   - options.f : fraction (0<f<=1) of good data points for identifying
%   significant PCs [f = 0.3]
%   - options.M : window length [length(X)/10]
%   - options.crit_name : criterion used to determine significant modes. 
% Possible values are:
%     * 'kaiser' = retain all modes whose eigenvalue exceeds the median of
%     all eigenvalues (ref. 4, page 540) [this is the default]
%     * 'mcssa'  = Monte Carlo SSA (ref. 2)
%     * 'pct_var'= fraction of total variance retained (0<P<=1)
%
%   - options.K : number of significant modes to be used (overrides crit_name)
%   - options.P : minimum fraction of signal to reconstruct if 'pct_var' option is used [0.7]
%   - options.MC : number of MC simulations in 'mcssa' case 
%   - options.iplot : boolean variable to plot the MCSSA eigenvalues [0]
%
% output :
%   spec : eigenvalue spectrum, in % variance
%   eig_vec : eigenvector matrix ("temporal EOFs")
%   PC      : matrix of principal components
%   RC      : matrix of RCs (N*M, K) (only if K>0)
%   RCp     : Reconstructed timeseries
%   modes   : index of modes retained by the specified criteria
%
%   USC Climate Dynamics lab, 2013. (Maud Comboul & Julien Emile-Geay)
%
%   References:
%   [1] Vautard, R. and Ghil, M. (1989). Singular spectrum analysis in nonlinear
%       dynamics, with applications to paleoclimatic time series. Physica D, 35:395?424.
%   [2] Allen, M. R. and Smith, L. A. (1996). Monte Carlo SSA: Detecting irregular
%       oscillations in the presence of coloured noise. J. Clim., 9:3373?3404.
%   [3] Schoellhamer, D. H. (2001). Singular spectrum analysis for time series with missing data.
%      Geophysical Research Letters, 28(16):3187?3190.
%   [4] Wilks, D. S. (2011). Statistical Methods in the Atmospheric Sciences:
%         an Introduction. Academic Press, San Diego.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error(nargchk(1, 2, nargin))     % check number of input arguments

% first make sure the vector size is 1*N
X = X(:)';
N = length(X);

% set default values
f = 0.3;
M = round(N/10);
iplot = 0; 
kaiser_crit = true;
mcssa = false;
pct_var = false;

% if options parameter is not empty
if nargin > 1 && ~isempty(options)
   if isfield(options,'M')
      M = options.M;
      if M > N
         error('M should be less than the signal length')
      end
   end
   
   if isfield(options,'f')
      f = options.f;
      if f<=0 || f>1
         error('f must be in (0,1]')
      end
   end
   
   if isfield(options,'K')
      K = options.K;
      if K>M
         error('K should be less than the window length M')
      end
      modes = 1:K;
   end
   
   if isfield(options,'crit_name')
      if strcmpi(options.crit_name,'mcssa')
         mcssa = true;
         kaiser_crit = false;
         MC = 1000;
         if isfield(options,'MC')
            MC = options.MC;
         end
      elseif strcmpi(options.crit_name,'pct_var')
         pct_var = true;
         kaiser_crit = false;
         P = 0.7;
         if isfield(options,'P')
            P = options.P;
         end
      end
   end
end

% compute autocorelation
[Xr,mu,sigma] = standardize(X);

for j1=1:M
   j=j1-1;
   % modify for nan
   prod=Xr(1:N-j).*Xr(j+1:N);
   c(j1)=sum(prod(~isnan(prod)))/(sum(~isnan(prod))-1);
end

% Fill in Covariance matrix
C = toeplitz(c(1:M));

% solve eigenvalue problem
[eig_vec, eig_val] = eig(C);

deval      = diag(eig_val);
[deval, I] = sort(deval, 'descend'); % sort in decreasing order
eig_vec    = eig_vec(:, I);
eig_val    = eig_val(:,I);
nmodes     = length(deval); 
% determine significant eigenvalues
% =================================
if kaiser_crit
   modes = find(deval > median(deval))';
end

if pct_var
   sum_eval=sum(deval);
   var_percent=deval(1)/sum_eval;
   K=1;
   while var_percent < P*sum_eval
      var_percent = var_percent + deval(K);
      K=K+1;   
   end
   modes = 1:min(K,nmodes);
end

if mcssa
   [w,a,var,SBC,FPE,th] = arfit(Xr(~isnan(Xr))',1,1); %fit AR(1) model
   s = sqrt(var);
   noise = zeros(N,MC);
   noise(1,:) = repmat(Xr(find(~isnan(Xr),1)),[1 MC]);
   for jt=2:N  %
      noise(jt,:) = a(1)*noise(jt-1,:)+ s*randn(1,MC);
   end
   
   for m = 1:MC
      noise(:,m) = (noise(:,m)-mean(noise(:,m)))/std(noise(:,m));      
      [Gn,ln]       = xcorr(noise(:,m),M-1,'unbiased');
      Cn            = toeplitz(Gn(M:2*M-1));%/M;
      Lambda_R(:,m) =  diag(eig_vec * Cn * eig_vec'); % noise "eigenvalues"
   end
   q95 = quantile(Lambda_R,0.95,2);
   modes = find(deval > q95)'; % index of modes rising above the background
   if isempty(modes)
      display('modes is empty ');
   end
   if iplot
      display(['MCSSA modes retained: ', int2str(signif)]);
      fig('MCSSA'),clf
      v = [1:M]'; ligr = [ 0.7000    0.7000    0.7000];
      lmin = min(Lambda_R,[],2); lmax = max(Lambda_R,[],2);
      area_fill(v',lmin',lmax',ligr,ligr,0,0.3),hold on
      plot(v,eig_val,'kx',v,q95,'r-','linewidth',[2])
   end
end

%
% determine principal component time series (eqn b.6) in [3]
%
PC=nan*ones(N-M+1,M);
for k=1:M
   for i=0:N-M
      %   modify for nan
      prod=Xr(i+1:i+M).*eig_vec(:,k);
      ngood=sum(~isnan(prod));
      %   must have at least m*f good points
      if ngood>=M*f
         PC(i+1,k)=sum(prod(~isnan(prod)))*M/ngood; % the columns of this matrix are Ak(t), k=1 to M (T-PCs)
      end
   end
end

% compute reconstructed timeseries if K > 0
Np=N-M+1;

if ~isempty(modes)
   RC=zeros(N,length(modes));
   % first M terms
   for t=1:M-1
      Av=flipud(PC(1:t,modes));
      eig_vec_red=eig_vec(1:t,modes);
      RC(t,:)=1/t*sum(Av.*eig_vec_red,1);
   end
   %middle of timeseries
   for t=M:Np
      Av=flipud(PC(t-M+1:t,modes));
      eig_vec_red=eig_vec(1:M,modes);
      RC(t,:)=1/M*sum(Av.*eig_vec_red,1);
   end
   % last M terms
   for t=Np+1:N
      Av=flipud(PC(t-M+1:Np,modes));
      eig_vec_red=eig_vec(t-N+M:M,modes);
      RC(t,:)=1/(N-t+1)*sum(Av.*eig_vec_red,1);
   end
   % sum modes and restore the mean and variance
   RCp = sigma*sum(RC,2)+ mu;
else
   RC = []; RCp = [];
end

spec = deval;
