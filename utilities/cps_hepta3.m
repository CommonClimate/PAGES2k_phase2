function [rec,nest,nproxy] = cps_hepta3(proxy,target,tp,tt)
% ================================================
%   CPS : Composite Plus Scale 
%
%     [rec,diag] = cps_hepta3(proxy,target, calib)
%
%    inputs : 	- proxy matrix  (nm x nr) (missing values as NaNs)        
%             	- target series (nt x 1)  
%   			- tp: time axis for proxies (nm x 1) 
%               - tt: time axis for target  (nt x 1) 
% 
%    outputs:   - reconstructed series (nm x 1)
%  
%					
% 	 NB : time is understood to runs backwards (i,e. first row of the matrix is the most recent)
%  
%     Hepta Technologies, USC, June 2015. 
%   uses : ncorrcoef.m, nmean.m, nstd.m, verif_stat.m, overlap_frac.m, standardize.m
% ===========================================================================================

 nt     = length(target);
[nm,nr] = size(proxy);
% Compute first 2 moments over entire interval
mu = mean(target); sig = std(target,0,1);

% Compute correlations with target timeseries over calibration interval
for k = 1:nr
	R    = ncorrcoef(target,proxy(1:nt,k)); 
	rhot(k)  = sign(R(1,2));
end	
% find number of predictors
nproxy = zeros(nm,1);
for n=1:nm
	nz = ~isnan(proxy(n,:));
	nproxy(n) = sum(nz);
end

% 2)  STANDARDIZE AND COMPOSITE
%  ===========================================

% find common period
common = find(nproxy==max(nproxy));
pmean = nmean(proxy(common,:),1);
pstds = nstd(proxy(common,:),1);
% standardize using calibration mean and variance
proxys   = (proxy - repmat(pmean,nm,1))./repmat(pstds,nm,1); 


% 3)  NESTING
%  ===========================================

% Begin iterative nesting
jumps  = find(diff(nproxy)<0); % Track jumps in the  number of predictors
jumps  = [1 jumps'];
Nnest  = max(1,length(jumps)-1); % Number of jumps
rec  =  zeros(nm,1);
if Nnest > 1						
	for n=1:Nnest
		nest{n} = jumps(n):jumps(n+1); lnest=length(nest{n});
		[pnest,availn{n},missn{n}]=overlap_frac(proxys(nest{n},:),1:lnest,1);
		% select subset of proxy matrix
		pnest = proxys(1:jumps(n+1),availn{n}); lnestc = size(pnest,1);
		% Average proxys weighted by their correlation
		pavg{n} = nmean(pnest.*repmat(rhot(availn{n}),lnestc,1),2);
		% Compute mean and variance of composite over calibration period
		%[~,xm,xs] = standardize(pavg{n}(1:nt));
        [~,xm,xs] = standardize(pavg{n});

		% Recale the average
		pavgs{n} = ((pavg{n}-xm)/xs)*sig + mu ;
		% Splice them all together ~
		rec(nest{n})=pavgs{n}(nest{n});
	end
else % If constant availability of proxies
	% Average proxys weighted by their correlation
	pavg = nmean(proxys.*repmat(rhot,nm,1),2);
	% Compute mean and variance of composite over calibration period
	[dum,xm,xs] = standardize(pavg(1:nt));
	% Recale the average
	rec = ((pavg-xm)/xs)*sig + mu ;
end

end
