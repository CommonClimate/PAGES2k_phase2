function [RE,CE,R2,MSE]=verif_stats(Z,Zhat,calib,verif);
% function [RE,CE,R2,MSE]=verif_stats(Z,Zhat,calib,verif);
%
% VERIF_STATS : Computes a recontruction's verification statistics 
% on the specified calibration and validation intervals
%	
%   inputs :  Z = original time series
% 		  Zhat = estimated time series
%		calib,verif : array indices of the corresponding periods	
%
%	Note : written in matrix form in case Z is a collection of timeseries	
%			time  is assumed to be the first index
%		
% ref :  Cook, E.R,  K. R. Briffa, and P. D. Jones, 1994: Spatial regression 
% methods in dendroclimatology: A review and comparison of two 
% techniques. Int. J. Climatol., 14, 379-402.
%
%   created by Julien Emile-Geay, GaTech, Nov 15 2007
%   edited Dec 6 2013 to include the MSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Zv      = Z(verif,:); 
Zv_hat  = Zhat(verif,:);
[nv,p] = size(Zv);
% Mean Squared Error
MSE = mean((Zv-Zv_hat).^2);

%  Reduction of Error
Zbar = mean(Z(calib,:),1);
RE=1 - MSE./mean((Zv-repmat(Zbar,nv,1)).^2);
%  Coefficient of efficiency
Zbar = mean(Zv,1);
CE= 1- MSE./mean((Zv-repmat(Zbar,nv,1)).^2);

% correlation coefficient
R2=zeros(1,p);
for p=1:p
	rho=corrcoef(Zv(:,p),Zv_hat(:,p));
	R2(p)=rho(1,2)^2; % square off-diagonal element
end

return





