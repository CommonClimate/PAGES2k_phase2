function cvx = cfspl(p,yrv,yrvn,xvn)
% cfspl:  curve fit by cubic spline
% CALL: cvx = cfspl(p,yrv,yrvn,xvn);
%
% D Meko  -- pre-1997
%
%*************** IN 
%
% p -- spline parameter, as defined in spline toolbox
% yrv -- years at which smoothed values are needed
% yrvn -- years at which original time series has data;  some might
%	be NaN
% xvn -- original time series;  some might be NaN
%
%
%************** OUT
%
% cvx -- values of smooth curve at years yrv
%
%************ NOTES
%
% cfspl.m  calls MATLAB function csaps.m to compute spline, and does only 
% minor bookkeeping
%
% Calling program might use user-written splinep.m to convert from "period with
% 0.5 amplitude of frequency response" to spline parameter p as required by 
% spline toolbox

  % Curve fitting
   yrvn(isnan(xvn))=[];  % delete rows with xvn as NaN
   xvn(isnan(xvn))=[];
   
  
   cvx=csaps(yrvn,xvn,p,yrv);
   cvx=cvx';
  