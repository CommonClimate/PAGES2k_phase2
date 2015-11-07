function [DT,T] = spl_dt(time, X, per,varargin)
%
%Usage:
%[DT,T] = spl_dt(time,X,per,rmvar);
%
%Inputs:
%time........Vector of time observations, may be subannual (i.e. 1984.833)
%X...........Data vector whos columns are to be detrended
%per.........Period, in years, above which 50% of the variance should be
%            removed
%
%rmvar.......(optional) percentage of variance to be removed at the period 
%            specified by "per". This value must be greater than zero, and 
%            less than one. If left blank, the default is 0.5 (50%).

%Outputs:
%DT..........Detrended version of X
%T...........Spline that was subtracted from each row of X
%
% 
%UW functions called:
% splinep (D. Meko)
% cfspl (D. Meko)
%
% Adapted by T Ault

[N,M] = size(X);
DT = nan(N,M);
T = nan(N,M);

if isempty(varargin);
    rmvar=0.5;
else
    rmvar=varargin{1};
end

for i = 1:M
    clear q x xdt
    q = find(~isnan(X(:,i)));
    x = X(q,i);
    yr = time(q);
    if ~isempty(q)
        p1 = splinep(per,rmvar);  % spline parameter for chosen spline
        cv1 = (cfspl(p1,yr,yr,x));
        xdt = x-cv1';
        DT(q,i) = xdt;
        T(q,i) = cv1;
    else
        DT(q,i) = NaN;
        T(q,i) = NaN;
    end
end

function p = splinep(per,amp)
% splinep:  spline smoothing parameter from Cook's period of spline response
% CALL: p = splinep(per,amp);
%
% Meko 3-4-92
%
%************  IN 
%
% per (1 X 1)r - period in years at which amplitude of freq response should be amp
% amp (1 x 1)r - desired amplitude of freq response at period per years
%
%************ OUT 
%
% p (1 x 1)r spline smoothing parameter, following Cook and Peters (1981) 
%  corresponding to a spline with the desired amplitude of frequency response
%  at period per years.  This p is same as the "spline smoothing parameter" P
%  required as input to the MATLAB spline toolbox function csaps.m.
%___________________________________________________________________

% Convert period per in frequency f

f=1.0/per;

% Use eq. (2) of Cook and Peters to find p

p=6*amp*(cos(2*pi*f)-1)^2/((cos(2*pi*f)+2)*(1-amp));


% End of function

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


