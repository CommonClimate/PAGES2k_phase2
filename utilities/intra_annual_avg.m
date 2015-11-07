function [da,ta,ts]=intra_annual_avg(d,t,ms,me,dir,thre);
% function [da,ta,ts]=intra_annual_avg(d,t,ms,me,dir,thre);
%
%  INTRA ANNUAL AVERAGE : averages climate data over some part 
%  of the annual cycle.
%
%  INPUT : 
%       d   : original data
%		t   : original timegrid
%		ms  : starting month
%		me  : end month
%		dir : direction of time grid : +1 means forward, -1 is backwards. (default=1)
%       frac: lowest acceptable fraction of monthly values [default = 0.75]  
%	OUTPUT :  
%		ta : annual time grid 
%		da : suitably averaged data
%       ts : datevec time grid
%
%  Assumes :  - fractional, EVENLY-SPACED timegrid (May 1978 is 1978.4166667)
%			  - months to be centered on the 15th day
%
%	===========================================================		
%  History : created Oct 16th 2007, J.E.G., GaTech.
%  Updated : 04-Dec-2014, J.E.G., USC.

if (nargin<5)
   dir =1; frac = 0.75;
elseif (nargin<6)
   frac = 0.75;
end	
nt=length(t); 
% fraction of available values


% convert time axis using datenum technology
year=floor(t);
month=floor((t-year)*12.0)+1;
ts=[year, month, repmat(15,nt,1) , repmat(0,nt,1), repmat(0,nt,1), repmat(0,nt,1)];
tn=datenum(ts);

if (dir==-1) % backwards   
	ta=[max(year):-1:min(year)]'; da=nan(length(ta),1);
	% do middle-chunk years
	for k = 1:length(ta)-1
		me_h=find(tn==datenum([ta(k),me,30,0,0,0]));
		ms_h=find(tn==datenum([ta(k),ms,1,0,0,0]));
        dm   = d(ismember(tn,me_h:ms_h));
        f = sum(~isnan(dm))/numel(dm);
        if f >= frac
    		da(k) = nmean(dm); 
        end
	end
else  % FORWARD
	ta=[min(year):max(year)]'; da= nan(length(ta),1);
	% do middle-chunk years
	for k = 1:length(ta)
		ms_h = datenum([ta(k),ms,1,0,0,0]);
		me_h = datenum([ta(k),me,30,0,0,0]);	
        dm   = d(ismember(tn,ms_h:me_h));
        f = sum(~isnan(dm))/numel(dm);
        if f >= frac
    		da(k) = nmean(dm); 
        end
    end
end



