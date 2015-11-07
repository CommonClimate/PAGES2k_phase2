function [da,ta,ts]=intra_annual_avg(d,t,ms,me,dir);
% function [da,ta,ts]=intra_annual_avg(d,t,ms,me,dir);
%
%  INTRA ANNUAL AVERAGE : averages climate data over some part 
%  of the annual cycle.
%
%  INPUT : 
%  	d: original data
%		t: original timegrid
%		ms: starting month
%		me : end month
%		dir= direction of time grid : +1 means forward, -1 is backwards. (default=1)
%
%	OUTPUT :  
%		ta : annual time grid 
%		da : suitably averaged data
%     ts : datevec time grid
%
%  Assumes :  - fractional timegrid (May 1978 is 1978.4166667)
%				  - months to be centered on the 15th day
%
%   IMPORTANT : end points are NaN by default : you have to do them manually as they don't 
%     always contain enough months to produce a meaningful 'seasonal' average
%	===========================================================		
%  History : created Oct 16th 2007, J.E.G., GaTech.

if (nargin<5)
   dir=1;
end	


nt=length(t);
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
		ms_h=find(tn==datenum([ta(k+1),ms,1,0,0,0]));
		da(k)=nmean(d(me_h:ms_h)); 
	end
else  % FORWARD
	ta=[min(year):max(year)]'; da= nan(length(ta),1);
	% do middle-chunk years
	for k = 1:length(ta)-1
		ms_h = datenum([ta(k),ms,1,0,0,0]);
		me_h = datenum([ta(k+1),me,30,0,0,0]);		
		da(k)=nmean(d(ismember(tn,ms_h:me_h))); 
    end
    ta = ta + 1;
end


