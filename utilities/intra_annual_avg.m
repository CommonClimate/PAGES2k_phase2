function [da,ta,ts]=intra_annual_avg(d,t,ms,me,frac);
% function [da,ta,ts]=intra_annual_avg(d,t,ms,me,frac);
%
%  INTRA ANNUAL AVERAGE : averages climate data over some part
%  of the annual cycle.
%
%  INPUT :
%       d   : original data
%		t   : original timegrid
%		ms  : starting month
%		me  : end month
%       frac : lowest acceptable fraction of monthly values [default = 0.75]
%
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
%               fractional averaging when the period is incomplete
%            02-Feb-2016, J.E.G., USC.  
%               fixed the year-straddling bug and refactored


t = t(:); %ensure column vector is used
dir = sign(median(diff(t))); % direction of time flow

if (nargin<5)
    frac = 0.66;
end
nt=length(t);

% convert time axis using datenum technology
year=floor(t);
month=floor((t-year)*12.0)+1;
ts=[year, month, repmat(15,nt,1) , repmat(0,nt,1), repmat(0,nt,1), repmat(0,nt,1)];
tn=datenum(ts);

if (dir==-1) % backwards
    ta=[max(year):-1:min(year)]'; da=nan(length(ta),1);
    for k = 1:length(ta)
        ms_h = datenum([ta(k),ms,1,0,0,0]);
        if me>ms  %  end month is after start month: stay in curren year.
            me_h = datenum([ta(k),me,30,0,0,0]);
        else  % jump to next year
            me_h = datenum([ta(k)-1,me,30,0,0,0]);
        end
        dm  = d(ismember(tn,me_h:ms_h));
        f   = sum(~isnan(dm))/numel(dm);
        if f >= frac
            da(k) = nmean(dm);
        end
    end
    % update year axis if the majority of the months are in year -1
    if me<ms & numel([1:me])<numel([ms:12]) 
        ta = ta-1;
    end
else  % FORWARD
    ta=[min(year):max(year)]'; da=nan(length(ta),1);
    for k = 1:length(ta)
        ms_h = datenum([ta(k),ms,1,0,0,0]);
        if me>ms  %  end month is after start month: stay in curren year.
            me_h = datenum([ta(k),me,30,0,0,0]);
        else  % jump to next year
            me_h = datenum([ta(k)+1,me,30,0,0,0]);
        end
        dm  = d(ismember(tn,ms_h:me_h));
        f   = sum(~isnan(dm))/numel(dm);
        if f >= frac
            da(k) = nmean(dm);
        end
    end
    % update year axis if the majority of the months are in year+1
    if me<ms & numel([1:me])>numel([ms:12]) 
        ta = ta+1;
    end
end










