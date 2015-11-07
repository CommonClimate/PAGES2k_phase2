function [X,t,Xr,tr,mr,yr,npy,mRes,dt] = coral_interp(t_raw,X_raw,mcorr,mRes)
%   [X,t,Xr,tr,mr,yr,npy,dt] = coral_interp(t_raw,X_raw,mcorr,mRes)
%           interpolates a coral time series onto a regular time grid
%
% inputs  - t_raw: raw chronology in decimal years
%         - X_raw: raw data for dO18, Sr/Ca, etc.
%         - [mRes] : resolution (in months)   (default = average time spacing)
%         - [mcorr] : correction for datenum  (default = 1 month)
%
% outputs - X, t: raw data and chronology without duplicate values
%			 -	tr: even-interval time in decimal years
%         - Xr: data interpolated over tr
%		    - mr: evenly spaced months
%			 - yr: evenly spaced years
%         - npy: number of pts per year
%			 - dt: raw resolution in months
%			 - mRes: rounded resolution in months
%
%  History : created by Julien Emile-Geay, USC, October 2011
% =====================================================================

% make sure time must goes forward
if mean(diff(t_raw)) < 0
    X_raw = flipud(X_raw(:)); t_raw = flipud(t_raw(:));
end

% remove duplicate values
[t,X,ind] = consolidator(t_raw,X_raw);
dt = 12*abs(mode(diff(t))); % avg time spacing between consecutive values

% Make default correction factor 1 and default mRes the rounded-up resolution
if nargin < 3
    mcorr = 1; % correction factor to get proper months (i.e. from 1 to 12)
end
if nargin < 4
    if dt<=3 & dt>=1.5
        mRes = 2;
    elseif dt <1.5
        mRes = 1;
    elseif dt > 3
        mRes = ceil(dt); % conservative time resolution
    end
end



% correctly identify months and years
year = floor(t);  month = floor((t - year)*12+mcorr);
year_i = min(year); year_f = max(year);

ny = year_f - year_i +1; % number of years
first_month = month(1);
last_month = month(end);
dtm = 1/12; % regular monthly spacing
npy = 12/mRes; % number of samples per year in interpolated coral series

%  Interpolate on regular grid if needed
if std(diff(t))~=0 % if time intervals are not uniform
    % define regular time grid
    mr1 = [first_month:mRes:12]; n1 = numel(mr1);
    if mRes > 1
        mr2 = repmat([1+~isodd(first_month):mRes:12],[1 ny-2]); n2 = numel(mr2);
        mr3 = [1+~isodd(first_month):mRes:last_month]; n3 = numel(mr3);
    else
        mr2 = repmat([1:mRes:12],[1 ny-2]); n2 = numel(mr2);
        mr3 = [1:mRes:last_month]; n3 = numel(mr3);
    end
    mr = [ mr1 mr2 mr3]';  % concatenate 3 monthly grids
    yr1 = repmat(year_i,[1 n1]); % first yearly grid
    y = [year_i+1:year_f-1];
    for k = 1:ny-2  % middle yearly grid
        yr2(npy*(k-1)+1:npy*k) = y(k);
    end
    yr3 = repmat(year_f,[1 n3]); % last yearly grid
    yr = [ yr1 yr2 yr3]'; nr = length(yr);
    dr = repmat(15, [nr 1]); % assumes samples were taken on the 15th of every month
    %[doy,tr] = date2doy(yr,mr,dr);  % yields uneven time increments
    tr = yr+(mr+0.5)/12.0-dtm;	 % simple, but efficient
    % make sure extrapolation is not going on
    tr = tr(tr<=max(t) & tr>=min(t));
    Xr = interp1(t,X,tr,'pchip');  % interpolate onto regular grid
else
    tr = t;
    yr = year;
    mr = month;
    Xr = X;
end
dt = 12*abs(mode(diff(tr)));


return
end

