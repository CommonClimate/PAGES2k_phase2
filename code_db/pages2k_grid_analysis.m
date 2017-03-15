function pages2k_grid_analysis(vers,options)
% Execution time is about 300 seconds on a 2014-vintage Mac Pro

addpath(genpath('../../utilities'))

if options.norm_p
    norm_string = 'normal';
else
    norm_string = 'raw';
end

% options
% method         = 'isospec'
% sample_thresh  = options.sample_thresh;

year = [];
fn = ['../data/PAGES2k_v' vers '_unpack.mat'];
load(fn)

temp_direc =cellfun(@lower,{S.climateInterpretation_interpDirection},'UniformOutput',false);

directionMultiplier = NaN(size(temp_direc));
for d = 1:size(temp_direc,2)
    if strcmp(temp_direc(d),'negative')
        directionMultiplier(d) = -1;
    elseif strcmp(temp_direc(d),'positive')
        directionMultiplier(d) = 1;
    end
end


%% Prepare data
tmp  = load('../data/had4med_graphem_sp70.mat');  % TO DO: make this more general: should be able to load any temperature dataset, if stored under a more generic name
ti   = unique(tmp.tvec(1:1968,1));  % TO DO: remove ad-hockeries
temp = tmp.Xf(1:1968,:); % 164 full calendar years from Jan 1850 - Dec 2013

t     = year(year>0);

% Overlapping period (proxy vs. temperature)
tc    = intersect(ti,t);
tcal  = ismember(ti,tc);
pcal  = ismember(t,tc);

% Annual
annual = (temp(1:12:end,:) + temp(2:12:end,:) + temp(3:12:end,:) + ...
    temp(4:12:end,:) + temp(5:12:end,:) + temp(6:12:end,:) + ...
    temp(7:12:end,:) + temp(8:12:end,:) + temp(9:12:end,:) + ...
    temp(10:12:end,:) + temp(11:12:end,:) + temp(12:12:end,:))./12;

% JJA
summer = (temp(6:12:end,:) + temp(7:12:end,:) + temp(8:12:end,:))./3;

% DJF
winter = (temp(12:12:end,:) + temp(1:12:end,:) + temp(2:12:end,:))./3;

% Restrict temperature and proxy to the overlapping period
annual = annual(tcal,:); % mat = annual;
summer = summer(tcal,:); % jja = summer;
winter = winter(tcal,:); % djf = winter;

% Remove duplicate records
if options.norm_p
    disp('Normalize proxy before calculating correlations')
    norm_string = 'normal';
    proxc_ann = gaussianize(proxy_ann(pcal,:));
    proxc_djf = gaussianize(proxy_djf(pcal,:));
    proxc_jja = gaussianize(proxy_jja(pcal,:));
else
    norm_string = 'raw';
    proxc_ann = proxy_ann(pcal,:);
    proxc_djf = proxy_djf(pcal,:);
    proxc_jja = proxy_jja(pcal,:);
end

avail = ~isnan(proxc_ann);
ns    = size(temp,2);
np    = size(proxc_ann,2);



%% Decide if any pre-processing of data is needed
per = {'ann','djf','jja','ama'};
if strcmp(options.dtype,'detrend')
    annual = detrend(annual);
    summer = detrend(summer);
    winter = detrend(winter);
    for v = 1:np
        for s = 1:3
            eval(['Xa = proxc_' per{s} '(:,v)']);
            ta = tc;
            if sum(avail(:,v))>10
                [~,~,segment] = find_breaks(ta,Xa);
                nseg = length(segment);
                dtd  = NaN*Xa;
                if nseg > 1  % if more than 1 segment
                    for i=1:nseg;
                        ds = Xa(segment{i});
                        dtd(segment{i}) = detrend(ds);
                    end
                else
                    ind = ~isnan(Xa);
                    ds = Xa(ind);
                    % Maker sure they are column vectors
                    dtd(ind) = detrend(ds);
                end
                eval(['proxc_' per{s} '(ismember(tc,ta),v) = dtd;']);
            end
        end
    end
elseif strcmp(options.dtype,'diff1')
    % or first-difference
    annual = diff(annual);  summer = diff(summer); winter = diff(winter);
    proxc_ann  = diff(proxc_ann);
    proxc_djf  = diff(proxc_djf);
    proxc_jja  = diff(proxc_jja);
    avail  = ~isnan(proxc_ann);
end

locs  = tmp.loc;
lon = p_lon;
lat = p_lat;


%% regrid for plotting
ulat = unique(locs(:,2));
ulon = unique(locs(:,1));
sites = [lon' lat'];

ftarget = annual;
ptarget = proxc_ann .* repmat(directionMultiplier,size(ftarget,1),1);
cuttoffValue = 0.05;

for i = 1:length(sites);

    ilat = findnearest(ulat,lat(i));
    ilon = findnearest(ulon,lon(i));

    pointer = find(locs(:,2)==ulat(ilat(1)) & locs(:,1)==ulon(ilon(1)));

    [R,P] = corrcoef([squeeze(annual(:,pointer)) ptarget(:,i)],'rows','pairwise');
    r.annual(i) = R(2,1); p.annual(i) = P(2,1);

    good = ~isnan(squeeze(annual(:,pointer))) & ~isnan(ptarget(:,i));
    ncorrs(i) = sum(good);

    if ncorrs(i)>1
        [r.isospec(i),signif.isospec(i),F] = corr_isospec(squeeze(annual(good,pointer)),ptarget(good,i));
    else
        r.isospec(i) = NaN;
        signif.isospec(i) = 0;
    end
end

localCorrelation = r.isospec;
localSignificance = signif.isospec;
c = pi/180;

for i = 1:ns
    %thisGrid = locs(i,:); %
    %siteDistances = earthDistances([thisGrid; sites]); % long, lat
    %pointDistance = siteDistances(2:end,1);
    pointDistance = greatCircleDistance(c*locs(i,2), c*locs(i,1), c*p_lat, c*p_lon);
    inRadius = find(pointDistance<=options.search_radius);
    nsites(i) = length(inRadius);

    if nsites(i) > 0
        [R,P] = corrcoef([squeeze(ftarget(:,i)) ptarget(:,inRadius)],'rows','pairwise');
        R0 = R(2:end,1); % correlation between proxies and gridpoint
        P0 = P(2:end,1);

        significant = find(P0<=cuttoffValue);

        gridMeanAll(i)   = nanmean(R0);
        gridMedianAll(i) = nanmedian(R0);
        gridMaxAll(i)    = nanmax(R0);

        if isempty(significant)
            gridMean(i)   = NaN;
            gridMedian(i) = NaN;
            gridMax(i)    = NaN;
        else
            gridMean(i)   = nanmean(R0(significant));
            gridMedian(i) = nanmedian(R0(significant));
            gridMax(i)    = nanmax(R0(significant));
        end
    else
        gridMeanAll(i)   = NaN;
        gridMedianAll(i) = NaN;
        gridMaxAll(i)    = NaN;

        gridMean(i)   = NaN;
        gridMedian(i) = NaN;
        gridMax(i)    = NaN;
    end
end

ggridAll = NaN(4,length(ulon),length(ulat));
ggrid    = NaN(4,length(ulon),length(ulat));

for ii = 1:length(ulon)
    for jj = 1:length(ulat)
        pointer = find(locs(:,1)==ulon(ii) & locs(:,2)==ulat(jj))
        ggridAll(1,ii,jj) = gridMeanAll(pointer);
        ggridAll(2,ii,jj) = gridMedianAll(pointer);
        ggridAll(3,ii,jj) = gridMaxAll(pointer);
        ggridAll(4,ii,jj) = nsites(pointer);

        ggrid(1,ii,jj) = gridMean(pointer);
        ggrid(2,ii,jj) = gridMedian(pointer);
        ggrid(3,ii,jj) = gridMax(pointer);
        ggrid(4,ii,jj) = nsites(pointer);
    end
end

fn = ['../data/pages2k_hadcrut4_gridcorr_' options.dtype '_' vers '_' norm_string '_annual.mat'];
save(fn)
