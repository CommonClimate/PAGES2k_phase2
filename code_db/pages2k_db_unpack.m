function  S = pages2k_db_unpack(vers,visual)
% function pages2k__db_unpack_flat(vers)
%  INPUT:   vers (version of the database, string)
%  OUTPUT:  saved to Pages2kPhase2Database_{vers}_unpack
%
%  Loads data from one of the pages2k_vX.Y.Z.mat files and extract
%  data and metadata for use in climate applications.
%
%  started Nov 23 2015 by Julien Emile-Geay (USC), based on code by by J. Wang
%  Updated Feb  6 2017 for public release (JEG)

% Prepare T-sensitive records
% load data and packages
addpath(genpath('../utilities'))
load(['../data/PAGES2k_v' vers '.mat'])

% Keep a record if it is interpreted as temperature AND
clim_interp = {TS.climateInterpretation_variable};
keep        = find((strcmp(clim_interp,'T')) & strcmp({TS.paleoData_useInGlobalTemperatureAnalysis}, 'TRUE'));
S     = TS(keep);   nr = length(S);
not_kept    = setdiff([1:nr],keep);  % (diagnostic purposes only)

% Archive type
archive = {S.archiveType};
archive = strrep(archive,'tree ring','tree'); %fix Laanila archive type 


yearMin = nan(nr,1);
yearMax = nan(nr,1);
resMed = nan(nr,1);
resAvg = nan(nr,1);
resMax = nan(nr,1);
resMin = nan(nr,1);

% Extract time infomation first
for r = 1:nr
    t = S(r).year;
    [yearMax(r),yearMin(r)] = maxmin(floor(t));
end
year_i = min(yearMin); year_f = max(yearMax);
year   = year_i:year_f; ny = length(year);

% Aggregate proxy timseries into 3 timeseries, f(seasonal sensitivity)
proxy_ann = nan(ny,nr);
proxy_djf = nan(ny,nr);
proxy_jja = nan(ny,nr);
proxy_ama = nan(ny,nr);

% parameters for SSA interpolation
op.crit_name = 'pct_var'; op.P = 1.00;
n = 0; m = 0; gappy = zeros(nr,1);
for r = 1:nr
    %for r = 2  %debug
    T = S(r);
    disp(['Importing record ' int2str(r) '/' int2str(nr)])

    % Exclude NaNs in chronology. Just in case
    t     = T.year;     noNaN = ~isnan(t);
    t     = t(noNaN);   dt    = median(diff(t));
    X     = T.paleoData_values(noNaN);

    %  Some records go backwards. If so, flip them so everyone is forward-looking
    if dt < 0
        t  = flipud(t(:));
        dt = -dt;
        X = flipud(X);
    end

    % check for missing values and refine yearMin, YearMax
    noNaNx = ~isnan(X);
    [yearMax(r),yearMin(r)] = maxmin(floor(t(noNaNx)));

    if sum(noNaNx)< 2
        % display(['record ',num2str(r),', site ' T.dataSetName])
        display(['record ',num2str(keep(r)),', site ' T.dataSetName])
        n = n+1;
        miss(n) = r;
    else
        % CLEAN UP AND ANNUALIZE RECORDS
        tn = t(noNaNx); Xn = X(noNaNx);
        resMed(r) = median(diff(tn));
        resAvg(r) = mean(diff(tn));
        resMax(r) = max(diff(tn));
        resMin(r) = min(diff(tn));

        % remove duplicates
        [tc,Xc,ind] = consolidator(t(:),X(:));
        % find breaks
        [~,~,segment] = find_breaks(tc,Xc);
        nseg = length(segment);
        % find first and last data points
        tmin = min(tc(segment{1})); tmax = max(tc(segment{nseg}));
        interval = find(tc >= tmin & tc <= tmax);
        ti = tc(interval); Xi = Xc(interval);
        % find breaks
        [~,~,segment] = find_breaks(ti,Xi);

        gfrac(r) = numel(Xn)/numel(X); % gap fraction

        % identify gappy records
        if nseg > 1 & resMax(r) > 10 & resMed(r) < 1 % if more than 1 segment, and large gaps
            gappy(r) = 1;
        end

        if resMed(r) < 1 % if subannual, coarse-grain to annual resolution
            if gappy(r)  % there are large gaps, proceed segment by segment
                ta = unique(floor(ti))+1;
                Xcal = NaN(size(ta)); Xdjf = Xcal; Xjja = Xcal;
                for s = 1:nseg
                    ts = tc(segment{s}); Xs = Xc(segment{s});
                    % interpolate onto regular time grid
                    [~,~,Xr,tr] = subannual_interp(ts,Xs,dt);
                    % infill using singular spectrum analysis
                    [spec,eig_vec,PC,RC,Xi,modes] = hepta_ssam(Xr,op);
                    % intra annual averages
                    [cal,ty,~]  = intra_annual_avg(Xi,tr,1,12,1); % calendar year
                    Xcal(ismember(ta,ty)) = cal;
                    [ama,tm,~]  = intra_annual_avg(Xi,tr,4,3,1);  % tropical year
                    Xama(ismember(ta,tm)) = ama;
                    djf         = intra_annual_avg(Xi,tr,12,2,1); % DJF
                    Xdjf(ismember(ta,tm)) = djf;
                    jja         = intra_annual_avg(Xi,tr, 6,8,1);  % JJA
                    Xjja(ismember(ta,ty)) = jja;
                end
            else
                % interpolate onto regular time grid
                mRes = round(resMed(r)*12);
                [~,~,Xr,tr] = subannual_interp(tc,Xc,1,mRes);
                % infill using singular spectrum analysis
                if sum(isnan(Xr))>0
                    [spec,eig_vec,PC,RC,Xi,modes] = hepta_ssam(Xr,op);
                else
                    Xi = Xr;
                end
                % intra annual averages
                [Xcal,ta,~]  = intra_annual_avg(Xi,tr,1,12);
                [Xama,tm,~]  = intra_annual_avg(Xi,tr,4,3);
                Xdjf         = intra_annual_avg(Xi,tr,12,2);
                Xjja         = intra_annual_avg(Xi,tr, 6,8);

                if sum(~ismember(ta,unique(floor(tc))))~=0
                    % make sure corals are not annualized where there were original NaNs
                    Xcal(~ismember(ta,unique(floor(tc)))) = nan;
                    Xdjf(~ismember(ta,unique(floor(tc)))) = nan;
                    Xjja(~ismember(ta,unique(floor(tc)))) = nan;
                    Xama(~ismember(ta,unique(floor(tc)))) = nan;
                end
            end

            %sanita check
            if visual
                fig('Annualization'), clf
                plot(tc,Xc,'marker','o','color',rgb('Gray'),'MarkerSize',8,'linestyle','none'); hold on;
                plot(tm,Xama,'k-','linewidth',2);
                if resMed(r) < 1
                    plot(tm,Xdjf,'b-',ta,Xjja,'r-');
                end
                hold off; lab{1} = 'Raw'; lab{2} = 'Annualized'; legend(lab{:}), legend boxoff
                ylab = [T.paleoData_variableName ' (' removeLeadingAndTrailingSpaces(T.paleoData_units) ')'];
                site_n = strrep(T.dataSetName,'_','\_');
                ttl = [int2str(r),') ', archive{r},', ', site_n];
                fancyplot_deco(ttl,'Year (CE)',ylab,14);
                filen=['./figs/annualize/pages_2k_phase2_record_', sprintf('%03d',r) '.pdf'];
                export_fig(filen,'-r100','-cmyk','-painters','-nocrop')
            end
        elseif resMed(r)>=1 && resMin(r) < 1
            ta = [min(floor(tc)):max(floor(tc))]';
            F = griddedInterpolant(tc,Xc,'linear'); Xcal = F(ta);
            Xdjf = Xcal; Xjja = Xcal;  Xama = Xcal;
        else
            [tc,Xc,ind] = consolidator(tn(:),Xn(:));
            ta = round(tc);
            Xcal = Xc; Xdjf = Xc; Xjja = Xc;  Xama = Xc;
        end
        if strcmp(archive{r},'ice core')
            proxy_ann(ismember(year,ta),r) = Xcal;  % As per Eric Steig: diffusion renders subannual averages meaningless
            proxy_djf(ismember(year,ta),r) = Xcal;
            proxy_jja(ismember(year,ta),r) = Xcal;
            proxy_ama(ismember(year,ta),r) = Xcal;
        else
            proxy_ann(ismember(year,ta),r) = Xcal;
            proxy_djf(ismember(year,ta),r) = Xdjf;
            proxy_jja(ismember(year,ta),r) = Xjja;
            proxy_ama(ismember(year,ta),r) = Xama;
        end
        clear Xcal Xdjf Xjja ta tm
    end
end

% special case of PacificGBRWei2009: resolution too darn variable to handle
% Wei, G., M.T. McCulloch, G. Mortimer, W. Deng, and L. Xie. 2009. Evidence for ocean acidification in the Great Barrier Reef of Australia. Geochimica et Cosmochimica Acta, vol. 73, pp. 2332-2346. doi:10.1016/j.gca.2009.02.009
names = {S.dataSetName};
iWei = find(strcmpi(names,'Ocean2kHR-PacificGBRWei2009')); %2 hits
for r = iWei
    T = S(r);
    tc = T.year;
    Xc = T.paleoData_values;
    resMed(r) = 5;
    resAvg(r) = 5;
    resMax(r) = 5;
    resMin(r) = 5;
    t5 = [min(tc):5:max(tc)]; % 5y
    Xm = nmean(Xc);
    X5 = pchip(tc,Xc-Xm,t5)+Xm;
    ty = [min(tc):1:max(tc)]; % yearly
    Xcal = pchip(t5,X5,ty);

    % save to proxy matrices
    proxy_ann(ismember(year,ty),r) = Xcal;
    proxy_djf(ismember(year,ty),r) = Xcal;
    proxy_jja(ismember(year,ty),r) = Xcal;
end
% restrict to Common Era
yearCE    = year>0;
proxy_ann = proxy_ann(yearCE,:);
proxy_djf = proxy_djf(yearCE,:);
proxy_jja = proxy_jja(yearCE,:);
year      = year(yearCE);

% Coordinates
p_lat = [S.geo_meanLat];
p_lon = [S.geo_meanLon];

clear TS D tr Xr tc Xc ts Xs spec eig_vec PC RC  Xi modes Xn tn n m r Xcal
clear gappy gfrac interval iWei noNan noNaNx nseg ny segment t T ta ti tmax
clear tmin ty X X5 Xm clim_interp F me ms

recordNames = {S.dataSetName};
vname       = {S.paleoData_variableName};
units       = {S.paleoData_units};


save(['../data/PAGES2k_v' vers '_unpack.mat'],'-v6')

%save with v6 because of this:
%http://undocumentedmatlab.com/blog/improving-save-performance
% and because -v7.3 is far too inefficient on Mac OS Sierra
