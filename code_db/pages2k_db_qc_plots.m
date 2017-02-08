function pages2k__db_qc_plots(vers,options)
% function pages2k__db_qc_plots(vers)
%  INPUT:   vers: (version of the database, string)
%           dtype: (pre-processing of timeseries, options are: detrending,
%               first-difference of the original time series, and the original)
%           n_neigh: %number of proxy neighbors to compare to (if no temperature calibration possible)
%           norm_p : boolean flag for inverse transform sampling to
%           normality
%  OUTPUT:  a whole of lot figures (if export = 1)
%
% Visual check of proxy records (especially its sensitivity to
% temeperature or similarity to proxy neighbors), and metadata.
fontname = 'Helvetica';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

dtype = options.dtype;

year = [];
addpath(genpath('../utilities'))
load JEG_graphics % load some graphical definitions
load(['../data/PAGES2k_v' vers '_unpack']) % load proxy database
load ../../data/distances.mat   % load P-T and P-P distances
C = m_shaperead('../../data/NaturalEarth/ne_50m_admin_0_countries'); % load country boundaries
P = m_shaperead('../../data/NaturalEarth/ne_50m_admin_1_states_provinces_lakes');  % load province boundaries for large countries

% load correlation file
method = options.method;
%load ../../data/corr_hadcrut/pages2k_hadcrut4_corr_original_1-1020
if options.norm_p
    fname = ['../../data/corr_hadcrut/pages2k_hadcrut4_corr_' dtype '_' vers '_normal_' method];
else
    fname = ['../../data/corr_hadcrut/pages2k_hadcrut4_corr_' dtype '_' vers '_raw_' method];
end
load(fname)

%% Coordinates and other variables in plot
x  = -180:5:180; nx = length(x);
y  = -90:5:90;   ny = length(y);

tmp  = load('../../data/temp/had4med_graphem_sp70');
locs = tmp.loc;
inds = (locs+2.5)./5;
ilat = inds(:,2)+18;
ilon = inds(:,1)+36;

ns    = size(locs,1);
np    = size(rho_mat,2);
nt    = length(year);

%% Data
% instrumental period
ti    = unique(tmp.tvec(1:1968,1)); % only take complete years
ti_min  = min(ti); % start of instrumental period
% overlap with temperature
pcal  = ismember(year,ti);

% proxy over the calibration period
proxc = proxy_ann(pcal,:);
avail = ~isnan(proxc);

% plotting parameters
xrange = [-0.6 0.6];  % correlation range
n      = options.n_neigh; % number of proxy neighbors
nc     = 21; % number of contours
samples_thresh = 20; %minimum number of samples for a correlation to be meaningful
res    = abs(resMed);
warning off

% clear all pdfs from subdirectories
regions =  unique({S.geo_pages2kRegion}); nreg = numel(regions);
if options.erase
    for k = 1:nreg
        region = regions{k};
        unix(['rm -rf ../../figs/qc/' region '/*.pdf']);
        unix(['rm -rf ../../figs/qc/' region '/*.eps']);
    end
end

% sort records by pages2k id.
ID=({S.paleoData_pages2kID})';
[IDs isort] = sort(ID);

%% Plot the living daylight out of all records

instCalib = false(np,1);
rho_n     = cell(np,1);
signif_n  = cell(np,1);

for s = 1:nr
    r = isort(s);
    T = S(r);
    %% extract raw data
    to = T.year;
    Xo = T.paleoData_values;
    if median(diff(to))<0
        to = flipud(to);
        Xo = flipud(Xo);
    end
    %% extract annualized data
    Xa    = proxy_a(:,r);
    ta    = year(~isnan(Xa));

    n_raw_samples(r) = sum(to >= min(ti_min));
    n_ann_samples(r) = sum(ta >= min(ti_min));

    %% metadata strings
    if ~iscell(T.paleoData_values) && isvector(T.paleoData_values)...
            && isempty(strfind(lower(T.paleoData_variableName),'depth')) ...
            && isempty(strfind(lower(T.paleoData_variableName),'distance')) ...
            && isempty(strfind(lower(T.paleoData_variableName),'year'))...
            && sum(~isnan(T.paleoData_values(:)))~=0
        ri=strfind(T.paleoData_variableName,'/');
        if ~isempty(ri)
            T.paleoData_variableName=T.paleoData_variableName(setdiff(1:length(T.paleoData_variableName),ri));
        end
        ri=strfind(T.paleoData_variableName,'%');
        if ~isempty(ri)
            T.paleoData_variableName=T.paleoData_variableName(setdiff(1:length(T.paleoData_variableName),ri));
        end

        fn=fieldnames(T);
        for i=1:length(fn)
            if isempty(T.(fn{i}))
                T.(fn{i})='?';
            end
        end
        if iscell(T.climateInterpretation_seasonality)
            seasonality = T.climateInterpretation_seasonality{1};
        elseif isnumeric(T.climateInterpretation_seasonality)
            seasonality = num2str(T.climateInterpretation_seasonality);
        else
            seasonality = T.climateInterpretation_seasonality;
        end

        if iscell(T.pub1_pubYear)
            pub1_year = num2str(T.pub1_pubYear{1});
        elseif isnumeric(T.pub1_pubYear)
            if length(T.pub1_pubYear) > 1
                pub1_year = num2str(T.pub1_pubYear(1));
            else
                pub1_year = num2str(T.pub1_pubYear);
            end
        elseif ischar(T.pub1_pubYear)
            pub1_year = T.pub1_pubYear;
        else
            pub1_year = '?';
        end

        % Only show the first DOI
        if iscell(T.pub1_DOI)
            if strfind(T.pub1_DOI{1},'http')
                pubDOI = T.pub1_DOI{2};
            else
                pubDOI = T.pub1_DOI{1};
            end
        elseif ischar(T.pub1_DOI)
            % if ~isempty(strfind(T.pub1_DOI,''''))
            pubDOI = char(regexp(T.pub1_DOI,'10.*','match'));
            if isempty(pubDOI)
                pubDOI = T.pub1_DOI;
            else
                if ~isempty(strfind(T.pub1_DOI,''''))
                    match  = regexp(pubDOI,'''','split');
                    pubDOI = match{1};
                else
                    match  = regexp(pubDOI,',\s','split');
                    pubDOI = match{1};
                end
            end
        end
        if isempty(pubDOI) || ~isempty(strfind(pubDOI,'NA'))
            pubDOI = 'NA';
        end

        if iscell(T.pub1_author)
            T.pub1_author = T.pub1_author{1};
        end
        % Actual plot
        fig('summary'); clf;
        set(gcf,'position',[229   0   800   600])

        %% Title
        % ht = suptitle({'Instructions: If metadata are missing or incorrect, please update on the metadata spreadsheet','Note: Additional metadata available on the spreadsheet'});
        dataSetName = T.dataSetName; %strrep(T.dataSetName,'_','\_');
        pages2kID   = T.paleoData_pages2kID; %strrep(T.paleoData_pages2kID,'_','\_');
        ht = suptitle([pages2kID ': ', dataSetName]);
        set(ht,'FontName','Helvetica','FontWeight','Bold','position',[0.5  -0.03 0])
        axes('position',[.8 .95 .2 .1])
        at = text(0,0,['\it{created: ' date '}']); set(at,style_a{:})
        set(gca,'color','none','XTick',[],'YTick',[],'Visible','off')

        %% first plot the raw data
        % Time series
        subplot(2,3,1:2)
        if ~isempty(ta) % If record is available over CE
            ho = plot(to,Xo,'marker','o','color',rgb('Gray'),'MarkerSize',5,'linestyle','none'); hold on;
            %ha = plot(ta,Xa(~isnan(Xa)),'color',Graph{p_code(r),1},'linewidth',2); hold off;
            ha = plot(year,Xa,'color',Graph{p_code(r),1},'linewidth',2); hold off;

            % if yearMin(r) ~= yearMax(r) again, relic from a dark age
            xlim([max(1,yearMin(r)),yearMax(r)]);
            lab{1}  = 'original'; lab{2} = 'annualized';
            %if resMed(r) <= 5 && resAvg(r) <=5
            legend([ho(1), ha],lab{:},'location','best'), legend boxoff
            %end
        else
            plot(to,Xo,'o-','color',Graph{p_code(r),1},'MarkerSize',8,'MarkerFaceColor',Graph{p_code(r),1},'linestyle','none','linewidth',2);
        end
        ylab    = [T.paleoData_variableName ' (' removeLeadingAndTrailingSpaces(T.paleoData_units) ')'];
        fancyplot_deco('','Year (CE)',ylab,14,'Helvetica')

        %% high resolution records
        if n_raw_samples(r) >= samples_thresh & n_ann_samples(r) >= samples_thresh
            r_jja = nan(ny,nx);
            r_djf = nan(ny,nx);
            r_mat = nan(ny,nx);
            instCalib(r) = 1;
            for l = idx_neigh{r}
                r_jja(ilat(l),ilon(l)) = rho_jja(l,r);
                r_djf(ilat(l),ilon(l)) = rho_djf(l,r);
                r_mat(ilat(l),ilon(l)) = rho_mat(l,r);
            end
            %============ meta data info
            axes('position',[.66 .55 .4 .4])
            ylim([5.5 11]), xlim([0 10])

            % text(0,10.5,['dataSetName: \bf{' T.dataSetName '}'], 'Interpreter', 'tex');
            text(0,10.5,['archiveType: \bf{' archive{r} '}'], 'Interpreter', 'tex')
            text(0,10,['PAGES2k region: ' T.geo_pages2kRegion], 'Interpreter', 'none')
            text(0,9.5,['Authors: ' T.pub1_author], 'Interpreter', 'none')
            text(0,9.0,['Year: ' pub1_year ', DOI: ' pubDOI], 'Interpreter', 'none')
            text(0,8.5,['Measured Variable: \bf{' strrep(T.paleoData_variableName,'_','\_') '}'], 'Interpreter', 'tex')
            text(0,8.0,['Units: ' T.paleoData_units], 'Interpreter', 'none')
            text(0,7.5,'climateIntepretation: ', 'Interpreter', 'none')
            text(1,7.0,['climateVariable: ' T.climateInterpretation_variable], 'Interpreter', 'none')
            text(1,6.5,['climateVariableDetail: ' T.climateInterpretation_variableDetail], 'Interpreter', 'none')
            text(1,6.0,['Seasonality: ' seasonality], 'Interpreter', 'none')
            text(1,5.5,['Direction: ' T.climateInterpretation_interpDirection], 'Interpreter', 'none')
            text(0,5.2,['useInGlobalTemperatureAnalysis: ' T.paleoData_useInGlobalTemperatureAnalysis], 'Interpreter', 'none')
            set(gca,'color','none','XTick',[],'YTick',[],'Visible','off')

            %%============  Correlation
            ha = subplot(234);
            m_proj('ortho','lat',round(p_lat(r)),'long',round(p_lon(r)));
            h = m_pcolor(x,y,r_mat); hold on; caxis(xrange);
            set(h,'edgecolor','none'),
            m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
            m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',Graph{p_code(r),1},'linewidth',2,'MarkerSize',12);
            m_line(locs(sig_mat(:,r)==0,1),locs(sig_mat(:,r)==0,2),'color',rgb('DimGray'),'Marker','x','linestyle','none','markersize',3)
            m_coast('color','k');
            colormap(cejulien2(nc))
            str_loc = sprintf('%6.2f',rho_mat(idx_loc(r,1),r));
            if sig_mat(idx_loc(r,1),r) == 1
                str_loc = ['{\bf ' str_loc '}'];
            end
            title(['Proxy vs. MAT ( R_{loc} =',str_loc,')'],'FontName' ,'Helvetica','FontSize',11);

            % ===JJA
            ha = subplot(235);
            m_proj('ortho','lat',round(p_lat(r)),'long',round(p_lon(r)));
            % m_proj('robinson');
            h = m_pcolor(x,y,r_jja); hold on; caxis(xrange);
            set(h,'edgecolor','none'),
            m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
            m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',Graph{p_code(r),1},'linewidth',2,'MarkerSize',12);
            m_line(locs(sig_jja(:,r)==0,1),locs(sig_jja(:,r)==0,2),'color',rgb('DimGray'),'Marker','x','linestyle','none','markersize',3)
            m_coast('color','k');
            str_loc = sprintf('%6.2f',rho_jja(idx_loc(r,1),r));
            if sig_jja(idx_loc(r,1),r) == 1
                str_loc = ['{\bf ' str_loc '}'];
            end
            title(['Proxy vs. JJA ( R_{loc} =',str_loc,')'],'FontName' ,'Helvetica','FontSize',11);

            % ===DJF
            ha = subplot(236); pos = get(ha,'position');
            m_proj('ortho','lat',round(p_lat(r)),'long',round(p_lon(r)));
            % m_proj('robinson');
            h = m_pcolor(x,y,r_djf); hold on; caxis(xrange);
            set(h,'edgecolor','none'),
            m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
            m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',Graph{p_code(r),1},'linewidth',2,'MarkerSize',12);
            m_line(locs(sig_djf(:,r)==0,1),locs(sig_djf(:,r)==0,2),'color',rgb('DimGray'),'Marker','x','linestyle','none','markersize',3)
            m_coast('color','k');
            str_loc = sprintf('%6.2f',rho_djf(idx_loc(r,1),r));
            if sig_djf(idx_loc(r,1),r) == 1
                str_loc = ['{\bf ' str_loc '}'];
            end
            ht = title(['Proxy vs. DJF ( R_{loc} =',str_loc,')'],'FontName' ,'Helvetica','FontSize',11);
            posT = get(ht,'position');
            hc = colorbar2('horiz','Correlation');
            set(hc,'position',[0.367 0.06 0.3 0.02],'Xtick',-0.5:0.25:0.5)
            set(ha,'position',pos)
            set(ht,'position',posT);
        else
            %% low resolution records
            % Find n closest hi-res proxy neighbors
            jj = setdiff([1:np]',r); % index of other proxies
            lat_p = p_lat(r);  lon_p = p_lon(r);
            lat_j = p_lat(jj); lon_j = p_lon(jj);
            res_j = res(jj); idx_r = idx_pp(r,:);
            idx_p = idx_r(res_j(idx_r)<=5);  n_p = numel(idx_p);
            prox_n = standardize(proxy_ann(:,jj(idx_p)));
            rho_n{r}  = zeros(n,1);
            % interpolate LR proxy record
            ti = [min(ta):max(ta)];
            Xi = pchip(ta,Xa(~isnan(Xa)),ti);
            k = 0;  in = 0; n_idx = [];
            while in < n && k < n_p
                k = k +1;
                noNaN = ~isnan(prox_n(:,k));
                tn = year(noNaN); Xn = prox_n(noNaN,k);
                tni = [min(tn):max(tn)]; Xni = interp1(tn,Xn,tni);
                % low_pass filter
                try
                    Xf = hepta_smooth(Xni,1/res(r));
                catch ME1
                    try
                        Xf = lowpassmin(Xni(:),1/res(r));
                    catch ME2
                        disp(['Record r=',int2str(r),' is too ghetto to lowpass'])
                        Xf = Xni;
                    end
                end
                %define overlap
                overlap = intersect(tni,ti);
                no = numel(overlap);
                if no >= samples_thresh
                    in = in +1; n_idx(in) = k;
                    display(['computing correlations for neighbor #', int2str(in), '/',int2str(n)]);
                    [rho_n{r}(in),signif_n{r}(in)] = corr_isospec(Xi(ismember(ti,overlap))',Xf(ismember(tni,overlap))');
                    if  no <50
                        thick(in) = 4;
                    elseif no>=50 && no <150
                        thick(in) = 8;
                    elseif no>=150 && no <500
                        thick(in) = 12;
                    else
                        thick(in) = 16;
                    end
                end
            end

            idx_pn = idx_p(n_idx); ds{r} = d_pp(r,idx_pn);
            %============ meta data info
            axes('position',[.66 .55 .4 .4])
            ylim([5.5 11]); xlim([0 10])
            text(0,10.5,['archiveType: \bf{' T.archiveType '}'], 'Interpreter', 'tex')
            text(0,10,['PAGES2k region: ' T.geo_pages2kRegion], 'Interpreter', 'none')
            text(0,9.5,['Authors: ' T.pub1_author], 'Interpreter', 'none')
            text(0,9,['Year: ' T.pub1_pubYear ', DOI: ' pubDOI], 'Interpreter', 'none')
            text(0,8.5,['parameter: \bf{' strrep(T.paleoData_variableName,'_','\_') '}'], 'Interpreter', 'tex')
            text(0,8,['units: ' T.paleoData_units], 'Interpreter', 'none')
            text(0,7.5,'climateIntepretation: ', 'Interpreter', 'none')
            text(1,7,['climateVariable: ' T.climateInterpretation_variable], 'Interpreter', 'none')
            text(1,6.5,['climateVariableDetail: ' T.climateInterpretation_variableDetail], 'Interpreter', 'none')
            text(1,6,['Seasonality: ' seasonality], 'Interpreter', 'none')
            text(1,5.5,['Direction: ' T.climateInterpretation_interpDirection], 'Interpreter', 'none')
            text(1,5.0,['QC Notes: ' T.paleoData_QCnotes], 'Interpreter', 'none')

            if strcmp(units{r},'deg C') || strcmp(units{r},'\u00b0C') || strcmp(units{r},'C') || strcmp(units{r},'(degrees C) ')| strcmp(units{r},'degC') | strcmp(units{r},'degrees Celsius');
                text(0,4,'{\bf Calibration Information}', 'Interpreter', 'tex')
                text(1,3.5,['calibration_equation: ' T.calibration_equation], 'Interpreter', 'none')
                text(1,3.0,['calibration_reference: ' T.calibration_reference], 'Interpreter', 'none')
                text(1,2.5,['calibration_uncertainty: ' T.calibration_uncertainty], 'Interpreter', 'none')
                text(1,2.0,'calibration_uncertaintyType: RMSE', 'Interpreter', 'none')
                text(1,1.5,['calibration_notes: ' T.calibration_notes], 'Interpreter', 'none')
            end

            set(gca,'color','none','XTick',[],'YTick',[],'Visible','off')
            %%============  Proxy location and correlation to HR neighbors
            ha=subplot(234); cla;

            if abs(lat_p)>60  % polar cases
                aminlat = round(min(abs(p_lat(idx_pn))));
                viewRadius = max([90-aminlat+5,20]);
                m_proj('stereographic','lat',round(lat_p),'long',round(lon_p),'radius',viewRadius); hold on
                m_grid('linewi',1,'tickdir','out','yticklabels',[]);
            else              % other cases
                maxAngle = max(ds{r}/110./cosd(p_lat(jj(idx_pn))));
                if maxAngle >= 60
                    m_proj('ortho','lat',round(lat_p),'long',round(lon_p)); hold on
                elseif maxAngle <= 10
                    m_proj('stereographic','lat',round(lat_p),'long',round(lon_p),'radius',10); hold on
                else
                    m_proj('stereographic','lat',round(lat_p),'long',round(lon_p),'radius',maxAngle); hold on
                end
                m_grid('box','off','xtick',9,'ytick',6,'xticklabels',[],'yticklabels',[]);
            end
            m_coast('patch',rgb('WhiteSmoke'));
            for k=1:length(P.ncst), % plot province boundaries
                m_line(P.ncst{k}(:,1),P.ncst{k}(:,2),'color',rgb('SlateGray'),'linewidth',.5);
            end
            for k=1:length(C.ncst), % plot country boundaries
                m_line(C.ncst{k}(:,1),C.ncst{k}(:,2),'color',rgb('SlateGray'),'linewidth',.5);
            end

            if in > 0
                % assign color based on correlations
                col = t2c(abs(rho_n{r}),51,'pmk',[0,1]);
                for k=1:in
                    rk = idx_pn(k);
                    [lat,lon,dis] = greatcircle(lat_p,lon_p,lat_j(rk),lon_j(rk),5);
                    if signif_n{r}(k)
                        m_line(lon,lat,'Color',col(k,:),'linewidth',3,'linestyle','-');
                    else % make it dashed if insignificant
                        %[H1,H2,H3] = m_line_fewer(lon,lat,'Color',col(k,:),'linewidth',1,'linestyle',':');
                        %m_scatter(lon,lat,'SizeData',25,'MarkerEdgeColor',col(k,:),);
                        m_line(lon,lat,'Color',col(k,:),'linewidth',2,'linestyle','--');
                    end
                    % plot proxy locations
                    m_line(lon_j(rk),lat_j(rk),'marker',Graph{p_code(jj(rk)),2},'MarkerEdgeColor','k','MarkerFaceColor',Graph{p_code(jj(rk)),1},'linewidth',1,'MarkerSize',9);
                end
                % plot proxy itself
                m_line(lon_p,lat_p,'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',Graph{p_code(r),1},'linewidth',1,'MarkerSize',12);
                colorbar2('horiz','Absolute correlation')
                title('Correlation to HR neighbors','FontName','Helvetica','FontSize',14,'Fontweight','bold')
                %text(-1,-1,'(solid=significant)','Interpreter', 'none')
                %  correlation vs distance
                subplot(235); cla; hold on;
                pc = p_code(jj(idx_pn));
                for k = 1:in
                    rk = idx_pn(k);
                    pb(k)=plot(ds{r}(k),abs(rho_n{r}(k)),'MarkerFaceColor',Graph{pc(k),1},'MarkerEdgeColor','none','Markersize',thick(k),'Marker',Graph{pc(k),2},'linewidth',2);
                end
                ylim([0 1]);
                if sum(signif_n{r}) > 0
                    set(pb(signif_n{r}),'MarkerEdgeColor','k')
                end
                hold off
                fancyplot_deco('Proximity Effect','Distance (km)','Absolute Correlation',14,'Helvetica')
            else
                % plot proxy itself
                m_line(lon_p,lat_p,'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',Graph{p_code(r),1},'linewidth',1,'MarkerSize',12);
                title('Proxy location','FontName','Helvetica','FontSize',14,'Fontweight','bold')
            end
            % create legend
            subplot(236);
            set(gca,'visible','off'),axis([0 1 0 1])
            text(.3,.2,'# years of overlap','fontsize',14)
            line(.05,.10,'marker','o','MarkerFaceColor',rgb('DimGray'),'MarkerEdgeColor','none','MarkerSize',4)
            text(0,0,'<50','fontsize',10)
            line(.35,.10,'marker','o','MarkerFaceColor',rgb('DimGray'),'MarkerEdgeColor','none','MarkerSize',8)
            text(.25,0,'50 to 150','fontsize',10)
            line(.65,.10,'marker','o','MarkerFaceColor',rgb('DimGray'),'MarkerEdgeColor','none','MarkerSize',12)
            text(.55,0,'150 to 500','fontsize',10)
            line(.95,.10,'marker','o','MarkerFaceColor',rgb('DimGray'),'MarkerEdgeColor','none','MarkerSize',16)
            text(.90,0,'>500','fontsize',10)
        end
        clear n_idx

        %% export
        if options.export
            region  = T.geo_pages2kRegion;
            dirname = ['../../figs/qc/' region];
            if ~exist(dirname,'dir') % if necessary
                system(['mkdir ', dirname]); % create regional directory
            end
            filen=[dirname '/pages_2k_phase2_record_', sprintf('%03d',r)];
            export_fig([filen '.pdf'],'-r200','-cmyk','-nocrop','-painters')
            %hepta_figprint([filen '.eps'],100)
        end
    end
    % pause;
end

fn = ['../data/PAGES2k_v' vers '_unpack'];

save(fn,'instCalib','rho_n','signif_n','ds','n_raw_samples','n_ann_samples','-append','-v6');
