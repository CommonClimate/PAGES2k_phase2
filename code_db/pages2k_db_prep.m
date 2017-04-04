function pages2k__db_prep(vers,options)
%  PAGES2k__db_PREP: gathers various data files into a structure
%   friendly for analysis and reconstruction. Also, plots proxy map
%  after 3 levels of screening (reg, fdr, local), and saves to EPS if export = 1.
%
%  INPUT:   vers (version of the database, string)
%           options, including:
%                 dtype (pre-processing of timeseries, options are:
%                       'detrend' : linear detrending,
%                       'diff1'   : first-difference of the original time series,
%                       'original': no detrending
%                 norm_p : boolean flag for transformation to normality
%                 method : screening method (string), options are:
%                       'isospec' : isospectral test on Person correlations
%                       'knockoff': Knockoff filter
%                 export : boolean flag for figure export to EPS.
%
%
%  OUTPUT:  data structure  {pages2k_hadcrut4_ 'dtype' _ 'norm_string' _ 'vers' .mat}
%      contains 4 sub-structures:
%        1) pages2k, containing the raw & processes PAGES 2k dataset, + a
%        few ancillary quantities useful for analysis
%        2) hadcrut4, containing temperature field for annual, jja and djf
%        periods, as well as area-weighted global mean and coordinates
%        (time, longitude, latitude).
%        3) gridcorr, containing information from the grid-centric correlation analysis (see pages2k_grid_analysis.m)
%        4) pcorr, containing information from the proxy-centric correlation analysis (see pages2k__db_screen.m)
%

year = [];
addpath(genpath('../utilities'));
%load JEG_graphics;

fopts = fieldnames(options);
if strmatch('norm_p', fopts)
    norm_p = options.norm_p;
else
    norm_p = 1;
end

if strmatch('dtype', fopts)
    dtype = options.dtype;
else
    dtype = 'original';
end

if strmatch('method', fopts)
    method = options.method;
else
    method = 'isospectral';
end


if norm_p
    norm_string = 'normal';
else
    norm_string = 'raw';
end
% unpacked database
load(['../data/PAGES2k_v' vers '_unpack.mat'])

% proxy-centric correlations
fname1 = ['../data/pages2k_hadcrut4_corr_' dtype '_' vers '_' norm_string '_' method '.mat'];
pcorr  = load(fname1);
%
fname2 = ['../data/pages2k_hadcrut4_gridcorr_' dtype '_' vers '_' norm_string '_annual.mat'];
gridcorr = load(fname2,'ggrid','ulat','ulon','localCorrelation','localSignificance') % only get the variables we need, to


%%  Apply basic quality control:
%  proxies need to be available at least prior to 1850AD, and with a total of at least 20
%  observations
avl      = ~isnan(proxy_ann);
keep_qchf = find(resMed' <= 5 & resAvg' <= 5 & yearMin' <= 1850 & sum(avl,1) >= 20 & yearMax' >= 1920);
keep_qclf = find(resMed' > 5 & resAvg' > 5& yearMin' <= 1850 & sum(avl,1) >= 20);%

%% Proxy

% First, remove calcification rate fields (calci_raw, calci_anom)
% TO DO: check if this is still needed
%calcification = find(~cellfun('isempty',strfind({S.paleoData_variableName},'calci'))...
% & cellfun('isempty',strfind({S.paleoData_variableName},'calci_mean')));
% all associated with 10.1126/science.1214570

% TAKE OUT BOREHOLES: should be reserved for validation.
ma = {S.archiveType};
noBorehole  = find(~strcmp(ma,'borehole'));

radius    = linspace(500,options.search_radius,4);
nrad = length(radius);
n_fdr = zeros(nr, nrad);
qBH   = 0.05;
%seas = {'mat','djf','jja'}; nseas = length(seas);

for r = 1:nr % loop over records
    for ir = 1:nrad % loop over radii
        pct = sprintf('%3.2g',r/nr*100);
        disp(['Applying FDR procedure for radii ' int2str(radius) ', ' pct '% done'])
        isWithin = (pcorr.d_pt(r,:) <= radius(ir))';
        % index of significant temperature neighbors
        pvals = pcorr.pval_mat(isWithin,r);
        pvals_noz = pvals(~isnan(pvals)); % subsample the appropriate set of pvals
        if ~isempty(pvals_noz)
            n_fdr(r,ir)  = fdr(pvals_noz,qBH,'original','mean');
            signif_fdr_mat(r,ir)  = (n_fdr(r,ir)>0);
        else
            signif_fdr_mat(r,ir)  = 0;
        end
    end
end

%% Export data matrix
pages2k.proxy_ann     = proxy_ann;
if norm_p
    pages2k.proxy_n   = proxy_n;
end
%
if  options.InterpSuperAnn
    pages2k.proxy_a   = proxy_a;
    if norm_p
        pages2k.proxy_na  = proxy_na;
    end
end

pages2k.proxy_djf = proxy_djf;
pages2k.proxy_jja = proxy_jja;
pages2k.year      = year(year>0);
pages2k.avail     = avail;
pages2k.archive   = archive;
pages2k.archiveType  = archiveType;
pages2k.loc        = [p_lon;p_lat]';
pages2k.keep_qchf  = keep_qchf;
pages2k.keep_qclf  = keep_qclf;
pages2k.screen_loc = pcorr.screen_loc;
pages2k.screen_reg = pcorr.screen_reg;
pages2k.screen_fdr = pcorr.screen_fdr;
pages2k.radius     = radius;
pages2k.version    = vers;
pages2k.detrend    = dtype;
pages2k.S          = S;
pages2k.yearMin    = yearMin;
pages2k.yearMax    = yearMax;
pages2k.resMed     = resMed;
pages2k.resAvg     = resAvg;
pages2k.resMax     = resMax;
pages2k.resMin     = resMin;
pages2k.signif_fdr_mat = signif_fdr_mat;
pages2k.signif_n   = signif_n;
pages2k.rho_n      = rho_n;
pages2k.radius     = radius;
pages2k.noBorehole = noBorehole;
pages2k.p_code     = p_code;
pages2k.Graph      = Graph;
pages2k.edgec      = edgec;
pages2k.instCalib  = instCalib;
pages2k.ds         = ds;

nr = length(S);
%% include results in S
for r = 1:nr
    S(r).resMed = resMed(r);
    S(r).resMin = resMin(r);
    S(r).resMax = resMax(r);
    S(r).resAvg = resAvg(r);
    S(r).signifHRneighbors = signif_n{r};
    S(r).screenLoc = ismember(r,pcorr.screen_loc{1});
    S(r).screenReg = ismember(r,pcorr.screen_reg{1});
    S(r).screenFDR = signif_fdr_mat(r,4);
    S(r).yearMin    = yearMin(r);
    S(r).yearMax    = yearMax(r);
    S(r).paleoData_valuesRange = range(S(r).paleoData_values);
end
disp('Saving results to unpack file')
save(['../data/PAGES2k_v' vers '_unpack.mat'],'S','-append')


% output basic stats
fileID = fopen(['../data/pages2k_db_stats.tex'],'w');
fprintf(fileID,'PAGES2k essential stats\n');
fprintf(fileID,'================ COMPOSITION =================\n');
for a = 1:na % loop over archive types
    nArch(a)   = sum(p_code == a);
    pct(a)     = 100*nArch(a)/nr;
    fprintf(fileID,'%s %i %s %3.4g %s\n',[archiveType{a} ': '],nArch(a), ' records (',pct(a),'%)')
end
fprintf(fileID,'============== RESOLUTION ===============\n');
fprintf(fileID,'%s %i %s\n','The  mean resolution of non-tree archives is ', round(mean(resMed(p_code < 11))), 'years');
fprintf(fileID,'%s %3.2g %s\n','The  median resolution of non-tree archives is ', median(resMed(p_code < 11)), 'years');
fprintf(fileID,'%s %3.2g %s %3.2g %s\n','For sedimentary archives these numbers are ', mean(resMed(p_code == 7 | p_code == 8)), ' and ', median(resMed(p_code == 7 | p_code == 8)), ', respectively');
fprintf(fileID,'============= LENGTH  ==========\n');
Span = yearMax - yearMin; Span(Span>2000) = 2000;
fprintf(fileID,'%s %g %s\n','Minimum span:', min(Span),' years')
fprintf(fileID,'%s %g %s\n','Maximum span:', max(Span),' years')
fprintf(fileID,'%s %g %s\n','Median span:', median(Span),' years')
fprintf(fileID,'%s %g %s\n','Mean span (arithmetic):', mean(Span),' years')
fprintf(fileID,'%s %g %s\n','Mean span (geometric):', geomean(Span),' years')
fprintf(fileID,'=========== CORRELATION TO TEMPERATURE ================\n');
fprintf(fileID,'%i %s\n',numel(pcorr.screen_reg{1}), ' records retained by REG screening')
fprintf(fileID,'%i %s\n',numel(pcorr.screen_fdr{1}), ' records retained by FDR screening')
fprintf(fileID,'%i %s\n',numel(pcorr.screen_loc{1}), ' records retained by LOC screening')
fprintf(fileID,'=========== CORRELATION TO PROXIES ================\n');
nonCalibN = sum(instCalib);
nonCalib = ~cellfun(@isempty,signif_n);
proxySig = signif_n(nonCalib);
nproxySig = cellfun(@sum,proxySig);
fprintf(fileID,'%i %s\n',nonCalibN, ' records are not calibratable. Of those:')
fprintf(fileID,'%i %s\n',sum(nproxySig>0),' records are correlated to HR neighbors')
%
fclose(fileID)  % close file


%% HadCRUT4
tmp     = load('../data/had4med_graphem_sp70');
tmp_ann = load('../data/had4med_graphem_sp70_annual.mat');
%ti   = unique(tmp.tvec(1:1968,1));
temp = tmp.Xf(1:1968,tmp.idx); % 164 full calendar years from Jan 1850 - Dec 2013

% Annual
annual = (temp(1:12:end,:) + temp(2:12:end,:) + temp(3:12:end,:) + ...
    temp(4:12:end,:) + temp(5:12:end,:) + temp(6:12:end,:) + ...
    temp(7:12:end,:) + temp(8:12:end,:) + temp(9:12:end,:) + ...
    temp(10:12:end,:) + temp(11:12:end,:) + temp(12:12:end,:))./12;
% JJA
summer = (temp(6:12:end,:) + temp(7:12:end,:) + temp(8:12:end,:))./3;
% DJF
winter = (temp(12:12:end,:) + temp(1:12:end,:) + temp(2:12:end,:))./3;
%
hadcrut4.mat = annual;
hadcrut4.jja = summer;
hadcrut4.djf = winter;
hadcrut4.loc = tmp_ann.loc;
hadcrut4.t   = tmp_ann.t;
hadcrut4.gmean = tmp_ann.gmean;

disp('Bundling results in one big file')
dn = ['../data/pages2k_hadcrut4_' dtype '_' norm_string '_' vers '.mat'];

if ismac & ~verLessThan('matlab','R2016a') % horrible hack for JEG's machine. Mathworks, you suck!
    save(dn,'hadcrut4','pages2k','gridcorr','pcorr','-v6')
else
    save(dn,'hadcrut4','pages2k','gridcorr','pcorr')
end



% Visualization of the screened networks
% ======================================
FontName = 'Helvetica';
set(0,'defaultAxesFontName', FontName)
set(0,'defaultTextFontName', FontName)
style_t = {'FontName',FontName,'Fontweight','Bold','Fontsize',14};
style_l = {'FontName',FontName};
screen = {'fdr','loc','reg'}; ns = length(screen);
for i = 1:ns
    keep = eval(['pages2k.screen_' screen{i} '{1}']); % 1 for MAT
    p_code_k = p_code(keep);  avec = unique(p_code_k);
    aType=archiveType(avec); nak = numel(aType);
    Graph_k = Graph(avec,:); avail_k = avail(:,keep);
    nproxy = zeros(ny,nak); pind = zeros(nak,1); clear h;
    for k=1:nak % loop over archive types
        a = avec(k);
        nproxy(:,a) = sum(~isnan(avail_k(:,p_code_k == a)),2);
        pind(a) = find(p_code_k == a,1,'first');
    end
    pind = pind(pind~=0);

    n1000 = ceil(sum(nproxy(year == 1000,:))/10)*10;
    nr    = length(unique({S(keep).dataSetName})); % # of sites
    versl = strrep(vers,'_','/');
    fig('PAGES 2K screened'), clf;
    set(gcf,'PaperPositionMode','auto')
    %set(gcf,'position',[10 10 791 550])
    orient landscape
    % plot spatial distribution
    hmap=axes('Position', [.05 0.45 0.75 0.5]);
    m_proj('Robinson','clong',10);
    m_coast('patch',[.9 .9 .9]);
    m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',4,'fontname','Helvetica');
    % loop over records
    for j = 1:length(keep)
        r = keep(j);
        hk(j) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',[7],'linestyle','none');
    end
    text(-2,1.75,['Screened PAGES2k network (',screen{i},', ', int2str(length(keep)) , ' records from ', int2str(nr), ' sites)'],style_t{:});
    % legend
    %hk = h(keep);
    hl = legend(hk(pind),aType,'location',[.84 .6 .1 .2],style_l{:});
    set(hl, 'FontName', 'Helvetica','box','off');
    % temporal availability
    hstack=axes('Position', [0.1 0.1 0.8 0.29]);
    cmap=cell2mat(Graph(:,1));
    colormap(cmap);
    area(year,nproxy,'EdgeColor','w'), set(gca,'YAxisLocation','Right');
    xlim([1 2000])
    fancyplot_deco('','Year (CE)','# proxies');
    title('Temporal Availability',style_t{:})
    % inset
    frac=.5;
    hstackin=axes('Position', [0.1 0.2 frac*.8 0.14]);
    area(year,nproxy,'EdgeColor','w')
    axis([1 1000 0 n1000])
    set(hstackin,'xtick',[],'box','off','TickDir','out','TickLength',[.02 .02],'YMinorTick','on', 'YGrid','on')
    set(hstackin,'YAxisLocation','Right')
    set(hstackin,'ytick',[0:10:n1000],'FontName','Helvetica','YColor', [.3 .3 .3])
    title('First Millennium',style_l{:})
    if options.export
        export_fig(['../figs/PAGES2k_' vers '_spacetime_' screen{i} '.pdf'],'-r200','-cmyk','-painters')
        %hepta_figprint(['../../figs/synopsis/PAGES2K_phase2_' vers '_spacetime_' screen{i}])
    end
    clear hk
end
