%
%   Reconstruct global mean temperature from proxy composite
%     and investigate robustness to various choices.
%    All output is saved in 'f_out' (defined in 'workflow' script);
%  

% cosmetic definitions
style_t = {'FontName','Helvetica','FontSize',14,'FontWeight','bold'};
style_l = {'FontName','Helvetica','FontSize',12,'FontWeight','bold'};

% define a string that encodes these options
if norm_p == 1
    Gauss = 'normal';
else
    Gauss = 'raw';
end
% if detrend == 1
%     dts = 'CoralO18Detrend';
% else
%     dts = 'NoCoralO18Detrend';
% end

opstring =  ['PAGES2k_v' vers '_' sifting_style '_' screenHR_style '_' Gauss];

% load data matrices 
% ==================

% load essential arrays
t = pages2k.year;
tce  = pages2k.year(t >= tStart & t <= tEnd); nce = length(tce);

T = pages2k.S; nr = length(T);
names = {T.dataSetName};
yearMin = pages2k.yearMin;
yearMax = pages2k.yearMax;
resMed  = pages2k.resMed;
resAvg  = pages2k.resAvg;
resMax  = pages2k.resMax;
p_code  = pages2k.p_code; 
Graph   = pages2k.Graph; 
p_lon   = pages2k.loc(:,1);
p_lat   = pages2k.loc(:,2);
edgec   = pages2k.edgec;
archive = pages2k.archive;
      S = pages2k.S; 


incl = [1:nr]; 
% account for various pre-processing choices
if norm_p && detrend
    proxy = pages2k.proxy_nda(ismember(t,tce),incl);
elseif norm_p && ~detrend
    proxy = pages2k.proxy_na(ismember(t,tce),incl);
elseif ~norm_p && detrend
    proxy = pages2k.proxy_da(ismember(t,tce),incl);
else
    proxy = pages2k.proxy_a(ismember(t,tce),incl);
end

[ny,nr] = size(proxy);

% exploit temperature interpretation
sgn = cellfun(@lower,{T.climateInterpretation_interpDirection},'UniformOutput', false);
sgn_vec = zeros(nr,1);
for r = 1:nr
    if strcmpi(sgn{r},'positive') || strcmpi(sgn{r},'p')  % KLUDGE
        sgn_vec(r) = +1;
    elseif strcmpi(sgn{r},'negative')
        sgn_vec(r) = -1;
    end
end

proxy_sgn = standardize(proxy).*repmat(sgn_vec',[nce, 1]);

%% apply basic filtering: 
%  - proxies need to be available at least prior to 1850AD, and with a total of at least 20 observations
%  - kick boreholes out (they are great for validation, but they are not proxies) 
navl      = sum(~isnan(proxy_sgn)); ma = {pages2k.S.archiveType};
idx_qchr = find(~strcmp(ma,'borehole') & resMed' <= 5 & resAvg' <= 5 & yearMin' <= 1850 & navl >= navlMin);
idx_qclr = find(~strcmp(ma,'borehole') & resMed'  > 5 & yearMin' <= 1850 & navl >= navlMin);%

%% screening for significant correlations.
% 1) Calibratable Proxies
scr_reg  = pages2k.screen_reg{1}; % MAT regional correlation screening
scr_fdr  = pages2k.screen_fdr{1}; % MAT regional correlation screening controlling for false discovery rate
scr_loc  = pages2k.screen_loc{1}; % MAT local correlation screening

idx_scr = eval(['scr_' screenHR_style]);
% merge indices of screened proxies
idx_qchr_scr = intersect(idx_scr,idx_qchr);

CalibFalse_sig = zeros(nr,1); signif_n = pages2k.signif_n;
for r = 1:nr
    if sum(signif_n{r})>1
        CalibFalse_sig(r) = 1;
    end
end
idx_qclr_screen = intersect(find(CalibFalse_sig),idx_qclr);

% %%
% fig('Screening or not'), clf
% set(gcf,'Position',[440   270   852   628])
% pmax = 600; % scale for # proxies
% ylims = [-1.5 2];
% % Unscreened HR
% ax1 = subplot(3,2,1)
% Xlab = ''; Ylab = {'# records','HR composite'};
% P = proxy_sgn(:,idx_qchr); ps = size(P,2); navl = sum(~isnan(P),2);
% col{1} = rgb('Silver'), col{2} = rgb('Blue');
% [ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,col);
% set(ax(1),'xLim',[tStart tEnd]);
% title('High Resolution (\Delta t <=5y), unscreened',style_t{:})
% set(ax1,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
% set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% % Regional screening
% ax2 = subplot(3,2,2)
% P = proxy_sgn(:,intersect(scr_reg,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
% [ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,col);
% set(ax(1),'xLim',[tStart tEnd]);
% title('High Resolution, regional screening (R < 2000km)',style_t{:})
% set(ax2,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
% set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% % Regional screening w/ FDR
% ax3 = subplot(3,2,3)
% P = proxy_sgn(:,intersect(scr_fdr,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
% [ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,col);
% set(ax(1),'xLim',[tStart tEnd]);
% title('High Resolution, FDR screening (R < 2000km)',style_t{:})
% set(ax3,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
% set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% % local screening
% ax4 = subplot(3,2,4) 
% P = proxy_sgn(:,intersect(scr_loc,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
% [ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,col);
% set(ax(1),'xLim',[tStart tEnd]);
% title('High Resolution, local screening',style_t{:})
% set(ax4,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
% set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% % Low-res unscreened
% ax5 = subplot(3,2,5), pmax = 100;
% Xlab = 'Year CE'; Ylab = {'# records','LR composite'};
% P = proxy_sgn(:,idx_qclr); ps = size(P,2); navl = sum(~isnan(P),2);
% col{1} = rgb('Silver'), col{2} = rgb('Red');
% [ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,col);
% set(ax(1),'xLim',[tStart tEnd]);
% title('Low resolution (\Delta t > 5y), unscreened',style_t{:})
% set(ax5,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
% set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% % Low-res, screened
% ax6 = subplot(3,2,6) 
% P = proxy_sgn(:,idx_qclr_screen); ps = size(P,2); navl = sum(~isnan(P),2);
% col{1} = rgb('Silver'), col{2} = rgb('Red');
% [ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,col);
% set(ax(1),'xLim',[tStart tEnd]);
% title('Low resolution, screened against HR neighbors',style_t{:})
% set(ax6,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
% set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% 
% %
% %linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'xy'), 
% %xlim([0 2010]); ylim([-1.5 2]); 
% %export_fig(['../figs/composite_ScreeningEffects_' opstring '_Unsmoothed.pdf'],'-r200','-cmyk','-nocrop')
% 
% hepta_figprint(['../figs/composite_ScreeningEffects_' opstring '_Unsmoothed'])
% 



%% produce proxy composite
% the signal appears much cleaner w/ screening. 
switch sifting_style
    case 'noSift'
        idx_q = [1:nr];
        p_lr = numel(resMed <= 5); p_hr = numel(resMed > 5);
    case 'qcOnly'
        idx_q    = unique(union(idx_qchr,idx_qclr,'stable'));
        p_lr = numel(idx_qclr); p_hr = numel(idx_qchr);
    case 'qcScreenAll'
        idx_q    = union(idx_qchr_scr,idx_qclr_screen,'stable');
        p_lr = numel(idx_qclr_screen); p_hr = numel(idx_qchr_scr);
    case 'qcScreenHR'
        idx_q    = union(idx_qchr_scr,idx_qclr,'stable');
        p_lr = numel(idx_qclr); p_hr = numel(idx_qchr_scr);
    case 'qcScreenLR'
        idx_q    = union(idx_qchr,idx_qclr_screen,'stable');
        p_lr = numel(idx_qclr_screen); p_hr = numel(idx_qchr);
end

% extract records of interest
proxy_q = proxy_sgn(:,idx_q);
% how many are there?
p_q  = size(proxy_q,2);

weights = cosd(p_lat(idx_q))';
proxy_qw = proxy_q.*repmat(weights,[nce 1]);
% ========== 
if lat_weight
   proxy_qs = standardize(proxy_qw);
else
   proxy_qs = standardize(proxy_q); 
end

% ADD WEIGHT BY LOCAL TEMP VARIANCE 


%% MAP THE RECORDS USED IN THE COMPOSITE
edgec = cell(nr,1);
edgec(1:nr) = {'none'};%{rgb('Black')};
edgec(p_code == 6) = {Graph{6,1}};
%
archiveType = pages2k.archiveType;
p_code_k = p_code(idx_q);  avec = unique(p_code_k);
aType  = archiveType(avec); nak = numel(aType);
Graph_k = Graph(avec,:); avail_k = pages2k.avail(t >= tStart & t <= tEnd,idx_q);
nproxy = zeros(ny,nak); pind = zeros(nak,1); 
for k=1:nak % loop over archive types
    a = avec(k);
    nproxy(:,a) = sum(~isnan(avail_k(:,p_code_k == a)),2);
    pind(a) = find(p_code_k == a,1,'first');
end
pind = pind(pind~=0);

n1000 = ceil(sum(nproxy(tce == 1000,:))/10)*10;
nsites    = length(unique({pages2k.S(idx_q).dataSetName})); % # of sites
versl = strrep(vers,'_','.');
fig('PAGES2k screened'), clf;
set(gcf,'position',[10 10 791 550])
orient landscape
% plot spatial distribution
hmap=axes('Position', [.05 0.45 0.75 0.5]);
m_proj('Robinson','clong',10);
m_coast('patch',[.9 .9 .9]);
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',8);
% loop over records
for j = 1:length(idx_q)
    r = idx_q(j);
    hk(j) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',[7],'linestyle','none');
end
text(-2,1.75,['Screened PAGES2k network (',sifting_style,' sifting , ', int2str(length(idx_q)) , ' records from ', int2str(nsites), ' sites)'],style_t{:});
% legend
hl = legend(hk(pind),aType,'location',[.84 .6 .1 .2]);
set(hl, 'FontName', 'Helvetica','box','off');
% temporal availability
hstack=axes('Position', [0.1 0.1 0.8 0.29]);
cmap=cell2mat(Graph(:,1));
colormap(cmap);
area(tce,nproxy,'EdgeColor','w'), set(gca,'YAxisLocation','Right');
xlim([1 2000])
fancyplot_deco('','Year (CE)','# proxies');
title('Temporal Availability',style_t{:})
% inset
frac=.5;
hstackin=axes('Position', [0.1 0.23 frac*.8 0.14]);
area(tce,nproxy,'EdgeColor','w')
axis([1 1000 0 n1000])
set(hstackin,'xtick',[],'box','off','TickDir','out','TickLength',[.02 .02],'YMinorTick','on', 'YGrid','on')
set(hstackin,'YAxisLocation','Right')
set(hstackin,'ytick',[0:50:n1000],'YColor', [.3 .3 .3])
title('First Millennium',style_l{:})
froot = ['../figs/' opstring '_mapRecordsIncluded'];
export_fig([froot '.pdf'],'-r200','-nocrop','-cmyk','-painters');
clear hk


%% mean composite
%p_comp = nmean(proxy_qs,2);
%p_std  =  nstd(proxy_qs,0,2);

save(f_out)







