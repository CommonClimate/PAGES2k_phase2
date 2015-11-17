%
%   Reconstruct global mean temperature from proxy composite
%     and investigate robustness to various choices.
%    All output is saved in 'fout' (defined in 'workflow' script);
%  

% cosmetic definitions
style_t = {'FontName','Helvetica','FontSize',14,'FontWeight','bold'};

% define a string that encodes these options
if norm_p == 1
    Gauss = 'Gaussianized';
else
    Gauss = 'NotGaussianized';
end
if detrend == 1
    dts = 'CoralO18Detrend';
else
    dts = 'NoCoralO18Detrend';
end

opstring =  [sifting_style '_'  Gauss '_' dts];

% load data matrices [TO DO: clean up so only the merged one is needed]
% ==================
% HadCRUT4.2 data interpolated as in https://www.authorea.com/users/18150/articles/20243/_show_article
hadcrut4 = load('./data/had4med_graphem_sp70_annual');
% PAGES2k database
pages2k = load(['./data/Pages2kPhase2Database_' vers '_unpack']);
% define output file
merged = load(['./data/pages2k_hadcrut4_original_' vers '.mat'])% merged proxy /temperature matrix
% load essential arrays
t = pages2k.year;
tce  = pages2k.year(t >= tStart & t <= tEnd); nce = length(tce);

T = pages2k.TS_temp; nr = length(T);
names = {T.paleoArchiveName};
yearMin = pages2k.yearMin;
yearMax = pages2k.yearMax;
resMed  = pages2k.resMed;
resMin  = pages2k.resMin;
resAvg  = pages2k.resAvg;
resMax  = pages2k.resMax;
p_code  = pages2k.p_code; 
Graph   = pages2k.Graph; 

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
sgn = lower({T.climateInterpretationInterpDirection});
sgn_vec = zeros(nr,1);
for r = 1:nr
    if strcmpi(sgn{r},'positive')
        sgn_vec(r) = +1;
    elseif strcmpi(sgn{r},'negative')
        sgn_vec(r) = -1;
    end
end
proxy_sgn = standardize(proxy).*repmat(sgn_vec',[nce, 1]);

%% apply basic filtering: 
%  - proxies need to be available at least prior to 1850AD, and with a total of at least 20 observations
%  - kick boreholes out (they are great for validation, but they are not proxies) 
navl      = sum(~isnan(proxy_sgn)); ma = {pages2k.TS_temp.archiveType};
idx_qchr = find(~strcmp(ma,'borehole') & resMed' <= 5 & resAvg' <= 5 & yearMin' <= 1850 & navl >= navlMin);
idx_qclr = find(~strcmp(ma,'borehole') & resMed'  > 5 & yearMin' <= 1850 & navl >= navlMin);%

%% screening for significant correlations.
% 1) Calibratable Proxies
radius   = merged.pages2k.radius;
scr_reg  = find(merged.pages2k.screen_reg{1}); % MAT regional correlation screening
scr_fdr  = find(merged.pages2k.signif_fdr_mat(:,radius == 2000)); % MAT regional correlation screening controlling for false discovery rate
scr_loc  = find(merged.pages2k.screen_loc{1}); % MAT local correlation screening

% merge indices of screened proxies
idx_qchr_scr_fdr = intersect(scr_fdr,idx_qchr);

CalibFalse_sig = zeros(nr,1); signif_n = pages2k.signif_n;
for r = 1:nr
    if sum(signif_n{r})>1
        CalibFalse_sig(r) = 1;
    end
end
idx_qclr_screen = intersect(find(CalibFalse_sig),idx_qclr);

%%
fig('Screening or not'), clf
set(gcf,'Position',[440   270   852   628])
pmax = 600; % scale for # proxies
ylims = [-1.5 2];
% Unscreened HR
ax1 = subplot(3,2,1)
Xlab = ''; Ylab = {'# records','HR composite'};
P = proxy_sgn(:,idx_qchr); ps = size(P,2); navl = sum(~isnan(P),2);
col{1} = rgb('Silver'), col{2} = rgb('Blue');
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution (\Delta t <=5y), unscreened',style_t{:})
set(ax1,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Regional screening
ax2 = subplot(3,2,2)
P = proxy_sgn(:,intersect(scr_reg,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution, regional screening (R < 2000km)',style_t{:})
set(ax2,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Regional screening w/ FDR
ax3 = subplot(3,2,3)
P = proxy_sgn(:,intersect(scr_fdr,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution, FDR screening (R < 2000km)',style_t{:})
set(ax3,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% local screening
ax4 = subplot(3,2,4) 
P = proxy_sgn(:,intersect(scr_loc,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution, local screening',style_t{:})
set(ax4,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Low-res unscreened
ax5 = subplot(3,2,5), pmax = 100;
Xlab = 'Year CE'; Ylab = {'# records','LR composite'};
P = proxy_sgn(:,idx_qclr); ps = size(P,2); navl = sum(~isnan(P),2);
col{1} = rgb('Silver'), col{2} = rgb('Red');
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('Low resolution (\Delta t > 5y), unscreened',style_t{:})
set(ax5,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Low-res, screened
ax6 = subplot(3,2,6) 
P = proxy_sgn(:,idx_qclr_screen); ps = size(P,2); navl = sum(~isnan(P),2);
col{1} = rgb('Silver'), col{2} = rgb('Red');
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('Low resolution, screened against HR neighbors',style_t{:})
set(ax6,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);

%
%linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'xy'), 
%xlim([0 2010]); ylim([-1.5 2]); 
%export_fig(['./figs/composite_ScreeningEffects_' opstring '_Unsmoothed.pdf'],'-r200','-cmyk','-nocrop')

hepta_figprint(['./figs/composite_ScreeningEffects_' opstring '_Unsmoothed'])

%% smooth things out
% =====================
% instrumental
tt    = hadcrut4.t;
gmean_f = hepta_smooth(hadcrut4.gmean(ismember(tt,tce)), f_hr);
tt = hadcrut4.t(ismember(tt,tce));
% proxies
proxy_f  = proxy_sgn;
op.crit_name = 'pct_var'; op.P = 1.00 ;

% manual fix for Wannamaker timeseries:
rWan = find(strcmp('AtlanticMaineWanamaker2008',pages2k.names));
resMax(rWan) = 100; 
resMin(rWan) = 0.99;
% HIGH RESOLUTION RECORDS
j = 0;
for r = idx_qchr
    j = j+1
    disp(['Smoothing record #',int2str(j) '/' int2str(numel(idx_qchr))])
    Z = proxy_sgn(:,r); tn = tce(~isnan(Z));
    interval = find(tce >= min(tn) & tce <= max(tn));
    ti = tce(interval); Zi = Z(interval);
    % find breaks
    [~,~,segment] = find_breaks(ti,Zi); nseg = numel(segment);
    gfrac = sum(isnan(Zi))/numel(Zi); % gap fraction
    Zf = NaN*Zi;   % dummy array
        
    if nseg > 1 && resMin(r)<1 && resMax(r)> 10 % if more than 1 segment, and large gaps
        for i=1:nseg;
            Zs = Zi(segment{i}); ts = ti(segment{i});
            try
                Zf(segment{i}) = hepta_smooth(Zs,f_hr); % LOWPASS FILTER
            catch
                Zf(segment{i}) = Zs; % don't smoothe
            end        
        end
    elseif sum(isnan(Zi)) > 0
        disp(['     SSAM infilling necessary for record #',int2str(r)])
        % pre-infill using singular spectrum analysis
        try
            [spec,eig_vec,PC,RC,Zssa,modes] = hepta_ssam(Zi,op);
            Zf = hepta_smooth(Zssa,f_hr);  % LOWPASS FILTER
        catch
            Zn = interp1(ti,Zi,ti,'pchip');
            Zf = hepta_smooth(Zn,f_hr);
        end
    else  
        Zf = hepta_smooth(Zi,f_hr);  % LOWPASS FILTER
    end   
    proxy_f(interval,r) = Zf;
    
    % plot things out
    if smoothExport
        fig('Smooth Sanity'), clf
        ho = plot(ti,Zi,'marker','o','color',rgb('Gray'),'MarkerSize',5,'linestyle','none'); hold on;
        ha = plot(ti,Zf,'color',Graph{p_code(r),1},'linewidth',2); hold off;
        xlim([max(1,yearMin(r)),yearMax(r)]);
        lab{1}  = 'annualized data'; lab{2} = [int2str(1/f_hr) 'y smoothed'];
        site_n = strrep(T(r).paleoArchiveName,'_',' ');
        ttl = [int2str(j),') record #',int2str(r),', ',T(r).archiveType,', ', site_n];
        ylab = [T(r).measurementShortName ' (' removeLeadingAndTrailingSpaces(T(r).measurementUnits) ')'];
        fancyplot_deco(ttl,'Year (CE)',ylab,14);
        legend(lab{:}), legend boxoff
        filen=['./figs/smoothing/pages_2k_phase2_record_' sprintf('%03d',r) '_' opstring '.pdf'];
        %pause
        export_fig(filen,'-r100','-cmyk','-painters','-nocrop')
    end
end

% LOW RESOLUTION RECORDS
% ======================
j = 0;
for r = idx_qclr
    j = j+1
    disp(['Smoothing record #',int2str(j) '/' int2str(numel(idx_qclr))])
    Z = proxy_sgn(:,r); tn = tce(~isnan(Z));
    interval = find(tce >= min(tn) & tce <= max(tn));
    ti = tce(interval); Zi = Z(interval);
    % LOWPASS FILTER
    Zf = hepta_smooth(Zi,f_lr);  
    proxy_f(interval,r) = Zf;
    % plot things out
    if smoothExport
        fig('Smooth Sanity'), clf
        ho = plot(ti,Zi,'marker','o','color',rgb('Gray'),'MarkerSize',5,'linestyle','none'); hold on;
        ha = plot(ti,Zf,'color',Graph{p_code(r),1},'linewidth',2); hold off;
        xlim([max(1,yearMin(r)),yearMax(r)]);
        lab{1}  = 'annualized data'; lab{2} = [int2str(1/f_lr) 'y smoothed'];
        site_n = strrep(T(r).paleoArchiveName,'_',' ');
        ttl = [int2str(j),') record #',int2str(r),', ',T(r).archiveType,', ', site_n];
        ylab = [T(r).measurementShortName ' (' removeLeadingAndTrailingSpaces(T(r).measurementUnits) ')'];
        fancyplot_deco(ttl,'Year (CE)',ylab,14);
        legend(lab{:}), legend boxoff
        filen=['./figs/smoothing/pages_2k_phase2_record_' sprintf('%03d',r) '_' opstring '.pdf'];
        %pause
        export_fig(filen,'-r100','-cmyk','-painters','-nocrop')
    end
end





%hepta_figprint(['./figs/composite_' opstring 'ScreeningEffects_Smoothed'])


%% produce proxy composite
% the signal appears much cleaner w/ screening. 
switch sifting_style
    case 'noSift'
        idx_q = [1:nr];
        p_lr = numel(resMed <= 5); p_hr = numel(resMed > 5);
    case 'qcOnly'
        idx_q    = union(idx_qchr,idx_qclr,'stable');
        p_lr = numel(idx_qclr); p_hr = numel(idx_qchr);
    case 'qcScreenAll'
        idx_q    = union(idx_qchr_scr_fdr,idx_qclr_screen,'stable');
        p_lr = numel(idx_qclr_screen); p_hr = numel(idx_qchr_scr_fdr);
    case 'qcScreenHR'
        idx_q    = union(idx_qchr_scr_fdr,idx_qclr,'stable');
        p_lr = numel(idx_qclr); p_hr = numel(idx_qchr_scr_fdr);
    case 'qcScreenLR'
        idx_q    = union(idx_qchr,idx_qclr_screen,'stable');
        p_lr = numel(idx_qclr_screen); p_hr = numel(idx_qchr);
end

% extract records of interest
proxy_q = proxy_sgn(:,idx_q);
% how many are there?
p_q  = size(proxy_q,2);

weights = cosd(pages2k.p_lat(idx_q));
proxy_fw = proxy_f(:,idx_q).*repmat(weights,[nce 1]);
% ========== 
if lat_weight
    proxy_fs = standardize(proxy_fw);%.*repmat(sgn_vec(idx_q)',[nce, 1]);
else
    proxy_fs = standardize(proxy_f(:,idx_q));%.*repmat(sgn_vec(idx_q)',[nce, 1]);  
end


%%
fig('Screening or not, after smoothing'), clf
set(gcf,'Position',[440   270   852   628])
pmax = 600; % scale for # proxies
ylims = [-1.5 2];
% Unscreened HR
ax1 = subplot(3,2,1)
Xlab = ''; Ylab = {'# records','HR composite'};
P = proxy_fs(:,idx_qchr); ps = size(P,2); navl = sum(~isnan(P),2);
col{1} = rgb('Silver'), col{2} = rgb('Blue');
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution (\Delta t <=5y), unscreened',style_t{:})
set(ax1,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Regional screening
ax2 = subplot(3,2,2)
P = proxy_fs(:,intersect(scr_reg,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution, regional screening (R < 2000km)',style_t{:})
set(ax2,'Ylim',[0 pmax],'Ytick',[0:100:pmax]); 
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Regional screening w/ FDR
ax3 = subplot(3,2,3)
P = proxy_fs(:,intersect(scr_fdr,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution, FDR screening (R < 2000km)',style_t{:})
set(ax3,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% local screening
ax4 = subplot(3,2,4) 
P = proxy_fs(:,intersect(scr_loc,idx_qchr)); ps = size(P,2); navl = sum(~isnan(P),2);
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('High Resolution, local screening',style_t{:})
set(ax4,'Ylim',[0 pmax],'Ytick',[0:100:pmax])
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Low-res unscreened
ax5 = subplot(3,2,5), pmax = 100;
Xlab = 'Year CE'; Ylab = {'# records','LR composite'};
P = proxy_fs(:,idx_qclr); ps = size(P,2); navl = sum(~isnan(P),2);
col{1} = rgb('Silver'), col{2} = rgb('Red');
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('Low resolution (\Delta t > 5y), unscreened',style_t{:})
set(ax5,'Ylim',[0 pmax],'Ytick',[0:100:pmax]);
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
% Low-res, screened
ax6 = subplot(3,2,6) 
P = proxy_fs(:,idx_qclr_screen); ps = size(P,2); navl = sum(~isnan(P),2);
col{1} = rgb('Silver'), col{2} = rgb('Red');
[ax,h1,h2] = yyplot(tce,navl,nmean(P,2),Xlab,Ylab,Ttl,col);
set(ax(1),'xLim',[tStart tEnd]);
title('Low resolution, screened against HR neighbors',style_t{:})
set(ax6,'Ylim',[0 pmax],'Ytick',[0:100:pmax]);
set(ax(2),'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
%
%linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'xy'), 
% export
hepta_figprint(['./figs/composite_ScreeningEffects_' opstring '_Smoothed'])

%% mean composite
p_comp = nmean(proxy_fs,2);
p_std  =  nstd(proxy_fs,0,2);
save(fout)


% FIGURE
fig('global composite'), clf
[ax, h1, h2] = plotyy(tt,gmean_f, tce,p_comp, 'plot');
set(h1,'color',rgb('Silver'),'linewidth',2); 
set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','Ylim',[-.5 .5])
set(get(ax(1),'Ylabel'),'String','Global mean temperature (C)');
set(get(ax(1),'Ylabel'),style_t{:})
%
set(h2,'color',rgb('Red')); 
set(ax(2),'Ycolor',rgb('Red'),'Ylim',[-1 2],'Xtick',[-1:2],'YAxisLocation','left','TickDir','out','YMinorTick','on')
set(get(ax(2),'Ylabel'),'String','Proxy composite')
set(get(ax(2),'Ylabel'),style_t{:})
linkaxes([ax(1),ax(2)],'x'), xlim(ax(1),[0 2010]);
set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
set(ax(2),'XTick',[])
ttl = ['PAGES2k scaled composite vs HadCRUT4, ' int2str(1/f_hr) 'y smoothed' ];
title(ttl,style_t{:});
%axes(ax(2)), hold on
%hs = stairs(tBin,proxyBinAvg,'k-','linewidth',2);
lab{1} = 'HadCRUT4.2 global mean, smoothed';
lab{2} = 'PAGES2k QC composite, smoothed';
%lab{3} = '30-year binnned composite';
hl = legend([h1 h2],lab{:}); set(hl,'FontSize',12,'box','off'); 
%hepta_figprint(['./figs/HadCRUT4_vs_composite_' opstring])

export_fig(['./figs/composite_vsHadCRUT4_' opstring '.pdf'],'-cmyk','-r200','-nocrop','painters')


% %%
% % split along HR/LR criteria
% proxy_hr = proxy_fs(:,1:p_hr); 
% proxy_lr = proxy_fs(:,p_hr+1:p_lr+p_hr); 
% %
% fig('composites'), clf
% ax1 = subplot(2,1,1)
% plot(tce,nmean(proxy_hr,2)), 
% fancyplot_deco('High Resolution (\Delta t <=5y) composite','Year','Composite')
% ax2 = subplot(2,1,2) 
% plot(tce,nmean(proxy_lr,2)), linkaxes([ax1,ax2],'x'), xlim([0 2010]);
% fancyplot_deco('Low Resolution (\Delta t > 5y) composite','Year','Composite')
% export_fig('./figs/split_composite.pdf','-r300','-cmyk')


% p_comp_mat = zeros(size(proxy_fs));
% 
% for k = 1:p_q
%     p_comp_mat(:,k) =  nmean(proxy_fs(:,setdiff([1:p_q],k)),2);
% end
% alph     = 0.05;
% quants   = [alph/2 0.5 (1-alph/2)]; % vector of quantiles
% p_comp_q = quantile(p_comp_mat,quants,2);
% p_std  = nstd(proxy_fs,0,2);
% plot(tce,p_comp_q) % almost no spread at all!!!




% Anaalyze variance in composites
% fig('Variance of composites'), clf
% for k = 1:nType
%     ax_arch(k) = subplot(3,2,k)
%     hs = scatter(n_arch{k},pstd_arch_n{k},36,Graph{u_code(k),1});
%     fancyplot_deco_lm([archType{k} 's'],'# series','Std dev of jackknifed mean (%)')
% end
% hepta_figprint(['./figs/compositeByArchiveJackknife_' opstring],400)








