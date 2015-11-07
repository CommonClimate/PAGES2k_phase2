%% Composite plus regression

% align proxies and temperature matrices.
nproxy = zeros(nce,1);
for j=1:nce
	nz = ~isnan(proxy(j,:));
	nproxy(j) = sum(nz);
end
iRecon = (nproxy>=10);  % 1 to 2010.

tr  = tce(iRecon); nt = length(tr);
proxy_r = proxy_fs(iRecon,:);
ind_t = find(hadcrut4.t <= max(tr));

inst = double(gmean_f(ind_t)); ti = hadcrut4.t(ind_t); ni = length(ti);
calib = ismember(tr,ti);  

% missingness patterns 
[np, kavl, kmis, prows, mp, iptrn] = missingness_patterns(proxy_r);
pmin = 1; pmax = 2;  ns = 50; 
alph = 0.05; quants   = [alph/2 0.5 (1-alph/2)]; % vector of quantiles
gmean_cpr = zeros(nt,3);
lab{1} = '95% prediction interval';
lab{2} = 'HadCRUT4.2 global mean';
lab{3} = 'PAGES2k reconstruction';
RE = zeros(nt,1); CE = RE; R2=RE;
yp_noise = NaN(nt,ns);

for j = 1:np
    disp(['Reconstructing for pattern ' int2str(j) '/' int2str(np)])
    p_nest{j} = nmean(proxy_r(:,kavl{j}),2); % proxy composite on nest j
    n_nest{j} = numel(kavl{j});
    % calibrate to temperature
    mdl = fitlm(p_nest{j}(calib),inst);
    % predict
    yp = predict(mdl,p_nest{j}); 
    % verification statistics
    pattern = (iptrn==j); lp = sum(pattern);
    [REn,CEn,R2n]=verif_stats(inst,yp(calib),1:ni,1:ni);
    RE(pattern,:) = repmat(REn,lp,1);
    CE(pattern,:) = repmat(CEn,lp,1);
    R2(pattern,:) = repmat(R2n,lp,1);
    
    %fig('Residuals'), clf
    %plotResiduals(mdl,'lagged')% looks AR(1)
    res = mdl.Residuals.Raw;
    % fit ar model to residuals
    [w,A,C,SBC,FPE,th] = arfit(res,1,pmax);
    
    nj = numel(prows{j});
    noise  = zeros(nj,ns);
    disp('   Generating noise')
    for s = 1:ns
        noise(:,s) = arsim(w,A,C,nj);
    end
    yp_noise(prows{j},:) = repmat(yp(prows{j}),[1 ns]) + noise;
    
    % update the global mean reconstruction with newer nest

    
%     fig('gmean cpr'),clf
%     ha = area_fill(tr,yp_quants{j}(:,3)',yp_quants{j}(:,1)',rgb('Silver'),rgb('Gray')); hold on
%     h1 = plot(tt,gmean_f,'color',rgb('Black'),'linewidth',2);
%     h2 = plot(tr,yp_quants{j}(:,2),'color',rgb('Red'));
%     axis([0 2010 -.8 .8])    
%     hl = legend([ha h1 h2],lab{:}); set(hl,'FontName','Palatino','FontSize',12,'box','off');
%     ttl = ['Reconstruction with PAGES2k pattern ' int2str(j) '/' int2str(np),', ',int2str(1/f_hr) 'y smoothed' ];
%     fancyplot_deco(ttl,'year CE','Global Mean Temperature (C)')
%     line(tr(prows{j}),0.4*ones(mp(j),1),'Marker','o','Markersize',6,'MarkerEdgeColor','none','MarkerFaceColor',rgb('DimGray'))
%    
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperPositionMode','manual')  
%     fn = ['./figs/pages2k_gmean_cpr_nest',sprintf('%03d',j),'.png'];
%     disp('   Exporting figure')
%     export_fig(fn,'-r150','-cmyk','-nocrop')
    % save frame!  
end

save(fout)


yp_quants = quantile(yp_noise,quants,2);

for q = 1:3
    gmean_cpr(:,q) = hepta_smooth(yp_quants(:,q),f_hr);
    %gmean_cpr(prows{j},q) = gmean_qs(prows{j});
end


fig('gmean cpr splice'),clf
[ax, h1, h2] = plotyy(tce,nproxy, tr,gmean_cpr(:,2),'bar','plot'); hold on
set(h1,'Edgecolor',rgb('Gainsboro')); 
set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','Ylim',[0 700],'Ytick',[0:100:700],'TickDir','out','Box','Off')
set(get(ax(1),'Ylabel'),'String','# records');
set(get(ax(1),'Ylabel'),style_t{:})
%
set(h2,'color',rgb('Red'),'linewidth',2); 
set(ax(2),'Ycolor',rgb('Red'),'Ylim',[-.8 .6],'YAxisLocation','left','TickDir','out')
set(ax(2),'YMinorTick','on','Xtick',[],'Ytick',[-.8:.2:.6])
set(get(ax(2),'Ylabel'),'String','Global Mean Temperature')
set(get(ax(2),'Ylabel'),style_t{:})
linkaxes([ax(1),ax(2)],'x'), xlim(ax(1),[0 2010]);
set(gcf, 'CurrentAxes', ax(2));
h0 = line(tt,gmean_f,'color',rgb('Black'),'linewidth',2);
ha = area_fill(tr,gmean_cpr(:,3)',gmean_cpr(:,1)',rgb('Red'),rgb('Red'),0.2,0.3); 
uistack(h0,'top')
%
set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
set(ax(2),'XTick',[])
%
hl = legend([ha h0 h2],lab{:}); set(hl,'FontName','Helvetica','FontSize',12,'box','off','location','northwest');
ttl = ['PAGES2k spliced reconstruction, smoothed'];
fancyplot_deco(ttl,'year CE','Global Mean Temperature (C)',16,'Helvetica')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperPositionMode','manual')
hepta_figprint(['./figs/pages2k_cpr_splice_' opstring])

save(fout)

%% k-fold cross validation
% define folds 
pcom_r = nmean(proxy_r(ismember(tr,ti),:),2);
gmean_kcv = nan(size(inst));
Kcv = 8;
[in, out, nin] = kcv_indices([1:ni], Kcv, 'blocks'); 
REns = zeros(ni,1); lev = 100*(1-alph); % percentile
CEns = REns; R2ns = REns; MSEns = REns; 
REk = zeros(Kcv,1); tmid = REk;
%cols = cubehelix(Kcv,[0.5,-1.5,1,1],[0.2,0.8]);
cols = brewermap(Kcv,'Set1');
for k = 1:Kcv
    % persistence-based prediction
    [REi,CEi,R2i,MSEi] = crossVal_ar1(double(inst),in{k},out{k},1000,lev);
    nout  = numel(out{k});
    REns(out{k}) = REi; CEns(out{k}) = CEi; R2ns(out{k}) = R2i; MSEns(out{k}) = MSEi;
    % proxy-based prediction
    mdl = fitlm(pcom_r(in{k}),inst(in{k}));
    % predict
    gmean_kcv(out{k}) = predict(mdl,pcom_r(out{k}));
end
for k = 1:Kcv
    [REk(k),CEk(k),R2k(k),MSEk(k)]=verif_stats(inst,gmean_kcv,in{k},out{k});
    tmid(k) = median(ti(out{k}));
end

fig('Validation'),clf
hold on
clear ax h1 h2
[ax, h1, h2] = plotyy(tmid,REk, ti,inst,'bar','plot'); hold on
set(h1,'Facecolor',rgb('Gainsboro'),'edgecolor','none'); 
set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','Ylim',[-.4 1])
set(ax(1),'TickDir','out','Box','Off','Ytick',[-.4:0.2:1],'YMinorTick','on')
set(get(ax(1),'Ylabel'),'String','RE');
set(get(ax(1),'Ylabel'),style_t{:})
%
set(h2,'color',rgb('Black'),'linewidth',2); 
set(ax(2),'Ycolor',rgb('Black'),'Ylim',[-.5 .5],'YAxisLocation','left','TickDir','out')
set(ax(2),'YMinorTick','on','Ytick',[-.5:.1:.5])
set(get(ax(2),'Ylabel'),'String','Global Mean Temperature')
set(get(ax(2),'Ylabel'),style_t{:})
linkaxes([ax(1),ax(2)],'x'), xlim(ax(1),[1850 2010]);
%
set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on')
set(ax(2),'XColor' , [.3 .3 .3], 'LineWidth', 1);
set(ax(1),'XTick',[])
%
for k = 1:Kcv
    hREcrit=line(ti(out{k}),REns(out{k}),'color',cols(k,:),'linewidth',2,'linestyle','--');
end
axes(ax(2))
for k = 1:Kcv
    hrec = line(ti(out{k}),gmean_kcv(out{k}),'color',cols(k,:),'linewidth',4);
end
hold off
lab{1} = 'HadCRUT4.2';
lab{2} = 'Proxy-predicted';
lab{3} = 'Reduction of error (RE)';
lab{4} = 'Critical RE (AR(1) benchmark)';
hl = legend([h2 hrec h1 hREcrit],lab{:}); set(hl,'FontName','Palatino','FontSize',12,'box','off','location','SouthEast');
%
fancyplot_deco('Instrumental vs proxy-predicted GMT','year','Global Mean Temperature (C)')
set(gca,'Ygrid','off')
hepta_figprint(['./figs/pages2k_cpr_kcv_' opstring])
