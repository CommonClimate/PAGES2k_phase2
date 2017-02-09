
% Purpose: investigate how the composite changes when binning sites by
% latitude
% Contributed by: Kevin J Anchukaitis
binStep = 50;
binEdges=flipud((2000:-binStep:0)');
binYear=mean([binEdges(2:end) binEdges(1:end-1)],2);
plotEdges=reshape([binEdges(1:end-1) binEdges(2:end)]',[],1);
tb = binEdges(2:end);
% analysis parameters
n_tresh = 10; nboot = 500; % # of bootstrap samples.
% plotting parameters
style_l = style_t; style_l{4} = 12;  xlims = [0 2000];  K = 1.1; % expansion factor for tick marks


ts = pages2k.S(idx_q);

cols = brewermap(5,'Set1');
lat = [ts.geo_meanLat]'; 

% -90 to -60
band(5).index = find(lat>=-90 & lat<=-60);
band(5).title = '60S to 90S';
band(4).index = find(lat>-60 & lat<=-30);
band(4).title = '30S to 60S';
band(3).index = find(lat>-30 & lat<30);
band(3).title = '30S to 30N';
band(2).index = find(lat>=30 & lat<60);
band(2).title = '30N to 60N';
band(1).index = find(lat>=60 & lat<=90);
band(1).title = '60N to 90N';

nmin = 5 ; %minimum number of records
ylims       = [-1 1.5]; dy = 0.5;
yticks      = ylims(1):dy:ylims(2);  % define tick marks

fig('By latitude band'), clf
hold on
for c = 1:5   
    subplot(5,1,c)
    proxy_b   = bin_x(tce',proxy_q(:,band(c).index),binEdges);
    pcom      = nmean(proxy_b,2);
    p_boot    = bootstrp(nboot,@nmean,proxy_b');
    ci_arch   = quantile(p_boot',[0.025 0.975],2);
    p_arch(c) = size(proxy_b,2);    
    n_arch    = sum(~isnan(proxy_b),2); 
    
    n_archPlot = reshape([n_arch n_arch]',[],1);
    pcom_archPlot = reshape([pcom pcom]',[],1);
    ciLo        = reshape([ci_arch(:,1) ci_arch(:,1)]',[],1);
    ciHi        = reshape([ci_arch(:,2) ci_arch(:,2)]',[],1);
    
    % plotting threshold
    thresh = (n_archPlot  >= nmin);
    col = cols(c,:);
    % plot solid line otherwise, and number of proxies
    [ax,h1,h2]  = plotyy(tb,n_arch,plotEdges,pcom_archPlot,@bar,@line); hold on
    set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
    set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right');
    set(ax(1),'yLim',[0 ceil(p_arch(c)/10)*10]);
    set(get(ax(1),'Ylabel'),'String','# records'), 
    set(get(ax(1),'Ylabel'),style_l{:}), set(ax(2),'yLim',ylims)
    %set(h2,'visible','off'); %make invisible for now;
    set(h2,'color',col,'linewidth',2,'linestyle',':');
    set(ax(2),'Ycolor',col,'YAxisLocation','left','Ytick',yticks)
    set(get(ax(2),'Ylabel'),'String','Composite');
    set(get(ax(2),'Ylabel'),style_l{:})
    set(gcf,'CurrentAxes',ax(2))
    set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    
    % plot bootstrap CI
    wide = (~isnan(ciLo) & ciHi-ciLo>range(ylims)/50.0);  %
    hci = area_fill(plotEdges(wide)',ciLo(wide)',ciHi(wide)',col,col,0.2);
    %plot transparent line for entire period
    hp = line(plotEdges(thresh),pcom_archPlot(thresh),'color',col,'linewidth',2,'linestyle','-','Parent',ax(2));
    % plot solid line when enough data are present
    uistack(h2,'top')
    % more cosmetics
    linkaxes([ax(1),ax(2)],'x'); set(ax(1),'xLim',xlims);
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XTick',[]) 
    ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
    title(band(c).title,style_l{:});
    
    
%     sampleSize = sum(~isnan(proxy_q(:,band(i).index)),2);
%     
%     insuffientSamples = find(sampleSize<5);
%     proxy_i = proxy_q(:,band(i).index);
%     proxy_i(insuffientSamples,:) = NaN; 
%     % bin it
%     bin_i = bin_x(tce',proxy_i,binEdges);
%     comp_i       = nmean(bin_i,2);
%     pcomp_i   = reshape([comp_i comp_i]',[],1);
%     bin_t = bin_x(tce',proxy_t,binEdges);
%     comp_t       = nmean(bin_t,2);
%     pcomp_t   = reshape([comp_t comp_t]',[],1);
%     subplot(5,1,i); 
%     plot(plotEdges,pcomp_i,'color',cols(i,:),'LineWidth',2); hold on    
%     ylim([floor(min(comp_i)) ceil(max(comp_i))])
%     xlim([1 2015])
%     title(band(i).title)
%     holdYlim = ylim;
%     line(plotEdges,pcomp_t,'color',[0.8 0.8 0.8],'LineWidth',2);
%     line(plotEdges,pcomp_i,'color',cols(i,:),'LineWidth',2); 
%     ylim(holdYlim)
%     set(gca,'xminortick','on','yminortick','on','box','off','tickdir','out','ygrid','on')
    if c==5
        xlabel('Year CE',style_l{:}) 
    else
        set(ax(1),'XTickLabel',[])
    end
end
suptitle(['Composite by latitude, n_min = ',int2str(nmin)],'Fontweight','Bold')

%set(gcf,'units','normalized','position',[0.39 0.18 0.33 0.75])
% export
froot = ['./figs/' opstring '_compositeByLatitude'];
hepta_figprint(froot,400)  %export to PDF as one needs to adjust transparency manually in Illustrator.
eps2pdfMac([froot '.eps'])

