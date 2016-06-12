
% Purpose: investigate how the composite changes when binning sites by
% latitude
% Contributed by: Kevin J Anchukaitis

ts = pages2k.S(idx_q);

cols = brewermap(5,'Paired');
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
fig('By latitude band'), clf
hold on
for i = 1:5   
    sampleSize = sum(~isnan(proxy(:,band(i).index)),2);
    proxy_t = proxy_sgn(:,band(i).index);
    insuffientSamples = find(sampleSize<5);
    proxy_i = proxy_t;
    proxy_i(insuffientSamples,:) = NaN; 
    % bin it
    bin_i = bin_x(tce',proxy_i,binEdges);
    comp_i       = nmean(bin_i,2);
    pcomp_i   = reshape([comp_i comp_i]',[],1);
    bin_t = bin_x(tce',proxy_t,binEdges);
    comp_t       = nmean(bin_t,2);
    pcomp_t   = reshape([comp_t comp_t]',[],1);
    subplot(5,1,i); 
    plot(plotEdges,pcomp_i,'color',cols(i,:),'LineWidth',2); hold on    
    ylim([floor(min(comp_i)) ceil(max(comp_i))])
    xlim([1 2015])
    title(band(i).title)
    holdYlim = ylim;
    line(plotEdges,pcomp_t,'color',[0.8 0.8 0.8],'LineWidth',2);
    line(plotEdges,pcomp_i,'color',cols(i,:),'LineWidth',2); 
    ylim(holdYlim)
    set(gca,'xminortick','on','yminortick','on','box','off','tickdir','out','ygrid','on')
    if i==5
        xlabel('YEAR')
    end
end
suptitle(['Composite by latitude, n_min = ',int2str(nmin)])

set(gcf,'units','normalized','position',[0.39 0.18 0.33 0.75])
% export
froot = ['./figs/' opstring '_compositeByLatitude.pdf'];
export_fig(froot,'-r300','-cmyk')

%hepta_figprint(['./figs/compositeByLatitude_' smoothString],400)
%export_fig(['./figs/compositeByLatitude_' opstring '_' smoothString  '.pdf'],'-r200')
