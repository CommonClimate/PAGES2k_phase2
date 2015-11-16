
% Purpose: investigate how the composite changes when binning sites by
% latitude
% Contributed by: Kevin J Anchukaitis

ts = pages2k.TS_temp(idx_q);

cols = brewermap(5,'Paired');
lat = [ts.geo_latitude_value]'; 

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

fig('By latitude band'), clf
hold on
for i = 1:5   
    sampleSize = sum(~isnan(proxy(:,band(i).index)),2);
    proxys = proxy_r(:,band(i).index);
    insuffientSamples = find(sampleSize<5);
    proxys(insuffientSamples,:) = NaN; 
    subplot(5,1,i); 
    sdx = plot(tce,nanmean(proxys,2),'color',cols(i,:))
    ylim([floor(min(nanmean(proxys,2))) ceil(max(nanmean((proxys),2)))])
    xlim([1 2015])
    title(band(i).title)
    holdYlim = ylim;
    line(tce,nanmean((proxy(:,band(i).index)),2),'color',[0.9 0.9 0.9])
    line(tce,nanmean((proxys),2),'color',cols(i,:))
    ylim(holdYlim)
    set(gca,'xminortick','on','yminortick','on','box','on')
    if i==5
        xlabel('YEAR')
    end
end

set(gcf,'units','normalized','position',[0.39 0.18 0.33 0.75])
%hepta_figprint(['./figs/compositeByLatitude_' options.source],800)
export_fig(['./figs/compositeByLatitude_' smoothString  '.pdf'],'-r200')