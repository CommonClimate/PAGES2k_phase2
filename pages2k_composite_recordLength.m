
tMin = tEnd-[250 500 1000 1999];
yearMinQ = yearMin(idx_q);
cols = brewermap(length(tMin),'Paired');

fig('Sensitivity to record length'), clf
set(gcf,'Position',[440   270   852   628])
pmax = 300; % scale for # proxies
ylims = [-1.5 2]; thresh = 10;
hold on
for i = 1:4
    ax1 = subplot(2,2,i);
    Xlab = ''; Ylab = {'# records','Composite'};
    P = proxy_r(:,yearMinQ<=tMin(i)); ps = size(P,2); navl = sum(~isnan(P),2);
    Pc = nmean(P,2); Pn = Pc; Pn(navl<5) = NaN;
    % plot solid line otherwise, and number of proxies
    [ax,h1,h2]  = plotyy(tce,navl,tce,Pc,@bar,@plot); hold on
    set(h1,'edgecolor',rgb('Gainsboro')); % set(h1,'edgealpha',0.5);
    set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','TickDir','out','YMinorTick','on')
    set(get(ax(1),'Ylabel'),'String','# records')
    set(get(ax(1),'Ylabel'),style_l{:})
    set(h2,'visible','off'); %make invisible for now;
    set(h2,'color',rgb('Blue'),'linewidth',1);
    set(ax(2),'Ycolor',rgb('Blue'),'YAxisLocation','left')
    set(get(ax(2),'Ylabel'),'String','Comp.');
    set(get(ax(2),'Ylabel'),style_l{:})
    set(gcf,'CurrentAxes',ax(2))
    set(h2,'XData',tce(navl>=thresh));
    set(h2,'YData',Pc(navl>=thresh));
    refreshdata(h2);
    set(h2,'visible','on');
    %plot transparent line for entire period
    hp = patchline(tce,Pc,'edgecolor',rgb('Blue'),'linewidth',2,'edgealpha',0.3);
    % plot solid line when enough data are present
    uistack(h2,'top')
    % more cosmetics
    linkaxes([ax(1),ax(2)],'x'); set(ax(1),'xLim',xlims,'ylim',[0 pmax]);
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XTick',[])
    ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
    title([int2str(tMin(i)) 'CE and before' ],style_l{:});
end

set(gcf,'units','normalized','position',[0.1 0.1 0.6 0.4])
%hepta_figprint(['./figs/compositeByLatitude_' options.source],800)
hepta_figprint(['./figs/compositeByStartDate_' opstring '_' smoothString])
