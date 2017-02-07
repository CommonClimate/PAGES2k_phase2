%%  SDmakeFig01
%  Generate Fig 01 of the SD paper

% conventions
year = pages2k.year;
FontName = 'Helvetica';
set(0,'defaultAxesFontName', FontName)
set(0,'defaultTextFontName', FontName)
style_t = {'FontName',FontName,'Fontweight','Bold','Fontsize',14};
style_l = {'FontName',FontName,'Fontweight','Bold','Fontsize',12};
%%
% DEFINE TIME AXIS AND AVAILABILITY
% avail = ~isnan(proxy);
ny = length(year); nr = length(pages2k.S);
ns    = length(unique({S.dataSetName})); % # of sites
na = numel(archiveType);
nproxy = zeros(ny,na); pind = zeros(na,1);
avail = pages2k.avail;
for a = 1:na % loop over archive types
    nproxy(:,a) = sum(~isnan(avail(:,p_code == a)),2);
    pind(a) = find(p_code == a,1,'first');
    nArch(a)   = sum(p_code == a);
    leg_lbl{a} = [archiveType{a} ' (' int2str(nArch(a)) ')'];
end

% define edge color (black except for hybrids)
edgec = cell(nr,1);
edgec(1:nr) = {'none'};%{rgb('Black')};
edgec(p_code == 6) = {Graph{6,1}};

% =============
% PLOT SPACETIME COVERAGE
% =============
n1000 = ceil(sum(nproxy(year == 1000,:))/10)*10;
versl = strrep(vers,'_','/');
%% define figure geometry
fig('PAGES2k'), clf
set(gcf, 'Units','centimeters', 'Position',[0 0 20 27])
set(gcf, 'PaperPositionMode','auto')
p = panel(); p.pack([2 96 2]); p.margin = [10 10 10 10]
p(2).pack({2/5 2/5 1/5},{4/5 1/5});
p.fontsize = 12; p.fontname = 'Helvetica';

% title
p(1).select(), set(gca,'Visible','off')
ttl= ['PAGES2k 2.0.0 (',int2str(nr) , ' records from ', int2str(ns), ' sites)'];
text(0.2,0.5,ttl,'Fontweight','bold','Fontsize',16)

% plot spatial distribution
%ax1=axes('Position', [.05 0.67 0.75 0.3]);
p(2,1,1).select()
m_proj('Robinson','clong',10);
m_coast('patch',rgb('WhiteSmoke'));
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',8,'fontname','Helvetica');
% loop over records
for r = 1:nr
    h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',0.5,'MarkerSize',7,'linestyle','none');
end
title('a) Archive types',style_t{:});

% legend
hl = legend(h(pind),leg_lbl,'location','EastOutside',style_l{:});
set(hl, 'FontName', 'Helvetica','box','off');
% RESOLUTION
%ax2 = axes('Position', [.05 0.32 0.75 0.3]);
p(2,2,1).select()
m_proj('Robinson','clong',10);
m_coast('patch',rgb('WhiteSmoke'));
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',8,'fontname','Helvetica');
nc = 11; % number of color contours
scheme = 'RdYlBu'; cx = [-1.2,2.2];
colR = t2c_brewer(log10(resMed),nc,scheme,cx);
% loop over records
for r = 1:nr
    h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',colR(r,:),'linewidth',0.5,'MarkerSize',7,'linestyle',':');
end
% colorbar
ogSize = get(gca, 'Position'); caxis(cx)
hc = colorbar2('vert','Resolution (years)'); dtm = 1/12;
ytl = [dtm,2*dtm,6*dtm, 1, 2,  5, 10, 20, 50, 100]';
yt = log10(ytl);
pos = get(hc,'Position');
set(hc,'YTick',yt,'YLim',cx,'Position',[pos(1)+0.1,pos(2),pos(3),pos(4)]);
set(hc,'YTickLabel',rats(ytl),'FontName','Helvetica');
set(gca,'Position',ogSize);
title('b) Resolution',style_t{:});
hold on
%pause; % manually adjust the size of the figure
%hepta_figprint(['./figs/PAGES2k_v' vers '_dbview']);

%% TEMPORAL AVAILABILITY
p(2,3,1).select()
plot(p(2,3,1).axis,tce,nmean(proxy_q,2))
area(p(2,3,1).axis,year,nproxy,'EdgeColor','w');  colormap(gca,cell2mat(Graph(:,1)));
set(gca,'YAxisLocation','Right','Ylim',[0 620],'XLim',[0 2000],'FontName','Helvetica');
fancyplot_deco('','Year (CE)','# proxies',12,'Helvetica');
title('c) Temporal Availability',style_t{:})
% inset
pos = p(2,3,1).axis.Position;
frac=.5;
hstackin=axes('Position', [pos(1) pos(1)+0.1 pos(3)*frac pos(4)*frac]);
area(year,nproxy,'EdgeColor','w'), colormap(hstackin,cell2mat(Graph(:,1)));
axis([1 1000 0 200])
set(hstackin,'xtick',[],'box','off','TickDir','out','TickLength',[.02 .02],'YMinorTick','on', 'YGrid','on')
set(hstackin,'YAxisLocation','Right','YTick',[0:50:200],'FontSize',10)
set(hstackin,'FontName','Helvetica','YColor', [.3 .3 .3])
title('First Millennium','Fontweight','bold')
display('Adjust figure size and press Enter')
pause; 
export_fig(['./figs/PAGES2k_v' vers '_dbview.pdf'],'-r300','-nocrop','-cmyk','-painters');

%export_fig(['./figs/PAGES2k_v' vers '_tempAvail.pdf'],'-r300','-cmyk','-painters','-transparent');
%hepta_figprint(['./figs/PAGES2k_v' vers '_tempAvail.eps']);
