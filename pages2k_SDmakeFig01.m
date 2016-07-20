%%  SDmakeFig01
%  Generate Fig 01 of the SD paper
% the figure has to be generated in 2 parts, then glued in Illustrator,
% because Matlab is a retarded language that cannot handle 2 different
% colormaps on the same figure.  Annoying? Mildly so.

% conventions
year = pages2k.year;
style_t = {'FontName','Helvetica','FontSize',16,'FontWeight','bold'};
style_l = {'FontName','Helvetica','FontSize',14};

%%
% DEFINE TIME AXIS AND AVAILABILITY
% avail = ~isnan(proxy);
ny = length(year);
avail = NaN(ny,nr);
p_code = nan(nr,1);
ns    = length(unique({S.dataSetName})); % # of sites
for r = 1:nr
    % define proxy code
    s = strfind(archiveType,lower(archive{r}));
    p_code(r) =  find(~cellfun('isempty', s), 1 ); % number between 1 and na
    avail(ismember(year,[yearMin(r):yearMax(r)]),r)=1;
    % edge color (black except for ice cores)
    if strcmpi(archive{r},'ice core')
        edgec{r} = Graph{5,1};
    else
        edgec{r} = 'k';
    end
end
na = numel(archiveType);
nproxy = zeros(ny,na); pind = zeros(na,1);
for a = 1:na % loop over archive types
    nproxy(:,a) = sum(~isnan(avail(:,p_code == a)),2);
    pind(a) = find(p_code == a,1,'first');
end

% =============
% PLOT SPACETIME COVERAGE
% =============
n1000 = ceil(sum(nproxy(year == 1000,:))/10)*10;
nr    = length(unique({S.dataSetName})); % # of sites
versl = strrep(vers,'_','/');
%%
fig('PAGES 2K'), clf
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [440 144 660 800])
% plot spatial distribution
hmap=axes('Position', [.05 0.67 0.75 0.3]);
m_proj('Robinson','clong',10);
m_coast('patch',[.9 .9 .9]);
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',8,'fontname','Helvetica');
% loop over records
for r = 1:nr
    h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',[7],'linestyle','none');
end
ttl= ['PAGES 2K (Phase 2) as of ' versl ' (',int2str(nr) , ' records from ', int2str(ns), ' sites)'];
%suplabel(ttl)
title('a) Archive types',style_t{:});
% legend
hl = legend(h(pind),archiveType,'location',[.83 .70 .1 .2],style_l{:});
set(hl, 'FontName', 'Helvetica','box','off');
% RESOLUTION
hres = axes('Position', [.05 0.32 0.75 0.3]);
m_proj('Robinson','clong',10);
m_coast('patch',[.9 .9 .9]);
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',8,'fontname','Helvetica');
nc = 11; % number of color contours
scheme = 'RdYlBu'; cx = [-1.2,2.2];
colR = t2c_brewer(log10(resMed),nc,scheme,cx);
% loop over records
for r = 1:nr
    h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',colR(r,:),'linewidth',[1],'MarkerSize',[7],'linestyle','none');
end
% colorbar
ogSize = get(gca, 'Position'); caxis(cx)
hc = colorbar2('vert','Resolution (years)'); dtm = 1/12;
ytl = [dtm,2*dtm,6*dtm, 1, 2,  5, 10, 20, 50, 100]';
yt = log10(ytl);
pos = get(hc,'Position');
set(hc,'YTick',yt,'YLim',cx,'Position',[pos(1)+0.1,pos(2),pos(3),pos(4)]);
set(hc,'YTickLabel',rats(ytl));
set(gca,'Position',ogSize);
title('b) Resolution',style_t{:});

%pause; % manually adjust the size of the figure
hepta_figprint(['./figs/PAGES2K_v' vers '_dbview']);


%% TEMPORAL AVAILABILITY
fig('Temporal avail'), clf
hstack=axes('Position', [0.1 0.1 0.8 0.3]); cmap=cell2mat(Graph(:,1));
colormap(cmap);
ha = area(year,nproxy,'EdgeColor','w'), set(gca,'YAxisLocation','Right','Ylim',[0 620]);
%xlim([1 2000])
fancyplot_deco('','Year (CE)','# proxies');
title('c) Temporal Availability',style_t{:})
% inset
frac=.5;
hstackin=axes('Position', [0.1 0.25 frac*.8 0.1]);
area(year,nproxy,'EdgeColor','w')
axis([1 1000 0 n1000])
set(hstackin,'xtick',[],'box','off','TickDir','out','TickLength',[.02 .02],'YMinorTick','on', 'YGrid','on')
set(hstackin,'YAxisLocation','Right')
set(hstackin,'ytick',[0:50:n1000],'FontName','Helvetica','YColor', [.3 .3 .3])
title('First Millennium',style_l{:})
%pause; % manually adjust the size of the figure

hepta_figprint(['./figs/PAGES2K_v' vers '_tempAvail.eps']);
