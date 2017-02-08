function pages2k_db_synopsis(vers,export)
% function pages2k_db_synopsis(vers,options)
%  INPUT:   vers (version of the database, string)
%           (optional) export. If true, exports figures to pdf. If not,
%           simply prints them to screen. [default =0]
%
% WHAT: Load data from one of the PAGES2k_v..._unpack files
%       defines graphical attributes, and plots a few synopsis figures. 
%
% WHO/WHEN started Dec 5 2014 by Julien Emile-Geay (USC)

% process options
if nargin < 2
    export = 0;
end
year = [];
% load data and packages
addpath(genpath('../utilities'));
load(['../data/PAGES2k_v' vers '_unpack.mat'])
FontName = 'Helvetica';
set(0,'defaultAxesFontName', FontName)
set(0,'defaultTextFontName', FontName)
style_t = {'FontName',FontName,'Fontweight','Bold','Fontsize',16};
style_l = {'FontName',FontName};
figpath = '../figs/';


%%
% extract archive types and assign colors
% =========================================

% Graphical Attributes of each proxy type
archiveType = unique(lower(archive)); na = length(archiveType);
Graph{1,1}=rgb('Gold');        Graph{1,2}= 'h'; % bivalve
Graph{2,1}=rgb('DarkKhaki');   Graph{2,2}= 'h'; % borehole
Graph{3,1}=rgb('DarkOrange');  Graph{3,2}= 'o'; % coral
Graph{4,1}=rgb('Black');       Graph{4,2}= 'p'; % documents
Graph{5,1}=rgb('LightSkyBlue');Graph{5,2}= 'd'; % glacier ice
Graph{6,1}=rgb('DeepSkyBlue'); Graph{6,2}= '*'; % hybrid
Graph{7,1}=rgb('RoyalBlue');   Graph{7,2}= 's'; % lake sediment
Graph{8,1}=rgb('SaddleBrown'); Graph{8,2}= 's'; % marine sediment
Graph{9,1}=rgb('Red');         Graph{9,2}= 'o'; % sclerosponge
Graph{10,1}=rgb('DeepPink');   Graph{10,2}= 'd'; % speleothem
Graph{11,1}=rgb('LimeGreen');  Graph{11,2}= '^'; %'tree'

%%
% DEFINE TIME AXIS AND AVAILABILITY
ny = length(year);
avail = NaN(ny,nr);

p_code = nan(nr,1);
for r = 1:nr
    % define proxy code
    s = strfind(archiveType,lower(archive{r}));
    p_code(r) =  find(~cellfun('isempty', s), 1 ); % number between 1 and na
    avail(ismember(year,[yearMin(r):yearMax(r)]),r)=1;
end

% define edge color (black except for hybrids)
edgec = cell(nr,1);
edgec(1:nr) = {'none'};%{rgb('Black')};
edgec(p_code == 6) = {Graph{6,1}};

nproxy = zeros(ny,na); pind = zeros(na,1);
for a = 1:na % loop over archive types
    nproxy(:,a) = sum(~isnan(avail(:,p_code == a)),2);
    pind(a) = find(p_code == a,1,'first');
    nArch(a)   = sum(p_code == a);
    leg_lbl{a} = [archiveType{a} ' (' int2str(nArch(a)) ')'];
    %disp(leg_lbl{a})
end

% =============
% PLOT SPACETIME COVERAGE
% =============
n1000 = ceil(sum(nproxy(year == 1000,:))/10)*10;
ns    = length(unique({S.dataSetName})); % # of sites
versl = strrep(vers,'_','.'); % legacy
fig('PAGES2K'), clf
orient landscape
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [440   144   896   654])
% plot spatial distribution
hmap=axes('Position', [.05 0.45 0.75 0.5]);
m_proj('Robinson','clong',10);
m_coast('patch',[.9 .9 .9]);
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',4,'fontname',FontName);
% loop over records
for r = 1:nr
    h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',[7],'linestyle','none');
end
text(-2,1.75,['PAGES2k network version ' versl ' (',int2str(nr) , ' records from ', int2str(ns), ' sites)'],style_t{:});
% legend
hl = legend(h(pind),leg_lbl,'location',[.84 .6 .1 .2],style_l{:});
set(hl, 'FontName', FontName,'box','off');
% TEMPORAL AVAILABILITY
hstack=axes('Position', [0.1 0.1 0.8 0.29]);
cmap=cell2mat(Graph(:,1));
colormap(cmap);
area(year,nproxy,'EdgeColor','w'), set(gca,'YAxisLocation','Right');
xlim([1 2000])
fancyplot_deco('','Year (CE)','# proxies',14,'Helvetica');
title('Temporal Availability',style_t{:})
% inset
frac=.5;
hstackin=axes('Position', [0.1 0.2 frac*.8 0.14]);
area(year,nproxy,'EdgeColor','w')
axis([1 1000 0 n1000])
set(hstackin,'xtick',[],'box','off','TickDir','out','TickLength',[.02 .02],'YMinorTick','on', 'YGrid','on')
set(hstackin,'YAxisLocation','Right')
set(hstackin,'FontName',FontName,'YColor', [.3 .3 .3])
title('First Millennium',style_l{:})

%pause; % manually adjust the size of the figure
if export
   %figpath = [strrep(userpath,':','/') 'tempfigs/'];
   fname = ['PAGES2k_' vers '_spacetime.pdf'];
   export_fig([figpath fname],'-r200','-cmyk','-painters');
   %hepta_figprint(['../../figs/synopsis/PAGES2K_phase2_' vers '_spacetime']);
   %saveFigure([figpath fname '_saveFigure.pdf']); %'painters', true
   %plot2svg([figpath fname '.svg']);
   %print([figpath fname],'-depsc','-r300','-cmyk')
end


%% polar projections

fig('PAGES2K projections'), clf
orient landscape
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [440   144   896   654])
p = panel();
p.pack('v', {1/2 []})
p(1).pack(1, 2);
rs = 30;
% Arctic records
p(1,1,1).select()
m_proj('stereographic','lat',90,'long',0,'radius',rs);
m_coast('patch',[.9 .9 .9]);
m_grid('xtick',6,'ytick',7,'xlabeldir','out','fontsize',10,'fontname',FontName,'Yticklabel',[]);
recs = find(p_lat >= 90-rs);
for r = recs
    ha(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',8,'linestyle','none');
end
% Antarctic records
p(1,1,2).select()
m_proj('stereographic','lat',-90,'long',0,'radius',rs);
m_coast('patch',[.9 .9 .9]);
m_grid('xtick',6,'ytick',7,'xlabeldir','out','fontsize',10,'fontname',FontName,'XaxisLocation','top','Yticklabel',[]);
% loop over records
recs = find(p_lat <= -90+rs);
for r = recs
    ha(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',8,'linestyle','none');
end
%  tropical records
p(2).select()
Slongs=[-100 43;-75 20; 20 145;43 100;145 295;100 295];
Slats= [  0  62;-62  0;-62   0; 0  62;-62   0;  0  62];
for l=1:6,
    m_proj('mollweide','long',Slongs(l,:),'lat',Slats(l,:));
    m_grid('fontsize',4,'xtick',[-180:30:360], 'xticklabels',[],...
        'ytick',[-80:20:80],'yticklabels',[],'linest',':','fontname',FontName);
    m_coast('patch',[.9 .9 .9]);
    plons = p_lon;
    if l>4
        plons(p_lon<0) = plons(p_lon<0)+360;
    end
    % loop over records
    recs = find(plons >= Slongs(l,1) & plons <Slongs(l,2));
    for r = recs
        h(r) = m_line(plons(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',[7],'linestyle','none');
    end
end;
set(gca,'xlimmode','auto','ylimmode','auto');
% titles
text(-1.3,2.8,'a) Arctic sites',style_t{:});
text(+1.1,2.8,'b) Antarctic sites',style_t{:});
text(0.5,0.9,'c) Non-polar sites',style_t{:});

if export
   fname = ['PAGES2k_' vers '_polarOrNot.pdf'];
   export_fig([figpath fname],'-r200','-cmyk','-painters');
end


%%  plot resolution of non-tree records
fig('Resolution'),clf
res_notree = resAvg(p_code < 11);
x = [1/12 1/6 1/4 1 2 5 10 20 50 100]'; hist(res_notree,x)
hp = findobj(gca,'Type','patch'); set(hp,'FaceColor',rgb('Silver'))
set(gca,'XScale','log');
fancyplot_deco('Average Resolution of non tree-ring proxies','Resolution (years)','Bin count')
if export
   fname = ['PAGES2k_' vers '_resolution.pdf'];
   export_fig([figpath fname],'-r300','-cmyk','-painters','-nocrop');
end

% Basic Stats 
fig('Span'),clf
Span = yearMax - yearMin;
Span(Span>2000) = 2000;
bins = [50 75 100:100:900 1000:500:2000];
hist(Span,bins)

save(['../data/PAGES2k_v' vers '_unpack.mat'],'-append','p_code','ny','Graph','archiveType','archive','na','recordNames','avail','edgec','-v6')
