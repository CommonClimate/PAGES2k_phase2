clear all;
addpath(genpath('./utilities'))
network = 'M08';  % choices: M08, PAGES 2013, PAGES 2016. 
vers = '1_13_1';  % version of the PAGES2k database to be use

% fonts
FontName = 'Helvetica';
set(0,'defaultAxesFontName', FontName)
set(0,'defaultTextFontName', FontName)
style_t = {'FontName',FontName,'Fontweight','Bold','Fontsize',16};
style_l = {'FontName',FontName};

% define appearance for every proxy class
icon{1,1}=rgb('Gold');        icon{1,2}= 'h';  icon{1,3} = 'bivalve';
icon{2,1}=rgb('DarkKhaki');   icon{2,2}= 'h';  icon{2,3} = 'borehole';
icon{3,1}=rgb('DarkOrange');  icon{3,2}= 'o';  icon{3,3} = 'coral';
icon{4,1}=rgb('DimGray');     icon{4,2}= 'p';  icon{4,3} = 'document';
icon{5,1}=rgb('LightSkyBlue');icon{5,2}= 'd';  icon{5,3} = 'glacier ice';
icon{6,1}=rgb('DeepSkyBlue'); icon{6,2}= '*';  icon{6,3} = 'hybrid';
icon{7,1}=rgb('RoyalBlue');   icon{7,2}= 's';  icon{7,3} = 'lake sediment';
icon{8,1}=rgb('SaddleBrown'); icon{8,2}= 's';  icon{8,3} = 'marine sediment';
icon{9,1}=rgb('Red');         icon{9,2}= 'o';  icon{9,3} = 'sclerosponge';
icon{10,1}=rgb('DeepPink');   icon{10,2}= 'd'; icon{10,3} = 'speleothem';
icon{11,1}=rgb('LimeGreen');  icon{11,2}= '^'; icon{11,3} = 'tree ring';
na = length(icon);

switch network
    case('M08')
        load('./data/mann08_proxymap');
        p_code(round(ptype/1000) == 9) = 11; % Tree Ring
        p_code(round(ptype/1000) == 8) = 5; % Ice
        p_code(round(ptype/1000) == 7) = 3; % Coral
        p_code(round(ptype/100) == 75) = 11; % Tree
        p_code(round(ptype/1000) == 6) = 10; % Speleo
        p_code(round(ptype/1000) == 5) = 4; % Doc
        p_code(round(ptype/1000) == 4) = 7; % Sediment
        p_code(round(ptype/1000) == 3 | round(ptype/1000) == 2) = 6; %Hybrid ("Composite" in Mann08)
        prox = proxy_nl;
        p_lon = plon;
        p_lat = plat;
        year = time; 
        [ny,nr]   = size(prox);             
        avail = nan(ny,nr); avail(~isnan(prox)) = 1;

    case{'PAGES 2013'}
        P1  = load('./data/pages2k_hadcrut4_data_v1');
        proxy = P1.proxy;
        keep      = setdiff(1:550,proxy.duplic);    % keep      = proxy.keep;
        prox      = proxy.raw(:,keep);  % no duplicate NA tree rings
        ptype     = proxy.type(keep);
        p_lon      = proxy.lon(keep);
        p_lat      = proxy.lat(keep);
        year      = proxy.t;
        [ny,nr]   = size(prox);
        
        p_code  = nan(nr,1); %na = 8; % 8 proxy types for this network
        p_code(strcmp('tree ring',ptype)) = 11;
        p_code(strcmp('instrumental',ptype)) = 6;
        p_code(strcmp('ice core',ptype)) = 5;
        p_code(strcmp('coral',ptype)) = 3;
        p_code(strcmp('speleothem',ptype)) = 10;
        p_code(strcmp('documentary',ptype)) = 4;
        p_code(strcmp('lake sediment',ptype)|strcmp('marine sediment',ptype)) = 7;
        p_code(strcmp('historic',ptype)) = 4;
        %
        avail = nan(size(prox)); avail(~isnan(prox)) = 1;
        
    case{'PAGES 2016'}
        load(['./data/pages2kTSv' vers '_unpack.mat']);
        [ny,nr]   = size(proxy_ann);       
end

%% compute essential variables


nproxy = zeros(ny,na); pind = nan(na,1);
for a = 1:na % loop over archive types
    nproxy(:,a) = sum(~isnan(avail(:,p_code == a)),2);
    if sum(p_code == a)>0
        pind(a) = find(p_code == a,1,'first');
        nArch(a)   = sum(p_code == a);
        lbl{a} = [icon{a,3} ' (' int2str(nArch(a)) ')'];
    end
end

% define edge color (black except for ice cores)
edgec = cell(nr,1);
edgec(1:nr) = {rgb('black')};
edgec(p_code == 6) = {icon{6,1}};

% define graphical parameters
axlim   = [0 2000 0 1200];
inaxlim = [0 1000 0 90];
ytick   = [25 50 75];

% =============
% PLOT SPACETIME COVERAGE
% =============
n1000 = ceil(sum(nproxy(year == 1000,:))/10)*10;
fig('Fixed Size'), clf
orient landscape
set(gcf,'PaperPositionMode','auto')
%set(gcf, 'Position', [440   144   896   654])
set(gcf, 'Position', [310   289   663   516])
% plot spatial distribution
hmap=axes('Position', [.05 0.45 0.75 0.5]);
m_proj('Robinson','clong',10);
m_coast('patch',rgb('WhiteSmoke'));
m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',4);
% loop over records
for r = 1:nr
    h(r) = m_line(p_lon(r),p_lat(r),'marker',icon{p_code(r),2},'MarkerEdgeColor',edgec{r},'MarkerFaceColor',icon{p_code(r),1},'linewidth',[.5],'MarkerSize',[6],'linestyle','none');
end
text(-.45,1.75,['Network: ' network ', ' int2str(nr) ' records'],'FontWeight','bold','FontName', FontName,'Fontsize',16);
% legend
hl = legend(h(pind(~isnan(pind))),lbl{find(~isnan(pind))},'location',[.84 .6 .1 .2]);
set(hl, 'FontName', FontName,'box','off','FontSize',12);
% TEMPORAL AVAILABILITY
hstack=axes('Position', [0.1 0.1 0.7 0.29]);
cmap=cell2mat(icon(:,1));
colormap(cmap);
area(year,nproxy,'EdgeColor','w'), set(gca,'YAxisLocation','Right');
xlim([0 2000])
%axis(axlim), %set(gca, 'Yscale', 'log','BaseValue', 1e0);
fancyplot_deco('','Year (CE)','# proxies',14,FontName);
title('Temporal Availability','FontName', FontName,'Fontsize',14)

%inset
frac=.5;
hstackin=axes('Position', [0.1 0.2 frac*.7 0.14]);
area(year,nproxy,'EdgeColor','w')
axis([1 1000 0 n1000])
set(hstackin,'xtick',[],'box','off','TickDir','out','TickLength',[.02 .02],'YMinorTick','on', 'YGrid','on')
set(hstackin,'YAxisLocation','Right')
set(hstackin,'ytick',[0:50:n1000],'FontName','Helvetica','YColor', [.3 .3 .3])
title('First Millennium','FontName','Helvetica','Fontsize',12,'Fontweight','bold')


%% export to PDF
export_fig(['./figs/' network '_dbviewNarrow.pdf'],'-r200','-cmyk','-painters');





% % Instrumental availability only
% fig('temporal availability'); set(gcf,'position',[233   394   970   324]); clf;
% hb = bar(ptime,nproxy,'stacked');
% colormap(proxycolor)
% axis([1850 2011 0 np])
% fancyplot_deco('Proxy availabity over the istrumental period (1850 - 2011) ','Time (Year A.D.)','# proxies',20);
% [l2,OBJH,OUTH,OUTM]=legend(hb(ind~=0),icon{p_code(ind(ind~=0)),3});
% set(l2,style_l{:});
% set(l2,'position',[0.8658    0.3648    0.1169    0.5485])

%export_fig([upper(network) '_avail_' num2str(np) '_instrumental.pdf'],'-cmyk','-r300')
