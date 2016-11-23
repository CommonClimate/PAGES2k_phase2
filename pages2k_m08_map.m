clear all; 

load JEG_graphics

icon{1,1}=Lawn_green;       icon{1,2}= '^'; icon{1,3} = 'Tree Ring';
icon{2,1}=rgb('DarkKhaki'); icon{2,2}= 'v'; icon{2,3} = 'Instrumental';  
icon{3,1}=skyblue;          icon{3,2}= '*'; icon{3,3} = 'Ice core';
icon{4,1}=ornj;             icon{4,2}= 'o'; icon{4,3} = 'Coral';
icon{5,1}=hotpink;          icon{5,2}= 'd'; icon{5,3} = 'Speleothem';
icon{6,1}=dkgr;             icon{6,2}= 'p'; icon{6,3} = 'Documentary';
icon{7,1}=maroon;           icon{7,2}= 's'; icon{7,3} = 'Sediment';
icon{8,1}=blue;             icon{8,2}= 'h'; icon{8,3} = 'Composite';

network = 'updated2k';
switch network
    case('m08')
        A        = load('itrdbmatrix');
        ptime    = A(4:end,1); nt = length(ptime);
        proxy    = A(:,2:end);
        % Remove Luterbacher reconstructions
        noluter  = proxy(3,:) ~= 2000; % 71 such series
        proxy_nl = proxy(:,noluter); np = size(proxy_nl,2); % define "noluter" proxy matrix
        pinfo    = proxy_nl(1:3,:);  ptype = pinfo(3,:);
        plon     = pinfo(1,:);       plat  = pinfo(2,:);
        clear proxy;
        %
        proxy    = proxy_nl(4:end,:);
        nproxy = zeros(nt,8);
        pcode(round(ptype/1000) == 9) = 1;
        pcode(round(ptype/1000) == 8) = 3;
        pcode(round(ptype/1000) == 7) = 4;
        pcode(round(ptype/100) == 75) = 1;
        pcode(round(ptype/1000) == 6) = 5;
        pcode(round(ptype/1000) == 5) = 6;
        pcode(round(ptype/1000) == 4) = 7;
        pcode(round(ptype/1000) == 3 | round(ptype/1000) == 2) = 8;
       

    case{'pages2k'}
        load ../pages2k_hadcrut4_data_v1.mat
        keep      = setdiff(1:550,proxy.duplic);
        % keep      = proxy.keep;
        prox      = proxy.raw(:,keep);  % no duplicate NA tree rings
        ptype     = proxy.type(keep);
        plon      = proxy.lon(keep);
        plat      = proxy.lat(keep);
        ptime     = proxy.t;
        [nt,np]   = size(prox);
        
        nproxy = zeros(nt,8);
        pcode  = nan(np,1);
        pcode(strcmp('tree ring',ptype)) = 1;
        pcode(strcmp('instrumental',ptype)) = 2;
        pcode(strcmp('ice core',ptype)) = 3;
        pcode(strcmp('coral',ptype)) = 4;
        pcode(strcmp('speleothem',ptype)) = 5;
        pcode(strcmp('documentary',ptype)) = 6;
        pcode(strcmp('lake sediment',ptype)|strcmp('marine sediment',ptype)) = 7;
        pcode(strcmp('historic',ptype)) = 8;

    case{'updated2k'}
        load ../../updated_2k/updated_2k_network.mat
        prox      = proxy.field;  % no duplicate NA tree rings
        ptype     = lower(proxy.ptype);
        plon      = proxy.lon;
        plat      = proxy.lat;
        ptime     = proxy.t;
        [nt,np]   = size(prox);
        
        nproxy = zeros(nt,8);
        pcode  = nan(np,1);
        pcode(strcmp('tree ring',ptype)) = 1;
        pcode(strcmp('instrumental',ptype)) = 2;
        pcode(strcmp('ice core',ptype) | strcmp('ice',ptype)) = 3;
        pcode(strcmp('coral',ptype)) = 4;
        pcode(strcmp('speleothem',ptype)|strcmp('speleo',ptype)) = 5;
        pcode(strcmp('documentary',ptype)) = 6;
        pcode(strcmp('lake sediment',ptype)|strcmp('marine sediment',ptype)|strcmp('sedi',ptype)|strcmp('sediment',ptype)) = 7;
        pcode(strcmp('historic',ptype)) = 8;
        
        
        
end

axlim   = [0 2000 0 1200];
inaxlim = [0 1000 0 90];
ytick   = [25 50 75];


fig('Proxy Availability'),clf
set(gcf,'position',[352    12   954   794])
hmap=axes('Position', [0.05 0.35 0.9 0.65]);
m_proj('Robinson','clong',0);
m_grid('xtick',[-180:60:180],'tickdir','out','ytick',[-90:30:90], 'color',dkgr, 'fontsize',8,'fontname','Times New Roman');

m_coast('color','k');

for k = 1:8
    nproxy(:,k) = sum(~isnan(prox(:,pcode == k)),2);
    if strcmp('pages2k',network)==1
        if sum(pcode == k)~=0
            ind(k)= find(pcode == k, 1, 'last' );
        else
            ind(k) = 0;
        end
    % elseif strcmp('m08',network)==1 && k~=2
    else
        ind(k)= find(pcode == k, 1, 'last' );
    end
        
end

% lon(lon<0) = lon(lon<0) + 360;
for k = 1:np
    hl(k)=m_line(plon(k),plat(k),'color','k','marker',icon{pcode(k),2},'MarkerFaceColor',icon{pcode(k),1},'MarkerSize',8,'linestyle','none');
end
style_t{4}= 30;
title([upper(network) ' network (' num2str(np) ' records)'],style_t{:});


for i = 1:8
    proxycolor(i,:) = icon{i,1};
end

hstack=axes('Position', [0.08 0.1 0.8 0.3]);
bar(ptime,nproxy,'stacked');
set(hstack,'YAxisLocation','Right')
colormap(proxycolor)
axis([0 2000 0 np])
fancyplot_deco('','Time','# proxies',20);

frac=.5;
hstackin=axes('Position', [0.08 0.175 frac*.8 0.18]);
hb = bar(ptime,nproxy,'stacked');
axis(inaxlim)

set(hstackin,'box','off')
set(hstackin,'YAxisLocation','Right','XMinorTick','on')
set(hstackin,'ytick',ytick,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3])


style_l{4} = 14;
% legend for the spatial distribution
[l1,~,OUTH,~]=legend(hl(ind(ind~=0)),icon{pcode(ind(ind~=0)),3}); 
set(l1,style_l{:});
set(l1,'position',[0.8450    0.7506    0.1543    0.2210])
% set(OUTH,'Color','k','MarkerFaceColor','k','Markersize',10);
legend boxoff

% legend for the temporal availability
[l2,OBJH,OUTH,OUTM]=legend(hb(ind~=0),icon{pcode(ind(ind~=0)),3});
set(l2,style_l{:});
set(l2,'position',[0.0791    0.3109    0.1543    0.2210])
legend boxoff

% set(gcf,'color','none')
export_fig([upper(network) '_avail_' num2str(np) '.pdf'],'-cmyk','-r300')


% Instrumental availability only
fig('temporal availability'); set(gcf,'position',[233   394   970   324]); clf;
hb = bar(ptime,nproxy,'stacked');
colormap(proxycolor)
axis([1850 2011 0 np])
fancyplot_deco('Proxy availabity over the istrumental period (1850 - 2011) ','Time (Year A.D.)','# proxies',20);
[l2,OBJH,OUTH,OUTM]=legend(hb(ind~=0),icon{pcode(ind(ind~=0)),3});
set(l2,style_l{:});
set(l2,'position',[0.8658    0.3648    0.1169    0.5485])

export_fig([upper(network) '_avail_' num2str(np) '_instrumental.pdf'],'-cmyk','-r300')
