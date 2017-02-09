%% Compare the difference of instrumental records on: annual vs seasonal scales.
addpath(genpath('../utilities'))
load JEG_graphics.mat
tmp  = load('../data/had4med_graphem_sp70');
temp = tmp.Xf(1:1968,:); % 164 full calendar years from Jan 1850 - Dec 2013

%% Prepare data
% Annual
annual = (temp(1:12:end,:) + temp(2:12:end,:) + temp(3:12:end,:) + ...
    temp(4:12:end,:) + temp(5:12:end,:) + temp(6:12:end,:) + ...
    temp(7:12:end,:) + temp(8:12:end,:) + temp(9:12:end,:) + ...
    temp(10:12:end,:) + temp(11:12:end,:) + temp(12:12:end,:))./12;
% JJA
summer = (temp(6:12:end,:) + temp(7:12:end,:) + temp(8:12:end,:))./3;
% DJF
winter = (temp(12:12:end,:) + temp(1:12:end,:) + temp(2:12:end,:))./3;
% or possibly seasons by hemisphere


%% compute correlations and make a map
ns      = size(temp,2);
rho_jja = diag(corr(annual,summer)); 
rho_djf = diag(corr(annual,winter));


% Assign coordinate information
x  = -180:5:180; nx = length(x);
y  = -90:5:90;   ny = length(y);

locs = tmp.loc;
inds = (locs+2.5)./5;
ilat = inds(:,2)+18;
ilon = inds(:,1)+36;

% Mapping
r_jja = nan(ny,nx);
r_djf = nan(ny,nx);
for i = 1:ns
    r_jja(ilat(i),ilon(i)) = rho_jja(i);
    r_djf(ilat(i),ilon(i)) = rho_djf(i);
end

% Plot 
pos = subplot_pos(2,2,[0.1 0.1]);
fig('Annual vs. Seasonal');clf;
h = findobj(gca,'Type','patch');  set(h,'FaceColor',rgb('Silver'));
% set(gcf,'position',[189   128   975   496])
set(gcf,'position',[277   100   965   611])
% ===JJA
ha = subplot(221);
pos(1,:) = [0.015,0.51,0.485,0.485];
set(ha,'position',pos(1,:))
m_proj('robinson');
h = m_pcolor(x,y,r_jja); hold on; caxis([0,1]);
set(h,'edgecolor','none'),
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
colormap(pmkmp)
title('HadCRUT4.2 \rho(MAT,JJA)',style_t{:})

h = subplot(222); 
set(h,'position',pos(2,:))
hist(rho_jja,20); xlim([0.2 1])
hs = findobj(gca,'Type','patch');
set(hs,'FaceColor',rgb('Silver'))
jja_m = mean(rho_jja);
hl    = line([jja_m jja_m],ylim);
set(hl,'linestyle','-.','linewidth',2,'color','k')
fancyplot_deco('Distribution (MAT vs. JJA)','','',14,'Times') % mean:',sprintf('%6.2f',jja_m),')'],'','')

% ===DJF
ha = subplot(223); 
%pos(3,:) = [0.015   0.025   0.485  0.485];
pos(3,:) = [0.015   0.0   0.485  0.485];

set(ha,'position',pos(3,:))
m_proj('robinson');
h = m_pcolor(x,y,r_djf); hold on; caxis([0,1]);
set(h,'edgecolor','none'),
m_grid('box','off','xtick',9,'ytick',6,'xlabeldir','end','xticklabels',[],'yticklabels',[]);
m_coast('color','k');
title('HadCRUT4.2 \rho(MAT,DJF)',style_t{:})

hc = colorbar2('horiz','\rho');
set(hc,'position',[0.0150    0.52    0.4850    0.0243])
set(ha,'position',pos(3,:))

h = subplot(224); 
set(h,'position',pos(4,:))
hist(rho_djf,20); xlim([0.2 1])
hs = findobj(gca,'Type','patch'); 
set(hs,'FaceColor',rgb('Silver'))
djf_m = mean(rho_djf);
hl    = line([djf_m djf_m],ylim);
set(hl,'linestyle','-.','linewidth',2,'color','k')
fancyplot_deco('Distribution (MAT vs. DJF)','','',14,'Times')%, mean:',sprintf('%6.2f',djf_m),')'],'','')

export_fig('../figs/hadcrut4_correlation.pdf','-r300','-cmyk','-painters')
