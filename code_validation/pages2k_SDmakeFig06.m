
%     Make Fig06 of the Scientific Data paper
% load arrays called multiple times
signif_n = pages2k.signif_n;
rho_n = pages2k.rho_n;

ulon   = gridcorr.ulon;
ulat   = gridcorr.ulat;
ggrid  = gridcorr.ggrid;

% now plot it
fig('Bete immonde'); clf;
radius = 2000;
nrad   = length(radius);
np     = length(S);
p_code = pages2k.p_code;


% compute summary correations
rMaxI = NaN(np,nrad);  % correlations to instrumental temperature
rMedI = rMaxI; rMeanI = rMaxI;
rMaxP = NaN(np,nrad);  % correlations to calibratable proxies
rMedP = rMaxP; rMeanP = rMaxP;

for r = 1:np
    % index of significant temparature neighbors
    isSignif = pcorr.sig_mat(:,r); isSignif(isnan(isSignif))=0; % turn NaNs to zeros

    % subset by radius
    for ir = 1:nrad
        if pages2k.instCalib(r) % 1) PROXIES THAT OVERLAP WITH THE INSTRUMENTAL RECORD
            isWithin = (pcorr.d_pt(r,:) <= radius(ir))';
            select = logical(isSignif.*isWithin);
            if ~isempty(find(select))
                rMaxI(r,ir)  =    max(abs(pcorr.rho_mat(select,r)));
                rMedI(r,ir)  = median(abs(pcorr.rho_mat(select,r)));
                rMeanI(r,ir) =   mean(abs(pcorr.rho_mat(select,r)));
            end

        else % 2) PROXIES THAT DO NOT
            isWithin = (pages2k.ds{r} <= radius(ir));

            % (there probably is a more elegant way to do this via logical indexing).
            if ~isempty(signif_n{r}) & sum(signif_n{r})>0;
                select = logical(signif_n{r}.*isWithin);
                if sum(select)>0
                    rMaxP(r,ir)  =    max(abs(rho_n{r}(select)));
                    rMedP(r,ir)  = median(abs(rho_n{r}(select)));
                    rMeanP(r,ir) =   mean(abs(rho_n{r}(select)));
                end
            end
        end
    end
end


% PLOT CORRELATION SUMMARIES
nc = 50; % number of color contours
scheme = 'Reds'; cx = [0,1]; bfac = 0.5; %brightening factor
pos = subplot_pos(2,2,[0.05 0.05]);
%
for ir = 1:nrad %
    colI = t2c_brewer(rMedI(:,ir),nc,scheme,cx,bfac);
    nNaN(ir) = sum(isnan(rMedI(:,ir)));
    nSig(ir) = sum(~isnan(rMedI(:,ir)));
    colI(isnan(rMedI(:,ir)),:) = repmat([0 0 0],[nNaN(ir) 1]);
    set(gcf,'PaperPositionMode','auto')

    % PROXY
    ax1 = subplot(2,2,1);
    set(ax1,'position',pos(1,:))
    % plot spatial distribution
    m_proj('Robinson','clong',0);
    m_coast('patch',rgb('WhiteSmoke'));
    m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',10,'fontname','Helvetica');

    % loop over records
    siz = 3*ones(np,1); siz(~isnan(rMedI(:,ir))) = 8;
    for r = 1:np
        h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',colI(r,:),'linewidth',1,'MarkerSize',siz(r),'linestyle','none');
    end
    %hcb1 = colorbar2('vert','|\rho|'); set(hcb1,'YTick',[0:0.2:1])
    title(['PROXY-FIELD: R(proxy,MAT), r \leq ',int2str(radius(ir)),' km'],'FontWeight','bold','FontSize',14,'FontName', 'Helvetica')
    y2 = ylim; x2 = xlim; y2 = y2(1,2); x2 = x2(1,1);
    Ht = text(x2,y2,'A');
    set(Ht,'FontWeight','bold','Units','normalized','Position',[-.05 .95],'fontsize',14,'FontName', 'Helvetica') % default, upper left

    %  NON ANNUAL PROXY
    ax3 = subplot(2,2,3);
    set(ax3,'position',pos(3,:))

    % color map
    colP = t2c_brewer(rMedP(:,ir),nc,scheme,cx,bfac);
    nNaNP(ir) = sum(isnan(rMedP(:,ir)));
    nSigP(ir) = sum(~isnan(rMedP(:,ir)));
    colP(isnan(rMedP(:,ir)),:) = repmat([0 0 0],[nNaNP(ir) 1]);

    % plot spatial distribution
    m_proj('Robinson','clong',0);
    m_coast('patch',rgb('WhiteSmoke'));
    m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',10,'fontname','Helvetica');

    % loop over records
    siz = 3*ones(np,1); siz(~isnan(rMedP(:,ir))) = 8;
    for r = 1:np
        h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',colP(r,:),'linewidth',1,'MarkerSize',siz(r),'linestyle','none');
    end
    %hcb3 = colorbar2('vert','|\rho|'), set(hcb3,'YTick',[0:0.2:1]);


    y2 = ylim; x2 = xlim; y2 = y2(1,2); x2 = x2(1,1);
    Ht = text(x2,y2,'C');
    set(Ht,'FontWeight','bold','Units','normalized','Position',[-.05 .95],'fontsize',14,'FontName', 'Helvetica') % default, upper left

    title(['PROXY-PROXY: R(LR proxy, HR neighbors), r \leq ',int2str(radius(ir)),' km'],'FontWeight','bold','FontSize',14,'FontName', 'Helvetica')

    % PROXY LOCAL
    ax2 = subplot(2,2,2);
    set(ax2,'position',pos(2,:))
    gridcorr.localCorrelation(gridcorr.localSignificance<1 | abs(gridcorr.localCorrelation)==1) = NaN;
    rMedL = abs(gridcorr.localCorrelation');
    colL= t2c_brewer(rMedL(:,ir),nc,scheme,cx,bfac);
    nNaNL(ir) = sum(isnan(rMedL(:,ir)));
    nSigL(ir) = sum(~isnan(rMedL(:,ir)));
    colP(isnan(rMedL(:,ir)),:) = repmat([0 0 0],[nNaNL(ir) 1]);

    m_proj('Robinson','clong',0);
    m_coast('patch',rgb('WhiteSmoke')); hold on
    siz = 3*ones(np,1);
    siz(~isnan(rMedL(:,ir))) = 8;
    for r = 1:np
        h(r) = m_line(p_lon(r),p_lat(r),'marker',Graph{p_code(r),2},'MarkerEdgeColor','k','MarkerFaceColor',colL(r,:),'linewidth',1,'MarkerSize',siz(r),'linestyle','none');
    end

    caxis([0 1])
    m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',10,'fontname','Helvetica') %'xtick',[-125 -122 -119 -116]);
    y2 = ylim; x2 = xlim; y2 = y2(1,2); x2 = x2(1,1);
    Ht = text(x2,y2,'B');
    set(Ht,'FontWeight','bold','Units','normalized','Position',[-.05 .95],'fontsize',14,'FontName', 'Helvetica') % default, upper left
    title(['LOCAL: R(proxy,local MAT)'],'FontWeight','bold','FontSize',14,'FontName', 'Helvetica')

    % GRID-CENTRIC
    ax4 = subplot(2,2,4);
    set(ax4,'position',pos(4,:))

    m_proj('Robinson','clong',0);
    ulonl = [ulon;ulon(end)+5];
    gggrid = squeeze(ggrid(2,:,:))';
    gggrid = [gggrid gggrid(:,1)]
    %
    mpxc = m_pcolor2(ulonl,ulat,gggrid); hold on
    m_coast('color','k'); hold on;
    m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',10,'fontname','Helvetica','backcolor',[0.85 0.85 0.85]) %'xtick',[-125 -122 -119 -116]);
    caxis([0 1])
    %
    title(['GRID: R(MAT,proxies), r \leq ',int2str(radius(ir)),' km'],'FontWeight','bold','FontSize',14,'FontName', 'Helvetica')
    set(ax4,'position',pos(4,:))

    y2 = ylim; x2 = xlim; y2 = y2(1,2); x2 = x2(1,1);
    Ht = text(x2,y2,'D');
    set(Ht,'FontWeight','bold','Units','normalized','Position',[-.05 .95],'fontsize',14,'FontName', 'Helvetica') % default, upper left
    %
    set(mpxc,'linestyle','none')
    set(gcf,'units','normalized','position',[0.39 0.53 0.52 0.40])
    %
    figure_name = ['./figs/PAGES2K_v' vers '_correlation_map' int2str(radius(ir))];
    %hepta_figprint(figure_name,400)
    pause
    export_fig([figure_name '.pdf'],'-r400','-nocrop','-cmyk','-painters');
end

% compute fraction of global area that is
[X,Y] = meshgrid(ulon,ulat);
area_weight = cosd(Y);
A = ~isnan(gggrid(:,2:end)).*area_weight;
area_fraction = sum(A(:))/sum(area_weight(:))  % 77% for v 1.9.0; 72% for v1.13.1

sprintf('%3.4g %s %i %s',100*area_fraction,'% of the planet is within ', radius,'km of a proxy')

% fig('colorbar'), clf
% m_proj('Robinson','clong',0);
% ulonl =[ulon;ulon(end)+5];
% gggrid = squeeze(ggrid(2,:,:))';
% gggrid = [gggrid gggrid(:,1)]
% %
% mpxc = m_pcolor2(ulonl,ulat,gggrid); hold on
% m_coast('color','k'); hold on;
% m_grid('xtick',6,'ytick',9,'xticklabel',[ ],'xlabeldir','middle', 'fontsize',10,'fontname','Helvetica','backcolor',[0.85 0.85 0.85]) %'xtick',[-125 -122 -119 -116]);
% caxis([0 1])
% %
% hcb4 = colorbar2('horiz','Absolute Correlation'), set(hcb4,'XTick',[0:0.2:1],'FontWeight','bold','FontSize',12,'FontName', 'Helvetica');
% hepta_figprint('./figs/correlation_map_colorbar.eps',400);
