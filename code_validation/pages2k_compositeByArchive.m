
% Purpose: investigate how the composite changes when binning sites by
%  archive type
% Contributed by: Julien Emile-Geay


% bin parameters
binStep = 50;
binEdges=flipud((2000:-binStep:0)');
binYear=mean([binEdges(2:end) binEdges(1:end-1)],2);
plotEdges=reshape([binEdges(1:end-1) binEdges(2:end)]',[],1);
tb = binEdges(2:end);


%% 1) composite by archive type
archive = {pages2k.S.archiveType};
archive = strrep(archive,'Tree Ring','tree');

arch = archive(idx_q);
proxy_r = proxy_sgn(:,idx_q);
% use only those archive types that passed QC
archType = unique(arch);
q_code = p_code(idx_q); 
u_code = unique(q_code);
q_code_freq = hist(q_code,u_code);
u_code_n =  u_code(q_code_freq >= 5); % only consider composites with >=5 series
select = find(ismember(q_code, u_code_n));

% % KLUDGE alert: manually exlude floating coral sequences (REMOVE once it's
% % done upstream)
% namesq = names(idx_q);
% Dam = find(~cellfun(@isempty,strfind(namesq,'Damassa'))); % could Matlab's string parsin possibly be any more bulky? 
% Pal = find(~cellfun(@isempty,strfind(namesq,'Palmyra')));
% select = setdiff(select,[Dam; Pal],'stable');
%
s_code = q_code(select);
proxy_s = proxy_r(:,select);
u_code = unique(s_code);
archType = unique(arch(select));
nType = length(archType); nType = length(u_code);


% analysis parameters
n_tresh = 10; nboot = 500; % # of bootstrap samples.

% plotting parameters
style_l = {'FontName','Helvetica','Fontweight','bold','Fontsize',12}
xlims = [0 2000];  c = 1.1; % expansion factor for tick marks
pnums = [50 20 60 50 60 500]; dp = [10 5  10 10 10 100]; % custom cosmetic ticks

fig('By archive type'), clf
set(gcf,'Position',[680   512   702   586])
%addpath('/Applications/MATLAB_R2013b.app/toolbox/stats/stats/'); % make sure it uses the right dir

for k = 1:nType
    if ismember(nType,[9 10]);
        ax_arch(k) = subplot(5,2,k);
    elseif ismember(nType,[7 8]);
        ax_arch(k) = subplot(4,2,k);
    elseif ismember(nType,[5 6]);
        ax_arch(k) = subplot(3,2,k);
    elseif ismember(nType,[3 4]);
        ax_arch(k) = subplot(2,2,k);
    end
    u = u_code(k); 
    
    proxy_arch = bin_x(tce',proxy_s(:,s_code == u),binEdges);
    
    p_arch(k)     = size(proxy_arch,2);    
    n_arch{k}     = sum(~isnan(proxy_arch),2); 
    pcom_arch{k}  = nmean(proxy_arch,2); 
    p_boot        = bootstrp(nboot,@nmean,proxy_arch');
    ci_arch{k}    = quantile(p_boot',[0.025 0.975],2);
    %ci_arch{k}  = bootci(nboot,@nmean,proxy_arch'); % cannot handle NaNs. Boooo!
    
    n_archPlot = reshape([n_arch{k} n_arch{k}]',[],1);
    pcom_archPlot = reshape([pcom_arch{k} pcom_arch{k}]',[],1);
    ciLo        = reshape([ci_arch{k}(:,1) ci_arch{k}(:,1)]',[],1);  
    ciHi        = reshape([ci_arch{k}(:,2) ci_arch{k}(:,2)]',[],1);
    ylims       = [floor(c*min(ciLo)), ceil(c*max(ciHi))];
    % set tick mark spacing
    if range(ylims) <= 2, 
        dy = 0.5;
    else
        dy = 1;
    end
    yticks      = ylims(1):dy:ylims(2);  % define tick marks
    % plotting threshold
    thresh = (n_archPlot  >= n_tresh);
    % plot solid line otherwise, and number of proxies
    [ax,h1,h2]  = plotyy(tb,n_arch{k},plotEdges,pcom_archPlot,@bar,@line); hold on
    set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
    set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right');%'TickDir','out','YMinorTick','on')
    set(get(ax(1),'Ylabel'),'String','# records')
    set(get(ax(1),'Ylabel'),style_l{:}), set(ax(2),'yLim',ylims)
    pticks = [0:dp(k):pnums(k)];
    set(ax(1),'yLim',[0 pnums(k)],'YTick',pticks);
    %set(h2,'visible','off'); %make invisible for now;
    set(h2,'color',Graph{u,1},'linewidth',2,'linestyle',':');
    set(ax(2),'Ycolor',Graph{u,1},'YAxisLocation','left','Ytick',yticks)
    set(get(ax(2),'Ylabel'),'String','Composite');
    set(get(ax(2),'Ylabel'),style_l{:})
    set(gcf,'CurrentAxes',ax(2))
    set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02])

    % plot bootstrap CI
    wide = (~isnan(ciLo) & ciHi-ciLo>range(ylims)/50.0);  %
    hci = area_fill(plotEdges(wide)',ciLo(wide)',ciHi(wide)',Graph{u,1},Graph{u,1},0.2);
    %plot transparent line for entire period
    hp = line(plotEdges(thresh),pcom_archPlot(thresh),'color',Graph{u,1},'linewidth',2,'linestyle','-','Parent',ax(2));
    % plot solid line when enough data are present
    uistack(h2,'top')
    % more cosmetics
    linkaxes([ax(1),ax(2)],'x'); set(ax(1),'xLim',xlims); 
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XTick',[])
    ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
    title([archType{k} ', ' int2str(p_arch(k)) ' records'],style_l{:});
end

froot = ['../figs/' opstring '_compositeByArchive'];
%export_fig([froot '.pdf'],'-r300','-nocrop','-cmyk','-painters'); % no transparency
%print(gcf, '-dpdf', '-cmyk','-r300',[froot '.pdf']);

hepta_figprint(froot,300)  %export to PDF as one needs to adjust transparency manually in Illustrator.
%eps2pdfMac([froot '.eps'])

save(f_out)



% % troubleshoot weird U-shape of coral composite
% % culprit: longest record
% coral = proxy_s(:,s_code == 2); 
% 
% % find the first row
% ncoral = nsum(~isnan(coral),2);
% earliest = find(ncoral == 1,1,'first')
% find(~isnan(coral(earliest,:))) 
% 
% %trace the litle scoundrel
% coral_ind = find(s_code == 2);
% scoundrel = pages2k.S(idx_q(select(coral_ind(1))));
% % it is : 'Gingerbreads Bahamas coral growth data', paleoData_TSid: 'Ocean2kHR_130'
% 
% fig('Coral composite'), clf 
% plot(tce,pcom_arch{1},'linewidth',2,'color',Graph{2,1}), hold on
% plot(tce,coral(:,1),'color',rgb('DimGray'))
% lab{1} = ['PAGE2k composite, v' versl]; 
% lab{2} = 'Gingerbreads Bahamas coral growth data';
% xlim([1550, 2000]); hl = legend(lab{:},style_l{:}), set(hl, 'box','off')
% fancyplot_deco('Coral composite U-ness investigation','Time','Z-score')
% hepta_figprint(['../figs/coral_composite' opstring '_' smoothString],400)




