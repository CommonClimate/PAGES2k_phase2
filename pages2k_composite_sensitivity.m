

% proxy matrix retained for plots
if strcmp(options.source,'unsmoothed')
    proxy_r = proxy(:,idx_q);
elseif strcmp(options.source,'smoothed')
    proxy_r = proxy_fs(:,idx_q);
else 
    disp('In your dreams! I don''t do that thing')
end

%% 1) composite by archive type
arch = pages2k.archive(idx_q);
% use only those archive types that passed QC
archType = unique(arch);
q_code = p_code(idx_q); 
u_code = unique(q_code);
q_code_freq = hist(q_code,u_code);
u_code_n =  u_code(q_code_freq > 2); % only consider composites with >2 series
select = find(ismember(q_code, u_code_n));

% KLUDGE alert: manually exlude floating coral sequences (REMOVE once it's
% done upstream)
namesq = names(idx_q);
Dam = find(~cellfun(@isempty,strfind(namesq,'Damassa'))); % could Matlab's string parsin possibly be any more bulky? 
Pal = find(~cellfun(@isempty,strfind(namesq,'Palmyra')));
select = setdiff(select,[Dam; Pal],'stable');
%
s_code = q_code(select);
proxy_fss = proxy_r(:,select);
u_code = unique(s_code);
archType = unique(arch(select));
% find 
nType = length(archType);
style_l = style_t; style_l{4} = 12;

fig('By archive type'), clf
set(gcf,'Position',[680   512   702   586])
addpath('/Applications/MATLAB_R2015b.app/toolbox/stats/stats/'); % make sure it uses the right dir
n_tresh = 10; xlims = [0 2000]; 
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
    proxy_arch = proxy_fss(:,s_code == u); 
    
    p_arch(k)  = size(proxy_arch,2);
    n_arch{k}  = sum(~isnan(proxy_arch),2); 
    pcom_arch{k} = nmean(proxy_arch,2); % composite
    %jack_mean = jackknife(@nmean,proxy_arch');
%     pstd_arch{k} = std(jack_mean)'; %  standard dev of jackknifed mean
%     pstd_arch_n{k} = std(jack_mean)'/nstd(pcom_arch{k})*100; % same, normalized
%     pcomJackMin{k} = min(jack_mean);
%     pcomJackMax{k} = max(jack_mean);
    % plotting threshold
    thresh = (n_arch{k}  >= n_tresh);
    
    % plot solid line otherwise, and number of proxies
    [ax,h1,h2]  = plotyy(tce,n_arch{k},tce,pcom_arch{k},@bar,@plot); hold on
    set(h1,'edgecolor',rgb('Gainsboro')); % set(h1,'edgealpha',0.5);
    set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','TickDir','out','YMinorTick','on')
    set(get(ax(1),'Ylabel'),'String','# records')
    set(get(ax(1),'Ylabel'),style_l{:})
    set(h2,'visible','off'); %make invisible for now;
    set(h2,'color',Graph{u,1},'linewidth',1);
    set(ax(2),'Ycolor',Graph{u,1},'YAxisLocation','left')
    set(get(ax(2),'Ylabel'),'String','Comp.');
    set(get(ax(2),'Ylabel'),style_l{:})
    set(gcf,'CurrentAxes',ax(2))
    set(h2,'XData',tce(thresh));  
    set(h2,'YData',pcom_arch{k}(thresh));
    refreshdata(h2); 
    set(h2,'visible','on');
    %plot transparent line for entire period
    hp = patchline(tce,pcom_arch{k},'edgecolor',Graph{u,1},'linewidth',2,'edgealpha',0.4); 
    % plot solid line when enough data are present
    uistack(h2,'top')
    % more cosmetics
    linkaxes([ax(1),ax(2)],'x'); set(ax(1),'xLim',xlims);
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XTick',[])
    ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
    title([archType{k} ', ' int2str(p_arch(k)) ' records'],style_l{:});
%     % plot area
%     axes(ax(2))
%     p_up = pcomJackMin{k};
%     p_dn = pcomJackMax{k};
    %ha = area_fill(tce,p_up,p_dn,Graph{u,1},Graph{u,1},0.3,0.3); 
end

hepta_figprint(['./figs/compositeByArchive_' opstring '_' options.source],800)


save(fout)
