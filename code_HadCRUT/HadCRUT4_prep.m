% annualizes the HadCRUT4.2 monthly data and makes Fig 03 of the Data Descriptor.
addpath(genpath('../utilities'));
% Load raw HadCRUT4 temperature data
H = load('../data/had4med')
% contains:
%H4med - raw median; structure contains unadulterated data as extracted from .nc file
%rawH - raw median; reshaped to 2D for infilling
%had4med - raw median w/ satellite data from Cowtan & Way 2014; 2D
% ta: time axis in datenum format
t_lon = H.loc(:,1); t_lat = H.loc(:,2);
HadCRUT4_med = H.had4med; [nt,pt] = size(HadCRUT4_med);
HadCRUT4_avail = reshape(~isnan(HadCRUT4_med),size(HadCRUT4_med));
completeness = sum(HadCRUT4_avail,1)/nt;
time = H.ta;
% define completeness threshold
thre = 0.3;

% Load GraphEM-interpolated HadCRUT4 temperature data
Hi = load('../data/had4med_graphem_sp70.mat')
idx = (completeness>=thre);
HadCRUT4_med_i = Hi.Xf(:,idx); pc = size(HadCRUT4_med_i,2);
fig('Grid boxes with enough data');
easy_scatter_map(H.loc(idx,:));

%% annualize
tv = datevec(time); 
ta = unique(tv(:,1));  na = length(ta);
HadCRUT4_med_annual = nan(na,pc);
ms = 1; me = 12;
% loop over grid points
for j = 1:pc
    for k = 1:na;
        ms_h = datenum([ta(k),ms,1,0,0,0]);
        me_h = datenum([ta(k),me,31,0,0,0]);
        dm   = HadCRUT4_med_i(ismember(time,ms_h:me_h),j);
        f = sum(~isnan(dm))/numel(dm);
        if f >= thre
            HadCRUT4_med_annual(k,j) = nmean(dm);
        end
    end
end

%% compute global averages
weight = repmat(cosd(t_lat)',[nt 1]);
weight_a = repmat(cosd(t_lat)',[na 1]);
HadCRUT4_med_avg = nmean(HadCRUT4_med.*weight,2);
HadCRUT4_med_avg_i = nmean(HadCRUT4_med_i.*weight(:,idx),2);
HadCRUT4_med_avg_a = nmean(HadCRUT4_med_annual.*weight_a(:,idx),2);
tan = datenum([ta, 6*ones(na,1),15*ones(na,1)]);

%% plot
fig('GMT'),clf
hr = plot(time,HadCRUT4_med_avg,'color',rgb('black')), hold on 
xlim([datenum(1845,1,1) datenum(2020,1,1)])
hi = plot(time,HadCRUT4_med_avg_i,'color',rgb('Crimson')), 
ha = plot(tan,HadCRUT4_med_avg_a,'ko','linewidth',1), hold off
lab{1} = 'Raw Median'; lab{2} = 'GraphEM, $sp = 0.7\%$ '; lab{3} = 'same, annual';
datetick('keeplimits'), hl = legend([hr hi ha],lab{:}), set(hl,'Interpreter','Latex','box','off','location','southeast','Fontsize',12);
fancyplot_deco('HadCRUT4 global mean','Time','Anomaly (K)');
export_fig('../figs/HadCRUT4_GM_monthly_GraphEMsp70.pdf','-cmyk','-r200','-painters');

%% export
loc = H.loc(idx,:);
field = HadCRUT4_med_annual;
gmean = HadCRUT4_med_avg_a;
t = ta;
save('../data/had4med_graphem_sp70_annual','field','t','gmean','loc')


 
