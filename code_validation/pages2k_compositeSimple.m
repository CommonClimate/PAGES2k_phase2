% This file generates a simple composite of the PAGES2k data
%  example: http://futureearth.org/blog/2017-jul-11/new-dataset-provides-most-complete-look-yet-climate-last-2000-years

%specify binning parameters
binStep= 30
resCutoff= 5

%specify normalization parameters
minTime=floor(2000./binStep*.2); %minimum number of timesteps over which to standardize
minPoints=floor(minTime.*.5); %minium number of timesteps within that interval, to be included

% Graphical definitions
Xlab   = 'Year AD';   Ylab = 'Composite (\circC)';
col{1} = rgb('SteelBlue'),  col{2} = rgb('Red');
nboot  = 2000; % # of bootstrap samples
ylims  = [-.5 .5]; %pmax = 600; % scale for # proxies

%calculate some bin parameters
if rem(2000,binStep)==0
    binEdges=flipud((2000:-binStep:0)');
else
    binEdges=flipud((2000:-binStep:(rem(2000,binStep)-binStep))');
end
binYear=mean([binEdges(2:end) binEdges(1:end-1)],2);
plotEdges=reshape([binEdges(1:end-1) binEdges(2:end)]',[],1);
%
% bin high and low-res records
comp = bin_x(tce',proxy_sgn,binEdges);

% composite matrices by averaging across columns, ignoring missing values
compMean = nmean(comp,2);  compBoot = bootstrp(nboot,@nmean,comp');
cicomp  = quantile(compBoot',[0.025 0.975],2);

%scale to instrumental...
%simply by matching the mean and variance of the bins
gBin=bin_x(hadcrut4.t,hadcrut4.gmean,binEdges);
good=find(~isnan(gBin)); gBM=mean(gBin(good)); gBS=std(gBin(good));

% find mean and standard deviation
GBM=nanmean(compMean(good)); GBS=nanstd(compMean(good));
% rescale composites
compMeanScaled = (((compMean./GBS)-GBM).*gBS)+gBM;
cicompScaled   = (((cicomp  ./GBS)-GBM).*gBS)+gBM;

%adjust for plotting
plotComp=reshape([compMeanScaled compMeanScaled]',[],1);
plotHAD=reshape([gBin gBin]',[],1);
cicomp_lo  = reshape([cicompScaled(:,1) cicompScaled(:,1)]',[],1);
cicomp_hi  = reshape([cicompScaled(:,2) cicompScaled(:,2)]',[],1);

fig('Global Index Bins'), clf;
%set(gcf,'Position',[440   270   852   628])
% plot bootstrap CIs
wide = (~isnan(cicomp_lo) & cicomp_hi-cicomp_lo>range(ylims)/50.0);  %
area_fill(plotEdges(wide)',cicomp_lo(wide)',cicomp_hi(wide)',col{1},col{1},0.2); hold on
% plot composites
line(plotEdges,plotComp,'color',col{1},'LineWidth',2); % Proxies
%line(plotEdges,plotHAD,'color',col{2},'LineWidth',2); % HadCRUT4

%plot N
%text(1300,.1,['nRecords= ' num2str(size(comp,2))],'color',col{1},style_l{:})
%text(1300,0,['HadCRUT4'],'color',col{2},style_l{:})
%
set(gca,'xLim',[tStart tEnd]); set(gca,'Ylim',[-0.4 0.2]);
ttl = ['PAGES 2k global composite, scaled to temperature (',int2str(binStep),'y averages)'];
fancyplot_deco(ttl,Xlab,Ylab,14,'Helvetica');

froot = ['../figs/compositeGlobalBins_web'];
export_fig([froot '.png'],'-r200','-nocrop','-cmyk','-painters');
%hepta_figprint(froot)
%eps2pdfMac([froot '.eps'])
