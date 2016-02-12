%
%   Reconstruct BINNED global mean temperature from proxy composite
%     and investigate robustness to various choices.
  




%---- SPECIFY PARAMETERS HERE!
%set up a sequence of six options

%specify binning parameters
binStep=[25 50 50 100 100 100];
resCutoff=[5 5 5 5 5 5];

%specify normalization parameters
minTime=floor(2000./binStep*.2); %minimum number of timesteps over which to standardize
minPoints=floor(minTime.*.5); %minium number of timesteps within that interval, to be included

%specify screening
%lo res. Options = 'none','basicFilter','hrNeighbors'
screenLoList={'none','basicFilter','none','basicFilter','none','hrNeighbors'};
%hi res. Options = 'none','basicFilter','local','regionalFDR','regional'
screenHiList={'none','regional','none','regional','none','regionalFDR'};
%---- end specification

for o=1:length(binStep)

%find the records that are low vs high res
loResI=find(resMed>=resCutoff(o));
hiResI=find(resMed<resCutoff(o));


%calculate some bin parameters
if rem(2000,binStep)==0
    binEdges=flipud((2000:-binStep(o):0)');
else
    binEdges=flipud((2000:-binStep(o):(rem(2000,binStep(o))-binStep(o)))');
end
binYear=mean([binEdges(2:end) binEdges(1:end-1)],2);
plotEdges=reshape([binEdges(1:end-1) binEdges(2:end)]',[],1);


screenLo=screenLoList{o};
screenHi=screenHiList{o};

%apply screening
%low res
if strcmpi(screenLo,'basicFilter')
    loBini=intersect(idx_qclr,loResI);
elseif strcmpi(screenLo,'hrNeighbors')
    loBini=intersect(idx_qclr_screen,loResI);
else
    loBini=loResI;
end
%hi res
if strcmpi(screenHi,'local')
    hiBini=intersect(intersect(scr_loc,idx_qchr),hiResI);
elseif strcmpi(screenHi,'regionalFDR')
    hiBini=intersect(intersect(scr_fdr,idx_qchr),hiResI);
elseif strcmpi(screenHi,'regional')
    hiBini=intersect(intersect(scr_reg,idx_qchr),hiResI);
elseif strcmpi(screenHi,'basicFilter')
    hiBini=intersect(idx_qchr,hiResI);
else
    hiBini=hiResI;
end

%calculate bins
loBins=nan(length(binYear),length(loBini));

for i=1:length(loBini)
    loBins(:,i)=bin_x(tce',proxy_sgn(:,loBini(i)),binEdges);  
end

hiBins=nan(length(binYear),length(hiBini));
for i=1:length(hiBini)
    hiBins(:,i)=bin_x(tce',proxy_sgn(:,hiBini(i)),binEdges);  
end

%standardize to common mean and variance
%this is a challenge


%--- low res first.

%find interval of common overlap
nBLo=sum(~isnan(loBins),2);
dum=sort(nBLo,'descend');
nStd=dum(minTime(o));

stdR=find(nBLo>=nStd);

%exclude records that don't have any overlap in that period
loS=loBins(stdR,:);
toKeep=find(sum(~isnan(loS))>=minPoints(o));

%
cMean=nanmean(loS(:,toKeep));%calculate mean of each column during interval
cStd=nanstd(loS(:,toKeep));% and std deviation
cMeanMat=repmat(cMean,length(binYear),1);
cStdMat=repmat(cStd,length(binYear),1);

%standardize the matrix
standardLoBins=(loBins(:,toKeep)-cMeanMat)./cStdMat;


%--- now hi res 

%find interval of common overlap
nBHi=sum(~isnan(hiBins),2);
dum=sort(nBHi,'descend');
nStd=dum(minTime(o));

stdR=find(nBHi>=nStd);

%exclude records that don't have any overlap in that period
hiS=hiBins(stdR,:);
toKeep=find(sum(~isnan(hiS))>=minPoints(o));

%
cMean=nanmean(hiS(:,toKeep));%calculate mean of each column during interval
cStd=nanstd(hiS(:,toKeep));% and std deviation
cMeanMat=repmat(cMean,length(binYear),1);
cStdMat=repmat(cStd,length(binYear),1);

%standardize the matrix
standardHiBins=(hiBins(:,toKeep)-cMeanMat)./cStdMat;

%calculate means
loResMean=nanmean(standardLoBins,2);
hiResMean=nanmean(standardHiBins,2);

%scale to instrumental...
%simply by matching the mean and variance of the bins
gBin=bin_x(hadcrut4.t,hadcrut4.gmean,binEdges);
good=find(~isnan(gBin));
gBM=mean(gBin(good));
gBS=std(gBin(good));

loGBM=nanmean(loResMean(good));
loGBS=nanstd(loResMean(good));
hiGBM=nanmean(hiResMean(good));
hiGBS=nanstd(hiResMean(good));

loResMeanScaled=(((loResMean./loGBS)-loGBM).*gBS)+gBM;
hiResMeanScaled=(((hiResMean./hiGBS)-hiGBM).*gBS)+gBM;


%adjust for plotting
plotLo=reshape([loResMeanScaled loResMeanScaled]',[],1);
plotHi=reshape([hiResMeanScaled hiResMeanScaled]',[],1);
plotHAD=reshape([gBin gBin]',[],1);




fig('Global Index Bins'),
if o==1;
    clf;
    set(gcf,'Position',[440   270   852   628])
end
pmax = 600; % scale for # proxies
ylims = [-.5 .5];
% Unscreened HR
ax1 = subplot(3,2,o)
Xlab = 'Year AD'; 
Ylab = 'Temperature (\circC)';
col{1} = rgb('DarkBlue'), col{2} = rgb('DarkGrey'); col{3} = rgb('Red');
plot(plotEdges,plotLo,'color',col{1},'LineWidth',2);
hold on
plot(plotEdges,plotHi,'color',col{2},'LineWidth',2);
hold on
plot(plotEdges,plotHAD,'color',col{3},'LineWidth',2);

%plot N
text(1300,.4,['nHR= ' num2str(size(standardHiBins,2))],'color',col{2},style_l{:})
text(1300,.3,['nLR= ' num2str(size(standardLoBins,2))],'color',col{1},style_l{:})
text(1300,.2,['HadCRUT4'],'color',col{3},style_l{:})



set(gca,'xLim',[tStart tEnd]);
ylabel(Ylab,style_l{:})
if o>=5
xlabel(Xlab,style_l{:})
end

title({['HR screen = ' screenHi '; LR screen = ' screenLo];[num2str(binStep(o)) '-yr bins, \Deltat ' num2str(resCutoff(o)) '-yr cutoff']},style_t{:})

set(gca,'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
end
hepta_figprint(['./figs/global_binned_composite_options'])
eps2pdfMac(['./figs/global_binned_composite_options.eps'])
