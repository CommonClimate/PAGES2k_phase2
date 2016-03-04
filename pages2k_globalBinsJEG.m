%
%   Reconstruct BINNED global mean temperature from proxy composite
%     and investigate robustness to various choices.

% Standardize proxies
proxy_s = standardize(proxy_sgn);
%---- SPECIFY PARAMETERS HERE!
% set up a sequence of six options

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
    
    % bin the proxy matrices
    loResBin = bin_x(proxy_s(:,loBini),tce',binEdges);
    hiResBin = bin_x(proxy_s(:,hiBini),tce',binEdges);
     
    %calculate means
    loResMean=nmean(loResBin,2);
    hiResMean=nmean(hiResBin,2);
    
    %scale to instrumental...
    
    % GLOBAL SCALING: simply by matching the mean and variance of the bins
    gBin=bin_x(hadcrut4.gmean,hadcrut4.t,binEdges);
    good=find(~isnan(gBin));
    gBM=mean(gBin(good));
    gBS=std(gBin(good));
    
    loGBM=nanmean(loResMean(good));
    loGBS=nanstd(loResMean(good));
    hiGBM=nanmean(hiResMean(good));
    hiGBS=nanstd(hiResMean(good));
    
    loResMeanScaled=(((loResMean./loGBS)-loGBM).*gBS)+gBM;
    hiResMeanScaled=(((hiResMean./hiGBS)-hiGBM).*gBS)+gBM;
    
    % also implement local scaling
        % TBD
   
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
    
    %plot legend
    text(1300,.4,['nHR= ' num2str(size(standardHiBins,2))],'color',col{2},style_l{:})
    text(1300,.3,['nLR= ' num2str(size(standardLoBins,2))],'color',col{1},style_l{:})
    text(1300,.2,['HadCRUT4'],'color',col{3},style_l{:})
    % limits 
    set(gca,'xLim',[tStart tEnd]);
    ylabel(Ylab,style_l{:})
    if o>=5
        xlabel(Xlab,style_l{:})
    end
    
    title({['HR screen = ' screenHi '; LR screen = ' screenLo];[num2str(binStep(o)) '-yr bins, \Deltat ' num2str(resCutoff(o)) '-yr cutoff']},style_t{:})
    
    set(gca,'Ylim',ylims,'Ytick',[ylims(1):ylims(2)]);
end
hepta_figprint(['./figs/binned_global_composite'])
eps2pdfMac(['./figs/binned_global_composite'])
