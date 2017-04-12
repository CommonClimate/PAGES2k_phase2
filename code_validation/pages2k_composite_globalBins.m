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

% Graphical definitions
Xlab   = 'Year AD';   Ylab = 'Temperature (\circC)';
col{1} = rgb('DarkBlue'), col{2} = rgb('DimGray'); col{3} = rgb('Red');
nboot  = 500; % # of bootstrap samples
ylims  = [-.5 .5]; %pmax = 600; % scale for # proxies

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
    %
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
    
    % bin high and low-res records
    loBins = bin_x(tce',proxy_sgn(:,loBini),binEdges);
    hiBins = bin_x(tce',proxy_sgn(:,hiBini),binEdges);
    
    % composite matrices by averaging across columns, ignoring missing values
    loResMean = nmean(loBins,2);  loResBoot = bootstrp(nboot,@nmean,loBins');
    ciLoRes  = quantile(loResBoot',[0.025 0.975],2);
    hiResMean = nmean(hiBins,2);  hiResBoot = bootstrp(nboot,@nmean,hiBins');
    ciHiRes  = quantile(hiResBoot',[0.025 0.975],2);
    
    %scale to instrumental...
    %simply by matching the mean and variance of the bins
    gBin=bin_x(hadcrut4.t,hadcrut4.gmean,binEdges);
    good=find(~isnan(gBin)); gBM=mean(gBin(good)); gBS=std(gBin(good));
    
    % find mean and standard deviation
    loGBM=nanmean(loResMean(good)); loGBS=nanstd(loResMean(good));
    hiGBM=nanmean(hiResMean(good)); hiGBS=nanstd(hiResMean(good));
    % rescale composites
    loResMeanScaled = (((loResMean./loGBS)-loGBM).*gBS)+gBM;
    ciLoResScaled   = (((ciLoRes  ./loGBS)-loGBM).*gBS)+gBM;
    hiResMeanScaled = (((hiResMean./hiGBS)-hiGBM).*gBS)+gBM;
    ciHiResScaled   = (((ciHiRes  ./hiGBS)-hiGBM).*gBS)+gBM;
     
    %adjust for plotting
    plotLo=reshape([loResMeanScaled loResMeanScaled]',[],1);
    plotHi=reshape([hiResMeanScaled hiResMeanScaled]',[],1);
    plotHAD=reshape([gBin gBin]',[],1);
    ciLoRes_lo  = reshape([ciLoResScaled(:,1) ciLoResScaled(:,1)]',[],1);  
    ciLoRes_hi  = reshape([ciLoResScaled(:,2) ciLoResScaled(:,2)]',[],1);
    ciHiRes_lo  = reshape([ciHiResScaled(:,1) ciHiResScaled(:,1)]',[],1);  
    ciHiRes_hi  = reshape([ciHiResScaled(:,2) ciHiResScaled(:,2)]',[],1);
 
    fig('Global Index Bins'),
    if o==1;
        clf;
        set(gcf,'Position',[440   270   852   628])
    end
    
    ax1 = subplot(3,2,o)
    % plot bootstrap CIs
    wide = (~isnan(ciLoRes_lo) & ciLoRes_hi-ciLoRes_lo>range(ylims)/50.0);  %
    area_fill(plotEdges(wide)',ciLoRes_lo(wide)',ciLoRes_hi(wide)',col{1},col{1},0.2); hold on
    wide = (~isnan(ciHiRes_lo) & ciHiRes_hi-ciHiRes_lo>range(ylims)/50.0);  %
    area_fill(plotEdges(wide)',ciHiRes_lo(wide)',ciHiRes_hi(wide)',col{2},col{2},0.2); 
    % plot composites
    line(plotEdges,plotLo,'color',col{1},'LineWidth',2); % LR   
    line(plotEdges,plotHi,'color',col{2},'LineWidth',2); % HR
    line(plotEdges,plotHAD,'color',col{3},'LineWidth',2); % HadCRUT4
    
    %plot N
    text(1300,.4,['nHR= ' num2str(size(hiBins,2))],'color',col{2},style_l{:})
    text(1300,.3,['nLR= ' num2str(size(loBins,2))],'color',col{1},style_l{:})
    text(1300,.2,['HadCRUT4'],'color',col{3},style_l{:})   
    %
    set(gca,'xLim',[tStart tEnd]);
    set(gca,'Ylim',ylims); %'Ytick',[ylims(1):ylims(2)]);
    ttl = {['HR screen = ' screenHi '; LR screen = ' screenLo];[num2str(binStep(o)) '-yr bins, \Deltat ' num2str(resCutoff(o)) '-yr cutoff']};
    %ylabel(Ylab,style_l{:})
    if o>=5
        fancyplot_deco(ttl,Xlab,Ylab,14,'Helvetica');
        %xlabel(Xlab,style_l{:})
    else
        fancyplot_deco(ttl,'',Ylab,14,'Helvetica');
    end   
    
end
froot = ['.../figs/' opstring '_compositeGlobalBins'];
export_fig([froot '.pdf'],'-r200','-nocrop','-cmyk','-painters');
%hepta_figprint(froot)
%eps2pdfMac([froot '.eps'])

