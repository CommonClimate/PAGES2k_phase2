binStep=[5 10 25 50];

for o=1:4
%calculate some bin parameters
if rem(2000,binStep)==0
    binEdges=flipud((2000:-binStep(o):0)');
else
    binEdges=flipud((2000:-binStep(o):(rem(2000,binStep(o))-binStep(o)))');
end
binYear=mean([binEdges(2:end) binEdges(1:end-1)],2);
plotEdges=reshape([binEdges(1:end-1) binEdges(2:end)]',[],1);

gBin=bin_x(hadcrut4.gmean,hadcrut4.t,binEdges);

plotHAD=reshape([gBin gBin]',[],1);

%plot 
subplot(2,2,o)
cla
plot(hadcrut4.t,hadcrut4.gmean,'k')
hold on
plot(plotEdges,plotHAD,'r')
title(['Raw vs binned HadCRUT4; bin = ' num2str(binStep(o))])
ylabel('\circC')
xlabel('Year (AD)')
xlim([1850 2000])
ylim([-.4 .6])
end

hepta_figprint('./figs/HadCRUT_test')
eps2pdfMac('./figs/HadCRUT_test.eps')