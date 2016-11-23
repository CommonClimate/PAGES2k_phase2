library(dplyr)
library(R.matlab)
library(fda)
library(ggplot2)
library(reshape2)
require(ggplot2)
require(reshape2)
require("matrixStats")

#change the path
basePath <- "/Users/wan/Documents/PAGES2k_phase2/"
# or use
basePath <- getwd()

data <- readMat(paste( basePath,"data/pages2k_hadcrut4_noDetrend_normal_1_12_2.mat",sep="/"))
fieldnames <- unlist(readMat(paste( basePath,"data/field_names.mat",sep="/")) )


norm_p <- 1 #should proxies be mapped to a standard normal ?    [boolean]
detrend <- 0 #do you want to detrend coral d18O proxies? [boolean]
navlMin <- 20 #what is your threshold for # samples over the Common Era?
tStart <- 1 #define start year (remember: the Common Era does not have a year 0).
tEnd   <- 2000 #define end year for the analysis

# define analysis options
lat_weight <- 0; # are we normalizing by the cosine of latitude? [boolean]
sifting_style <-'qcOnly'; # possible choices: noSift, qcOnly, qcScreenHR, qcScreenLR, qcScreenAll
screenHR_style <-"loc" # 'loc' = local; 'reg' = regional (within 2000km radius), 'fdr' = regional accounting for false discovery rate

#define function %!in% which is same as ~ismember
'%!in%' <- function(x,y)!('%in%'(y,x))

#define rfind function as similar to find in matlab
rfind <- function(x)seq(along=x)[as.logical(x)] 


t<-unlist(data$pages2k[rfind(fieldnames%in%'year')])
tce  <- tStart:tEnd
nce <- length(tce)
T   <- data$pages2k[rfind(fieldnames%in%'S')][[1]]
nr <- length(T[1,1,])
names <- T[18,,]
yearMin <- unlist(data$pages2k[rfind(fieldnames%in%'yearMin')])
yearMax <- unlist(data$pages2k[rfind(fieldnames%in%'yearMax')])
resMed  <- unlist(data$pages2k[ rfind(fieldnames%in%'resMed')])
resAvg  <- unlist(data$pages2k[ rfind(fieldnames%in%'resAvg')])
resMax  <- unlist(data$pages2k[ rfind(fieldnames%in%'resMax')])
p_code  <- unlist(data$pages2k[rfind(fieldnames%in%'p_code')])
Graph   <- data$pages2k[32][[1]][12:22]
p_lon   <- unlist(data$pages2k[rfind(fieldnames%in%'loc')][[1]][,1])
p_lat   <- unlist(data$pages2k[rfind(fieldnames%in%'loc')][[1]][,2])
edgec   <- unlist(data$pages2k[rfind(fieldnames%in%'edgec')])
archive <- unlist(data$pages2k[rfind(fieldnames%in%'archive')])
S <- data$pages2k[rfind(fieldnames%in%'S')][[1]]
incl <- 1:nr

#define multiplot function cite:http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/)
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# account for various pre-processing choices
if (norm_p & detrend){
  proxy = data$pages2k[rfind(fieldnames%in%'proxy_nda')][[1]][which(t%in%tce),incl]
} else if (norm_p & !detrend) {
  proxy = data$pages2k[rfind(fieldnames%in%'proxy_na')][[1]][which(t%in%tce),incl]
} else if (!norm_p & detrend) {
  proxy = data$pages2k[rfind(fieldnames%in%'proxy_da')][[1]][which(t%in%tce),incl]
}else {proxy = data$pages2k[rfind(fieldnames%in%'proxy_a')][[1]][which(t%in%tce),incl]}

ny=dim(proxy)[1]
nr=dim(proxy)[2]

#exploit temperature interpretation
sgn = T[11,,]

sgn_vec = matrix(1, nr, 1) 
for (i in 1:nr)
  {if (sgn[i]%in%"positive" || sgn[i]%in%'p')  #KLUDGE
        {sgn_vec[i] = +1}
    else if (sgn[i]%in%"negative")
        {sgn_vec[i] = -1}
  }

proxy_sgn = matrix(1, ny, nr) 


#This could take a few minutes
for (i in 1:ny)
  {  
  proxy_sgn[i,] = scale(proxy)[i,]*sgn_vec
}

navl = matrix(1, 1, nr) 
for (i in 1:nr)
{  
  navl[i] = sum(is.finite(proxy_sgn[,i]))
}

ma <- S[4,,]
rfind <- function(x)seq(along=x)[as.logical(x)] 
idx_qchr = rfind("borehole"%!in%ma & resMed<=5 & resAvg<= 5 & yearMin<= 1850 & navl >= navlMin)
idx_qclr = rfind("borehole"%!in%ma & resMed>5 & yearMin<= 1850 & navl>= navlMin)
idx_qc  = union(idx_qclr,idx_qchr)


# screening for significant correlations.
# 1) Calibratable Proxies
scr_reg  = unlist(data$pages2k[rfind(fieldnames%in%'screen_reg')][[1]][1]) #[15] MAT regional correlation screening
scr_fdr  = unlist(data$pages2k[rfind(fieldnames%in%'screen_fdr')][[1]][1]) #[16] MAT regional correlation screening controlling for false discovery rate
scr_loc  = unlist(data$pages2k[rfind(fieldnames%in%'screen_loc')][[1]][1]) #[14] MAT local correlation screening

idx_scr = get(paste("scr_",screenHR_style,sep=""))



# myplot return a plot of column proxies against time. The input argument shall be number of the column that is desired to plot.
myplot<- function(x){
  df <- data.frame(time = tStart:tEnd,proxy[,idx_qchr[x]])
  df <- melt(df, id.vars = "time", variable.name='series')
  result=ggplot(df, aes(time,value)) + 
    geom_line() + 
    stat_smooth() + 
    ggtitle( paste(x," archive=",archive[x],"   dataSetName=",names[x])) +
    theme(plot.title = element_text(size = 8))
  return(result)
}


## An example of using multiplot and myplot to generate plots of certain kinds of proxies
multiplot(myplot(1), myplot(2), myplot(3), myplot(4),myplot(5),myplot(6), cols=2)

# merge indices of screened proxies. Get the wanted indices stored so as to facilitate the filterplot function
HR <- idx_qchr
HR_local <- intersect(scr_loc,idx_qchr)
HR_fdr   <- intersect(scr_fdr,idx_qchr)
HR_regional <- intersect(scr_reg,idx_qchr)
LR <- idx_qclr
LR_local <- intersect(scr_loc,idx_qclr)
LR_fdr   <- intersect(scr_fdr,idx_qclr)
LR_regional <- intersect(scr_reg,idx_qclr)
QC=idx_qc
QC_local = intersect(scr_loc,idx_qc)
QC_fdr = intersect(scr_fdr,idx_qc)
QC_regional = intersect(scr_reg,idx_qc)

##Input the filter that one is willing to apply, i.e."HR_local", "QC_fdr" or "HR_regional" . myplot2 returns a plot of filtered mean proxies against time.
filterplot <- function(x){
  
  df <- data.frame(time = tStart:tEnd,rowMeans(proxy[,get(paste("",x,sep=""))],na.rm=TRUE))
  df <- melt(df, id.vars = "time", variable.name='series')
  result=ggplot(df, aes(time,value)) + 
    geom_line() + 
    stat_smooth() + 
    ggtitle(paste(x,"  n=",length(get(x)),sep="")) +
    labs(x="year",y="proxy",size=1)+
    theme(plot.title = element_text(size=4),
          axis.text.x = element_text(size=4))
  return(result)
}



#Generating and then saving the plot to the destinated location

#Contrasting HR only and HR+LR
png(paste( basePath,"figures/HRandQConly_Contrast.png",sep="/"), width = 4, height = 4, units = 'in', res = 1600)
multiplot(filterplot("HR"), filterplot("HR_regional"),filterplot("HR_fdr"), filterplot("HR_local"), filterplot("QC"), filterplot("QC_regional"),filterplot("QC_fdr"), filterplot("QC_local"),cols=2)
dev.off()


#HR contrast
png(paste( basePath,"figures/HR_Contrast.png",sep="/"), width = 4, height = 3, units = 'in', res = 1600)
multiplot(filterplot("HR_regional"),filterplot("HR_fdr"), filterplot("HR_local"), filterplot("QC_regional"), filterplot("QC_fdr"),filterplot("QC_local"),cols=2)
dev.off()


#LR contrast
png(paste( basePath,"figures/LR_Contrast.png",sep="/"), width = 4, height = 3, units = 'in', res = 1600)
multiplot(filterplot("LR_regional"),filterplot("LR_fdr"), filterplot("LR_local"), filterplot("QC_regional"), filterplot("QC_fdr"),filterplot("QC_local"),cols=2)
dev.off()

