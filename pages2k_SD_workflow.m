clear all; close all;  % clean slate
vers = '1_12_0';  % version of the database to be use
addpath(genpath('./utilities'))% load code utilities 

% define analysis options
lat_weight = 0; % are we normalizing by the cosine of latitude? [boolean]
sifting_style = 'qcOnly'; % possible choices: noSift, qcOnly, qcScreenHR, qcScreenLR, qcScreenAll
screenHR_style = 'reg';  % 'loc' = local; 'reg' = regional (within 2000km radius), 'fdr' = regional accounting for false discovery rate

norm_p = 1;  % should proxies be mapped to a standard normal ?    [boolean] 
detrend = 0; % do you want to detrend coral d18O proxies? [boolean]  
navlMin = 20; % what is your threshold for # samples over the Common Era?
tStart = 1; % define start year (remember: the Common Era does not have a year 0). 
tEnd   = 2000; %  define end year for the analysis

% define I/O files
f_out = ['./data/pages2k_composite_' vers '.mat'];

if detrend
    d_str = 'detrend';
else
    d_str = 'noDetrend';
end

if norm_p 
    g_str = 'normal';
else
    g_str = 'raw';
end

% LOAD  merged proxy /temperature data + output of correlation analyses
load(['./data/pages2k_hadcrut4_' d_str '_' g_str '_' vers])

% define graphical defaults
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultAxesFontWeight','normal')  
set(0,'defaultaxesfontsize',12); set(0,'defaulttextfontsize',12);

%% prep work
pages2k_SD_prep                  % makes SD Fig S2

% Now let the wild rumpus begin

%% STAGE 1: visualize dataset and produce Table 1
pages2k_SDmakeFig01   
pages2k_SDmakeFig02
pages2k_SDmakeTable

%% STAGE 2: create global binned composites, and test their sensitivity to methodological choices
pages2k_composite_globalBins      % makes SD Fig 03

%% STAGE 3: sensitivity analysis of the composite

% stratification by archive type
pages2k_compositeByArchive         % makes SD Fig 04
% stratification by latitude band
pages2k_compositeByLatitudeBand    % makes SD Fig S8 
% sensitivity to Screening Criterion
pages2k_compositeByScreenCrit      % makes SD Fig S4 
% sensitivity to record length
pages2k_composite_recordLength     % makes SD Fig S3 



% stratification however you want !
  % write your own and share it with the PAGES2k community

  
  
  
  
  
  
  
  
  
  
% extra diagnostics  
resAvg_nontrees = resAvg(p_code <10);
mean(resAvg_nontrees)
resAvg_seds = resAvg(p_code < 8 & p_code > 5);
mean(resAvg_seds)

