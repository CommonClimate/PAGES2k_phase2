clear all; close all;  % clean slate
vers = '1_9_0';  % version of the database to be use
addpath(genpath('./utilities'))% load code utilities 

% define analysis options
lat_weight = 0; % are we normalizing by the cosine of latitude? [boolean]
sifting_style = 'qcScreenHR'; % possible choices: noSift, qcOnly, qcScreenHR, qcScreenLR, qcScreenAll
screenHR_style = 'reg';  % 'loc' = local; 'reg' = regional (within 2000km radius), 'fdr' = regional accounting for false discovery rate

norm_p = 1;  % should proxies be mapped to a standard normal ?    [boolean] 
detrend = 0; % do you want to detrend coral d18O proxies? [boolean]  
navlMin = 20; % what is your threshold for # samples over the Common Era?
tStart = 1; % define start year (remember: the Common Era does not have a year 0). 
tEnd   = 2000; %  define end year for the analysis

% define I/O files
f_out = ['./data/pages2k_composite_' vers '.mat'];
% replace with Cowtan & Way median
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

f_merged = ['./data/pages2k_hadcrut4_' d_str '_' g_str '_' vers];

% define graphical defaults
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultAxesFontWeight','normal')  
set(0,'defaultaxesfontsize',12); set(0,'defaulttextfontsize',12);
% Now let the wild rumpus begin

%% STAGE 1: load relevant data files and prepare the data
pages2k_composite_prep
SDmakeTable 

%% STAGE 2: create global binned composites, and test their sensitivity to methodological choices
pages2k_composite_globalBins

%% STAGE 3: sensitivity analysis of the composite

% stratification by archive type
pages2k_compositeByArchive
% stratification by latitude band
pages2k_compositeByLatitudeBand
% sensitivity to record length
pages2k_composite_recordLength

% stratification however you want !
  % write your own and share it with the PAGES2k community
  




