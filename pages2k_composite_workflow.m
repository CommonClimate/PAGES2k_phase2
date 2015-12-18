clear all; close all;  % clean slate
vers = '1_7_1';  % version of the database to be use
addpath(genpath('./utilities'))% load code utilities 

% define analysis options
lat_weight = 0; % are we normalizing by the cosine of latitude? [boolean]
sifting_style = 'qcScreenHR'; % possible choices: noSift, qcOnly, qcScreenHR, qcScreenLR, qcScreenAll
norm_p = 0;  % do you want proxies Gaussianized ?    [boolean] 
detrend = 0; % do you want to detrend coral d18O proxies? [boolean]  
navlMin = 20; % what is your threshold for # samples over the Common Era?
smoothProxies = 0 ; % do you want to smooth proxies or use raw data?
smoothExport = 0; % should we graphically export result of the smoothing procedure? [boolean]
tStart = 1; % define start year (remember: the Common Era does not have a year 0). 
tEnd   = 2000; %  define end year for the analysis
f_hr = 1/10;  % smoothing cutoff frequency for high-resolution records [in 1/years]
f_lr = 1/30; % smoothing cutoff frequency for low-resolution records  [in 1/years]


% define I/O files
f_out = ['./data/pages2k_composite_' vers '.mat'];
f_db = ['./data/pages2kTS_' vers '_unpack'];
f_temp = './data/had4med_graphem_sp70_annual'; % replaced with Cowtan & Way median
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

%% STAGE 1a: load relevant data files and prepare the data
pages2k_composite_prep


%% STAGE 1b: create global binned composites, and test their sensitivity to methodological choices
pages2k_globalBins

%% STAGE 2: sensitivity analysis of the composite

% proxy matrix retained (purpose: comparison between raw and processed data)
if smoothProxies
    proxy_r = standardize(proxy_f(:,idx_q));   
    smoothString = 'Smoothed';
else
    proxy_r = proxy_sgn(:,idx_q);
    smoothString = 'Unsmoothed';
end

% stratification by archive type
pages2k_compositeByArchive
% stratification by latitude band
pages2k_compositeByLatitudeBand
% sensitivity to record length
pages2k_composite_recordLength

% stratification however you want !
  % write your own and share it with the PAGES2k community
  

%% STAGE 3: temperature reconstruction per se

% TO DO: fix # of records
pages2k_composite_regression





