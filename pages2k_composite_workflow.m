clear all; close all;  % clean slate
vers = '2015_08_20';  % version of the database to be use
addpath(genpath('./utilities'))% load code utilities 
% define graphical defaults
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultTextFontWeight','bold')  
set(0,'defaultaxesfontsize',12); set(0,'defaulttextfontsize',12);

% define output file
fout = ['./data/pages2k_composite_' vers '.mat'];

% define options
lat_weight = 0; % are we normalizing by the cosine of latitude? [boolean]
sifting_style = 'qcScreenHR'; % possible choices: noSift, qcOnly, qcScreenHR, qcScreenLR, qcScreenAll
norm_p = 1;  % do you want proxies Gaussianized ?    [boolean] 
detrend = 0; % do you want to detrend coral d18O proxies? [boolean]  
navlMin = 20; % what is your threshold for # samples over the Common Era?
smoothProxies = 1 ; % do you want to smooth proxies or use raw data?
smoothExport = 0; % should we graphically export result of the smoothing procedure? [boolean]
tStart = 1; % define start year (remember: the Common Era does not have a year 0). 
tEnd   = 2000; %  define end year for the analysis
f_hr = 1/10;  % smoothing cutoff frequency for high-resolution records
f_lr = 1/30; % smoothing cutoff frequency for low-resolution records




% Now let the wild rumpus begin

%% STAGE 1: load relevant data files and prepare the data
pages2k_composite_prep

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

pages2k_composite_regression





