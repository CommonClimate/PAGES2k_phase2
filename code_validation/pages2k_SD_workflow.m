clear all; close all;  % clean slate
vers = '2.0.0';  % version of the database to be use
addpath(genpath('../utilities'))% load code utilities

% define analysis options
lat_weight = 0; % are we normalizing by the cosine of latitude? [boolean]
sifting_style = 'qcOnly'; % possible choices: noSift, qcOnly, qcScreenHR, qcScreenLR, qcScreenAll
screenHR_style = 'reg';  % 'loc' = local; 'reg' = regional (within 2000km radius), 'fdr' = regional accounting for false discovery rate

norm_p = 1;  % should proxies be mapped to a standard normal? [boolean]
detrend = 0; % do you want to detrend coral d18O proxies? [boolean]
navlMin = 20; % what is your threshold for # samples over the Common Era?
tStart = 1; % define start year (remember: the Common Era does not have a year 0).
tEnd   = 2000; %  define end year for the analysis

% define I/O files
f_out = ['../data/pages2k_composite_' vers '.mat'];

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
load(['../data/pages2k_hadcrut4_' d_str '_' g_str '_' vers '.mat'])

% define graphical defaults
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultAxesFontWeight','normal')
set(0,'defaultaxesfontsize',12); set(0,'defaulttextfontsize',12);

%% prep work
pages2k_SD_prep                  % makes SD Fig S2

% Make Figure 1 (database view)
pages2k_SDmakeFig01
% Make parts for Figure 2 (comparison to previous efforts)
network = 'PAGES2k 2016';  % choices: M08, PAGES2k 2013, PAGES2k 2016. 
pages2k_SDmakeFig02_parts(network,vers)
% Make Figure 6 (correlation maps) 


%% create global binned composites (Fig 7)
pages2k_composite_globalBins      

%% stratification by archive type (Fig 8)
pages2k_compositeByArchive     

