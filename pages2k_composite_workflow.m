clear all; close all;  % clean slate
vers = '2015_08_20';  % version of the database to be use
addpath(genpath('./utilities'))% load code utilities 
% define graphical defaults
set(0, 'DefaultAxesFontName', 'Helvetica','DefaultAxesFontWeight','bold')
set(0, 'DefaultTextFontName', 'Helvetica','DefaultTextFontWeight','bold')  
set(0,'defaultaxesfontsize',12); set(0,'defaulttextfontsize',12);

% define output file
fout = ['./data/pages2k_composite_' vers '.mat'];

% Now let the wild rumpus begin

%% STAGE 1: load relevant data files and prepare the data
pages2k_composite_prep

%% STAGE 2: sensitivity analysis of the composite
options.source = 'smoothed'; % 'unsmoothed' or 'smoothed'
pages2k_composite_sensitivity2

%% STAGE 3: temperature reconstruction per se
