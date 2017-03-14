clear all; close all;
addpath(genpath('../utilities')) % load code
vers = '2.0.0';   % database version

% global options
options.InterpSuperAnn = 1; % re-interpolate super-annual proxies or not?
options.method = 'isospectral'; %'isospectral'; %'knockoff';
options.dtype  = 'noDetrend'; %'noDetrend'; % detrending method for correlation calculation
%  'detrend' : linear detrending,
%  'diff1'   : first-differencing the original time series,
%  'noDetrend': no detrending
options.norm_p = 0; % use normalized proxies (1) or not (0)

% screening/correlations options
options.sample_thresh = 20; %minimum number of samples for a correlation to be meaningful
options.nsim = 1000; % # of surrogate timeseries in non-parametric tests
options.search_radius  = 2000; % in km


% ================================================
% 1. DATABASE UNPACKING, PLOTTING AND PRE-PROCESSING
% ================================================

% 1.1 unpack database into a more usable form and annualize data
visual = 0; % (set to 1 for debug purposes)
pages2k_db_unpack(vers, visual);

% 1.2  map it all
export = 0;
pages2k_db_synopsis(vers, export);

% 1.3 Data pre-processing
options.export = 0;
pages2k_db_process(vers,options);

% 1.4 Screen for temperature predictors

% name of correlation file
if options.norm_p
    corr_file = ['../data/PAGES2k_hadcrut4_corr_' options.dtype '_' vers '_normal_' options.method];
else
    corr_file = ['../data/PAGES2k_hadcrut4_corr_' options.dtype '_' vers '_raw_' options.method];
end

if ~exist([corr_file, '.mat'],'file')  % Only run this if necessary
    pages2k_db_screen(vers,options);
end
%
% % grid centric analysis
pages2k_grid_analysis(vers,options)

% ================================================
% 2. QC figures
% ================================================
options.export  = 1;
options.erase   = 0; %if true, clears existing figures.
options.n_neigh = 10; %number of proxy neighbors to compare to (if no temperature calibration possible)
pages2k_db_qc_plots(vers,options);

% ================================================
% 3. EXPORT
% ================================================
% export all necessary matrices for reconstructions and SD figures
pages2k_db_prep(vers,options);
