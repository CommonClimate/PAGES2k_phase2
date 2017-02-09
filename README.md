
# PAGES2k SD code
Matlab code for the PAGES2k synthesis of temperature-sensitive proxies, version 2.0.0.
The code and data herein allow to reproduce all PAGES2k-related figures from the main text and SOM.
The code is provided "as is" and without warranty of running on your system. We welcome any and all feedback to make it run more smoothly.

## Input data
the code needs the files `PAGES2k_v2.0.0.mat` and `had4med_graphem_sp70.mat` from the FigShare repository in the **data/** directory. (or at least, an alias pointing to such files)

## Directory code_db

The code in this directory takes  `PAGES2k_v2.0.0.mat`, and `had4med_graphem_sp70.mat` and produces several synopsis plots, including all the quality control plots (Figs 4 and 5 of the main paper, as well as components of all the "QCfig_bundles" PDF files).
In most cases, additional visualizations are made if certain flags (visual, export) are set to 1 (TRUE).

The master code, **pages2k_db_workflow.m** is made up of the following pieces:

- **pages2k_db_unpack.m** unpacks the database, annualizes subannual data (mostly corals) and forms data matrices. if visual ==1, plots the visualizations in figs/annualize/

- **pages2k_db_synopsis.m** defines graphical attributes and makes a few synopsis plots (under figs)

- **pages2k_db_process.m** tests data for normality, applies inverse transform sampling to the data matrix `proxy`, so each column of `proxy_n` is distributed as a standard normal). If the flag InterpSuperAnn is set, the code also interpolates super-annual proxies to annual resolution. (proxy -> proxy_a, and proxy_n -> proxy_na)

- **pages2k_db_screen.m** carries out the record-based correlation analysis used in QC figures, as well as Fig 6a-c.
(NB: if nsim is large, this part of the code may take many hours to run, even on a decent computer. If its output file already exists in the data directory, this step will be skipped).

- **pages2k_grid_analysis.m** carries out the grid-based correlation analysis shown in Fig 6d.

- **pages2k_db_qc_plots.m** produces the QC figures, based on the saved output of **pages2k_db_screen.m**.

- **pages2k_db_prep** % exports all necessary matrices for subsequent analyses.

## Directory code_HadCRUT

- **HadCRUT4_prep.m** annualizes monthly HadCRUT4 data and generates Fig 3
- **HadCRUT4_seasonal.m** generates Fig S6

## Directory code_validation
This code assumes that **pages2k_db_workflow.m** has run to completion and produced a file called `pages2k_hadcrut4_noDetrend_raw _2.0.0.mat`, placed in the /data/ folder. (note: for this to happen, the code must have been run under default settings)

The master code, **pages2k_SD_workflow.m** is made up of the following pieces:

- **pages2k_SD_prep.m** prepares the data for analysis.
- **pages2k_SDmakeFig01** makes the components of Fig 01 (it needs to be assembled in Illustrator or a similar software to be publication-quality).
- **pages2k_SDmakeFig02_parts** maps the PASGES2k phase 2, M08 and PAGES2k phase 1 databases. The latter 2 form Fig 2.
- **pages2k_SDmakeFig06** assembles Fig 6 based on the parts computed in **pages2k_db_screen.m** and **pages2k_grid_analysis.m**
- **pages2k_composite_globalBins** creates global binned composites, makes Fig 7
- **pages2k_compositeByArchive**  stratification by archive type, makes Fig 8

All the data are in "data". All figures are exported to "figs" (sometimes in subfolders thereof).
Code dependencies needed to run the above codes are included in "utilities". If you are missing anything, we will be happy to add it. To comply with the license, some third-party packages may need to be installed separately:

- [m_map](http://www.eos.ubc.ca/~rich/map.html)
