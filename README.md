
# PAGES2k SD code
Matlab code for the PAGES2k synthesis of temperature-sensitive proxies, version 2.0.0.
The code and data herein allow to reproduce all PAGES2k-related figures from the main text and SOM. (though you might have to dig through the figs directory to find exactly what you want).

## Input data
the code needs the files `PAGES2k_v2.0.0.mat` and `had4med_graphem_sp70.mat` from the FigShare repository in the **data/** directory. (or at least, an alias pointing to such files)

## Directory code_db

This code takes  `PAGES2k_v2.0.0.mat`, and `had4med_graphem_sp70.mat` and produces several synopsis plots, including all the quality control plots (Figs 4 and 5 of the main paper, as well as components of all the "QCfig_bundles" PDF files).
The master code, **pages2k_db_workflow.m** is made up of the following pieces:

- **pages2k_db_unpack.m** unpacks the database, annualizes subannual data (mostly corals) and forms data matrices
- **pages2k_db_synopsis.m**

## Directory code_validation
This code assumes that **pages2k_db_workflow.m** has run to completion and produced a file called `pages2k_hadcrut4_noDetrend_raw _2.0.0.mat`, placed in the /data/ folder. 

% LOAD  merged proxy /temperature data + output of correlation analyses
load(['../data/pages2k_hadcrut4_' d_str '_' g_str '_' vers])


The master code, **pages2k_SD_workflow.m** is made up of the following pieces:
- **pages2k_SD_prep.m** prepares the data for analysis. The following codes then slice and dice the living daylight out of them:

- **pages2k_SDmakeFig01** makes the components of Fig 01 (it needs to be assembled in Illustrator or a similar software to be publication-quality).
- **pages2k_SDmakeFig02** same for Fig 02
- **pages2k_composite_globalBins**     create global binned composites, makes SD Fig 03
- **pages2k_compositeByArchive**       stratification by archive type, makes SD Fig 04

Al the data are in "data". The figures are all exported to "figs"
Code dependencies needed to run the above codes are included in "utilities". If you are missing anything, we will be happy to add it. To comply with the license, some third-party packages may need to be installed separately:

- [m_map](http://www.eos.ubc.ca/~rich/map.html)
