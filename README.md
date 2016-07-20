
# PAGES2k SD code
Data and Matlab code for the latest version of the PAGES2k synthesis of temperature-sensitive proxies.
The code and data herein allow to reproduce all figures in the main text, as well as all the supplementary figures having to do with the composite, plus Figure S2.

## Organization

The master code, **pages2k_SD_worflow.m** is made up of the following pieces:
- **pages2k_SD_prep.m** prepares the data for analysis. The following codes then slice and dice the living daylight out of them:

- **pages2k_SDmakeFig01** makes the components of Fig 01 (it needs to be assembled in Illustrator or a similar software to be publication-quality).
- **pages2k_SDmakeFig02** same for Fig 02
- **pages2k_SDmakeTable**              outputs a csv file, which can easily be converted in Excel to Table 1.
- **pages2k_composite_globalBins**     create global binned composites, makes SD Fig 03
- **pages2k_compositeByArchive**       stratification by archive type, makes SD Fig 04
- **pages2k_compositeByLatitudeBand**  stratification by latitude band, makes SD Fig S8
- **pages2k_compositeByScreenCrit**    sensitivity to Screening Criterion, makes SD Fig S4
- **pages2k_composite_ByRecordLength**  sensitivity to record length, makes SD Fig S3

Al the data are in "data". The figures are all exported to "figs"
Code dependencies needed to run the above codes are included in "utilities". If you are missing anything, we will be happy to add it. To comply with the license, some third-party packages may need to be installed separately:

- [m_map](http://www.eos.ubc.ca/~rich/map.html)
- [latextable](http://www.mathworks.com/matlabcentral/fileexchange/44274-latextable)
