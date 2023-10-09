# Energy_Balance_Climate_Feedback
Diagnosing climate feedback from latitudinal variation in observed climatology

To use this code, first download all code and files. 

Also download: the CRUTEM absolute temperature record (available https://crudata.uea.ac.uk/cru/data/temperature/); 
the HadCRUT5 temperature anomaly dataset (available here: version HadCRUT.5.0.1.0: https://www.metoffice.gov.uk/hadobs/hadcrut5/data/current/download.html); 
the CERES EBAF Edition 4.1 satellite observational dataset for outgoing radiation (available here: https://ceres.larc.nasa.gov/data/); 
The CLARA v2.1 dataset providing cloud amount data, (available at: https://doi.org/10.24381/cds.68653055); 
the climatological annual- and zonal-mean surface relative humidity monthly climatology for 1991 to 2000 from the European Centre for Medium Range Weather Forecasts (ECMWF) ERA5 reanalysis product, as supplied through the Copernicus Climate Change Service (C3S) through the ‘Essential climate variables for assessment of climate variability from 1979 to present’ (avalailbe at: https://cds.climate.copernicus.eu/cdsapp#!/dataset/ecv-for-climate-change?tab=overview ; downloaded 31st March 2023); 
and the code for calculating the height of the tropopause, from Mateus et al. (2022) (available for download here: https://github.com/pjmateus/global_tropopause_model); where this study uses the options for a bilinear interpolation model of the tropopause, and a surface at 3.0 potential vorticity units, where 1 potential vorticity unit is equal to 10-6 K kg-1 m2 s-1. 

(1) Keep all files ending “.m” and the “HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc” in a main folder
From that main folder:
(2) Place the 15 CMIP6 files (ending .nc) provided into a new folder named “abrupt4xCO2kerneldecomposition”
(3) Place the CLARA v2.1 cloud data files you have downloaded separately into a new folder named “Cloud_data”. The files should be named “CFCmm20181101000000219AVPOSE1GL.nc” to CFCmm19820201000000219AVPOS01GL.nc”
(4) In Matlab, run the script “lambda_zonal_observations.m” by typing “lambda_zonal_observations”.


Note that to perform the analysis the observational data must be separately downloaded from the sources above. The analysis in this study is conducted by extracting 5° latitudinal resolution data from each dataset, to match the resolution of the surface temperature datasets.

