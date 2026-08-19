# Large-sample synthetic hydrology and forecasting with skill modification and FIRO implications codebase
This repository contains the codebase to generate synthetic hydrologic sequences from Stochastic Weather Generator (SWG) inputs, simulate synthetic forecasts against these sequences using the methodology from Brodeur et al. (2025), impart skill modifications to the generated synthetic forecasts (skill improvements or degradations), and evaluate a Forecast Informed Reservoir Operations (FIRO) strategy against these synthetic sequences of hydrology and forecasts. The dataset and model uses a case-study of Lake Oroville, CA and supports the manuscript currently (as of June 2026) in review at Water Resources Research: 'A Large-Sample Framework for Evaluating Forecast Informed Reservoir Operations (FIRO) Under Future Hydroclimatic Extremes and Evolving Forecast Skill '
   
Brodeur, Z. P., Taylor, W., Herman, J. D., & Steinschneider, S. (2025). Synthetic Ensemble Forecasts: Operations‐Based Evaluation and Inter‐Model Comparison for Reservoir Systems Across California. Water Resources Research, 61(e2024WR039324). https://doi.org/10.1029/%25202024WR039324   
   

The workflow described below utilizes code (stored in this repository and archived via formal releases at link below) and data (stored in Zenodo repository) to:   
1. Generate synthetic hydrology for Lake Oroville (Location ID: YRS, Site ID: ORDC1) from 1000-year hydrologic sequences developed from SWG simulations and a process-based hydrologic model (SAC-SMA). The synthetic hydrology process is stochastic watershed modeling (SWM) procedure and generates a total of 100 samples that are used in various parts of the analysis _Note: The SWG and SAC-SMA models were developed in previous work and described in the draft manuscript_
2. Fit and generate synthetic forecasts for Lake Oroville against multiple samples of the synthetic hydrologic sequences
3. Post-process the synthetic forecasts to impart differing degrees of skill modification
4. Train a FIRO model to various combinations of synthetic hydrology and forecasts
5. Evaluate FIRO model against out-of-sample synthetic hydrology and forecasts
6. Plot results to support the formal analysis in the draft manuscript




Raw data to support this code are stored here:  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20801662.svg)](https://doi.org/10.5281/zenodo.20801662)  
Releases of this software are stored permanently here:   
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20803036.svg)](https://doi.org/10.5281/zenodo.20803036)   

---
#### Note: Once data have been downloaded from Zenodo repository and unzipped, the entire contents of the folder should be placed in a folder labeled 'data' in the root repository. A repository labeled 'out' will archive all the generated data.

## Dependencies
- R package 'future'
- R package 'DEoptim'
- R package 'twosamples'
- R package 'ncdf4'
- Python library 'numpy'
- Python library 'pandas'
- Python library 'xarray'
- Python library 'scipy.optimize'
- Python library 'numba'
   
## Workflow
The workflow is arranged in consecutive steps, where in many cases the data from the previous step is required to run the current step.
Note: Current settings of the scripts support processing specific to the 'YRS' location and 'ORDC1' site; ORDC1 will be reference here on out. 
### 1. Data Processing 
Clean and process the ensemble forecast data needed to support the synthetic forecasting routine
1) ./src/data_processing.R
   - Prepares the raw Hydrologic Ensemble Forecast Service (HEFS) ensemble forecast .csv's for the ORDC1 site to be processed by follow-on routines
   
### 2. Synthetic Hydrology Generation
Generates the 1000-year Stochastic Watershed Model (SWM) hydrologic sequences
1) ./src/swm_gen.R
   - generates stochastic, synthetic hydrologic sequence samples based on the hydrologic simulations from a couple SWG and SAC-SMA model; routine generates synthetic hydrology under both historical conditions and +4C warmed conditions in the SWG. 
  
### 3. Synthetic Forecast Generation and Skill Modification
Fit and simulate synthetic ensemble forecasts of differing skill levels against the synthetic hydrology sequences
1) ./src/optimize_synthetic_forecasts.R   
   - optimization routine for the synthetic forecasting model against the HEFS training data; calls the 'synthetic_forecast_opt-fun.R' function
2) ./src/create_synthetic_forecasts_swm.R
   - generates synthetic forecasts against both the the historical and 4C warmed synthetic hydrology scenarios; calls the 'syn_gen_swm.R' function
3) ./src/calc-climo-forecast.R
   - calculates a climatological forecast for each synthetic hydrology sequence to enable the skill degradation procedure 
4) ./src/gen_skill-mod_swm.R   
   - postprocesses generated synthetic forecasts to impart differing degrees of specified skill improvement between 0 and 1   
5) ./src/gen_skill-mod-negative_swm.R   
   - postprocesses generated synthetic forecasts to impart differing degrees of specified skill degradation between 0 and -1

#### 3a. Synthetic forecast helper functions
   - ./src/synthetic_forecast_opt-fun.R   
      - function used for the optimization; calls the 'syn_gen_opt.R' function
   - ./src/syn_gen_opt.R   
      - synthetic generation function used by the optimization procedure
   - ./src/syn_gen_swm.R   
      - synthetic generation function used to generate synthetic forecasts against the synthetic hydrologic sequences

### 4. Training and Simulation of FIRO model
Fit and simulate a FIRO model for Lake Oroville against different sequences of synthetic hydrology and forecasts
1) ./src/train_param_skill-mod_swm.py

#### 4a. FIRO model helper functions
   - ./src/model.py   
      - Lake Oroville FIRO simulation model
   - ./src/syn_util.py  
      - various helper functions for the FIRO simulation routine
   - ./src/util.py  
      - additional helper functions for the FIRO simulation routine

### 5. Analysis and Plotting
Simulate and analyze FIRO outcomes, calculate metrics of interest, and plot results
1) ./src/calc_ecrps-rankhist_swm.R   
   - calculates eCRPS and Rank Histogram metrics for synthetic forecast simulations for verification purposes
2) ./src/calc_skill-mod_display-stats.R  
   - calculates ensemble forecast metrics for synthetic forecast simulations to support plotting analyses
3) ./src/calc_swm-topx-metrics.py   
   - simulates FIRO model against synthetic hydrology and forecast sequences and evaluates outcomes in an aggregated sense to support plotting analyses
4) ./src/calc_swm-topx-evts-metrics.py   
   - simulates FIRO model against synthetic hydrology and forecast sequences and evaluates outcomes in an event-based manner to support plotting analyses
5) ./plot/*
   - plotting functions once previous analysis is complete arranged in the order presented in the draft manuscript

#### 5a. Analysis and plotting helper functions
   - ./src/forecast_verification_functions.R 
      - helper functions for the forecast verification routines written in R
   - ./src/ensemble_verification_functions.py
      - helper functions to support evaluation of FIRO outcomes
   - ./src/mm-cfs_conversion.R
      - initial conversion from baseline simulation in mm to cfs

## Contact

Zach Brodeur: zpbrodeur@ucsd.edu
