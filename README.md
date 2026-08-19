# Large-sample synthetic hydrology and forecasting with skill modification and FIRO implications codebase
This repository contains the codebase to generate synthetic hydrologic sequences from Stochastic Weather Generator (SWG) inputs, simulate synthetic forecasts against these sequences using the methodology from Brodeur et al. (2025), impart skill modifications to the generated synthetic forecasts (skill improvements or degradations), and evaluate a Forecast Informed Reservoir Operations (FIRO) strategy against these synthetic sequences of hydrology and forecasts. The dataset and model uses a case-study of Lake Oroville, CA and supports the manuscript currently (as of June 2026) in review at Water Resources Research: 'A Large-Sample Framework for Evaluating Forecast Informed Reservoir Operations (FIRO) Under Future Hydroclimatic Extremes and Evolving Forecast Skill '
   
Brodeur, Z. P., Taylor, W., Herman, J. D., & Steinschneider, S. (2025). Synthetic Ensemble Forecasts: Operations‐Based Evaluation and Inter‐Model Comparison for Reservoir Systems Across California. Water Resources Research, 61(e2024WR039324). https://doi.org/10.1029/%25202024WR039324   
   

The workflow described below utilizes code (stored in this repository and archived via formal releases at link below) and data (stored in Zenodo repository) to:   
1. Generate synthetic hydrology for Lake Oroville (Location ID: YRS, Site ID: ORDC1) from 1000-year hydrologic sequences developed from SWG simulations and a process-based hydrologic model (SAC-SMA). The synthetic hydrology process is stochastic and generates a total of 100 samples that are used in various parts of the analysis _Note: The SWG and SAC-SMA models were developed in previous work and described in the draft manuscript_
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
### 1. Data Processing 
Note: Current settings of the scripts support processing specific to the 'YRS' location and 'ORDC1' site; ORDC1 will be reference here on out. 
1) ./src/data_processing.R
   - Prepares the raw Hydrologic Ensemble Forecast Service (HEFS) ensemble forecast .csv's for the ORDC1 site to be processed by follow-on routines
   
### 2. Synthetic Hydrology Generation
1) ./src/optimize_synthetic_forecasts.R
   - main optimization script; calls the 'synthetic_forecast_opt-fun.R' function
2) ./src/synthetic_forecast_opt-fun.R
   - function used for the optimization; calls the 'syn_gen_opt.R' function
3) ./src/syn_gen_opt.R
   - synthetic generation function used by the optimization procedure.
  
### 3. Synthetic Forecast Generation and Skill Modification
1) ./src/optimize_synthetic_forecasts.R
   - main optimization script; calls the 'synthetic_forecast_opt-fun.R' function
2) ./src/synthetic_forecast_opt-fun.R
   - function used for the optimization; calls the 'syn_gen_opt.R' function
3) ./src/syn_gen_opt.R
   - synthetic generation function used by the optimization procedure.

### 4. Training and Simulation of FIRO model
Note: Must specify the 'loc' variable for the overall location and the 'keysite_name' variable for the site used to condition the kNN sampling. Other user defined parameters are set to the default used in the study, including '5fold-test' for the fully out-of-sample 5-fold procedure 
The user needs to run the following scripts in this order for the model to produce the synthetic forecasts. _Note that the ./src/create_synthetic_forecasts.R script calls the function ./src/syn_gen.R, which holds the actual synthetic forecast model_:
1) ./src/create_synthetic_forecasts.R
2) ./src/syn_gen.R

## 5. Analysis and Plotting

Detailed forecast verification for the synthetic forecasts is located at this repo: [https://github.com/zpb4/Synthetic-Forecast_Verification](https://github.com/zpb4/Synthetic-Forecast_Verification). This repo is designed to integrate seamlessly with the Synthetic Forecast generation output if located on the same root directory.


## Contact

Zach Brodeur: zpb4@cornell.edu
