# Large-sample synthetic hydrology and forecasting with skill modification and FIRO implications codebase
This repository contains the codebase to generate synthetic hydrologic sequences from Stochastic Weather Generator (SWG) inputs, simulate synthetic forecasts against these sequences using the methodology from Brodeur et al. (2025), impart skill modifications to the generated synthetic forecasts (skill improvements or degradations), and evaluate a Forecast Informed Reservoir Operations (FIRO) strategy against these synthetic sequences of hydrology and forecasts. The dataset and model uses a case-study of Lake Oroville, CA and supports the manuscript currently (as of June 2026) in review at Water Resources Research: 'A Large-Sample Framework for Evaluating Forecast Informed Reservoir Operations (FIRO) Under Future Hydroclimatic Extremes and Evolving Forecast Skill '
   
Brodeur, Z. P., Taylor, W., Herman, J. D., & Steinschneider, S. (2025). Synthetic Ensemble Forecasts: Operations‐Based Evaluation and Inter‐Model Comparison for Reservoir Systems Across California. Water Resources Research, 61(e2024WR039324). https://doi.org/10.1029/%25202024WR039324   
   

The workflow described below 



---
#### Note: After downloading and extracting data from Hydroshare resources above, ensure local directory path for HEFS data is configured: './Synthetic-Forecast-v2-FIRO-DISES/data/_main_hindcast_location_/...', where '...' are the site specific sub-repos defined in 'Data' section below. Unzipping the files can result in duplication in the data path and this must be corrected for the code to function.

## Dependencies
- R package 'stringr'
- R package 'lubridate'

   
Information below describes setup and execution of the model:   
## Data


## Workflow
### 1. Data Processing 
Note: Specification of the 'loc' variable in these scripts will process all associated sites
Processes the forecast and observation data to standard R data structures for Period of Record (POR) hindcast:
1) ./src/data_processing.R
   - all HEFS and observations data files for POR
   
### 2. Optimization
Note: Must specify the 'loc' variable for the overall location and the 'keysite_name' variable for the site used to condition the kNN sampling. Other user defined parameters are set to the default used in the study, including '5fold-test' for the fully out-of-sample 5-fold procedure outlined in the manuscript.  
Optimizes the threshold curve that constrains the scaling procedure for synthetic forecast generation:
1) ./src/optimize_synthetic_forecasts.R
   - main optimization script; calls the 'synthetic_forecast_opt-fun.R' function
2) ./src/synthetic_forecast_opt-fun.R
   - function used for the optimization; calls the 'syn_gen_opt.R' function
3) ./src/syn_gen_opt.R
   - synthetic generation function used by the optimization procedure.

### 3. Generation
Note: Must specify the 'loc' variable for the overall location and the 'keysite_name' variable for the site used to condition the kNN sampling. Other user defined parameters are set to the default used in the study, including '5fold-test' for the fully out-of-sample 5-fold procedure 
The user needs to run the following scripts in this order for the model to produce the synthetic forecasts. _Note that the ./src/create_synthetic_forecasts.R script calls the function ./src/syn_gen.R, which holds the actual synthetic forecast model_:
1) ./src/create_synthetic_forecasts.R
2) ./src/syn_gen.R

The output of the first two steps is an R array that is saved as an R data structure file (.rds). In order to further post-process data for transfer to other models, languages, etc, there are two output options:   

3) ./src/5fold-val_collate.R
   - combines the folds to a single array for the 5-fold generation schemes ('5fold', '5fold-test')
4) ./src/data_writeout_ncdf.R
   - writes both HEFS and synthetic HEFS files to a netCDF file
5) ./src/slice_plot-ens.R
   - slices a 10x sample subset from the generated synthetic forecast array for plotting; the raw arrays for large sample runs (e.g. 100 samples) require too much RAM for typical personal computers

All scripts create and output metadata to the ./out/_main_hindcast_location_/ subdirectory. For sites with separate 1986 data, there are separate scripts with a '_86.R' suffix to process those specific data subsets.  

## 4. Verification

Detailed forecast verification for the synthetic forecasts is located at this repo: [https://github.com/zpb4/Synthetic-Forecast_Verification](https://github.com/zpb4/Synthetic-Forecast_Verification). This repo is designed to integrate seamlessly with the Synthetic Forecast generation output if located on the same root directory.


## Contact

Zach Brodeur: zpb4@cornell.edu
