# Variable selection for balance regression

This repository contains the R code for reproducing the results in our paper titled ``Variable selection for balance regression with applications to microbiome compositional data'' which are organized in three folders:
* Functions/: Functions for various methods, e.g., our supervised log-ratio method, and utility functions for their comparison.
* Simulations/: Simulated data and results from comparing the methods on them.
* Data/: Results from comparing the methods on two case studies.

If you are interested in fitting our supervised log-ratios method for your own analyses, please see the GitHub repository [drjingma/slr](https://github.com/drjingma/slr) for instructions on how to install the most updated version of our R package.

## Functions

This folder contains:
* `slrs.R`: An implementation of our supervised log-ratio method.
* `util.R`: Helper functions that are used to facilitate data simulation and the comparison of the methods on metrics of interest.
* `codalasso.R`: An implementation of codalasso (i.e., compositional lasso for classification problems), borrowed from this [GitHub repository](https://github.com/cunningham-lab/codacore) which incorporates cross-validation. We edited their original code to enable cross-validation based on AUC.

## Experiments

This folder contains:
* R scripts for generating the simulated data. 
    * `sim_data_IBD.R`: simulate data based on the IBD study
    * `sim_data_CRC.R`: simulate data based on the CRC study
* R scripts for comparing different methods on the simulated data. 
    * `run_IBD.R`: run different methods on data simulated based on the IBD study
    * `run_CRC.R`: run different methods on data simulated based on the CRC study
* R scripts for summarizing the results: 
    * `metrics.R`: need to adjust the `condition` argument to toggle between different 
settings
* Parameters/: a folder where each file contains the parameters (`X`, `X.test`, `sbp`) used in one simulation setting.

## Data

This folder contains:
* `Franzosa_PRISM_UC.rda` and `Franzosa_Validation_UC.rda`: processed IBD data after removing rare taxa. For raw counts, please see [here](https://github.com/borenstein-lab/microbiome-metabolome-curated-data/tree/main/data/processed_data/FRANZOSA_IBD_2019).  
* `Franzosa_IBD_preprocessing.R`: R script for preprocessing the IBD data. 
* `Franzosa_IBD.R`: R script for analyzing the IBD data.
* `Yachida_CRC_2019.rda`: processed CRC data after removing rare taxa. For raw counts, please see [here](https://github.com/borenstein-lab/microbiome-metabolome-curated-data/tree/main/data/processed_data/YACHIDA_CRC_2019).
* `Yachida_CRC_preprocessing.R`: R script for preprocessing the CRC data.
* `Yachida_CRC.R`: R script for analyzing the CRC data.
* `metrics_data.R`: R script for summarizing the results; need to adjust the `dataname` argument to toggle between different data sets. 

 
