# About

This repo provides code, data, and references for the supervised log-ratio project. 

## To-dos 02/24

- Explore the correlation matrix in several real data examples and suggest ways of determining the number of clusters in step 1. If we need to select the number of clusters, what is the best way of doing so? Cross-validation?

## Issues 02/24

- How to address the unbalanced balance issue? 
- The number of clusters in the first step of spectral clustering need to be first determined using some optimization criterion, because when variables in the active set are correlated with variables in the inactive set, we may also observe similarity. However, these will not be as strong as the true active variables. We can visualize the similarity and determine what is the best way of partitioning the variables in the first step. 

**RCode**

This folder consists of the code used in simulations.


**COAT-master**

This folder consists of code for implementing the COAT method (see [here](https://doi.org/10.1080/01621459.2018.1442340) for the paper).

**Refs**

A few relevant papers:

  - Pawlowsky-Glahn, Vera, Juan José Egozcue, and Raimon Tolosana Delgado. "Principal balances." (2011).
  - Lin, Wei, et al. "Variable selection in regression with compositional covariates." Biometrika 101.4 (2014): 785-797.
  - Rivera-Pinto, J., et al. "Balances: a new perspective for microbiome analysis." mSystems 3.4 (2018).
  - Bates, Stephen, and Robert Tibshirani. "Log‐ratio lasso: Scalable, sparse estimation for log‐ratio models." Biometrics 75.2 (2019): 613-624.
  - Quinn, Thomas P., and Ionas Erb. "Interpretable Log Contrasts for the Classification of Health Biomarkers: a New Approach to Balance Selection." mSystems 5.2 (2020).
  - Gordon-Rodriguez, et al. "Learning sparse log-ratios for high-throughput sequencing data." Bioinformatics. 38(1):157--163 (2022) 

**Data**

The first data set comes from the following papers. The response variable is BMI, and predictors are microbiome otu data obtained with 16S. There are 98 samples and 87 variables. Note the predictors (`X`) are sparse. In addition, the total number of reads varies greatly across samples. To remove the artificial differences in total reads, we transformed the raw counts into compositional data after replacing zero counts by the maximum rounding error 0.5. This is a common practice as seen in Lin et al. (2014).  

  - Wu, Gary D., et al. "Linking long-term dietary patterns with gut microbial enterotypes." Science 334.6052 (2011): 105-108.
  - Lin, Wei, et al. "Variable selection in regression with compositional covariates." Biometrika 101.4 (2014): 785-797.