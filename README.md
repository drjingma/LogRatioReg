# About

This repo provides code, data, and references for the supervised log-ratio project. 

## To-dos 04/18

- When variables are correlated, two iterations of spectral clustering may not correctly identify active variables. There may be too many false positives. Explore the correlation matrix in several real data examples and suggest ways of determining the number of binary splits needed. Can we use cross-validation to determine the number of binary splits?

- Inform a more realistic simulation design based on the correlation observed in real data.

- Consider an example with 2 latent variables.

- Regularization helps in spectral clustering using K-means. Explore whether different regularization can help deal with the unbalanced balance case.  
        
        * Current results suggest that high-degree regularization does not seem to help with unbalanced balance regression, but could this be due to the choice of weights? How do the current weights impact the similarity matrix? 
        * Uniform regularization to all entries of the similarity matrix and low-rank approximation together yields better performance than using either alone.

- Spectral clustering is only optimal on population similarity matrix. What is the most appropriate low-rank regularization? Should we apply low-rank approximation to the Laplacian instead?

## Other Issues

- How to address the unbalanced balance issue? 

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