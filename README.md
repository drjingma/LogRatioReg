# Supervised Log-Ratios

This repository contains the R code for reproducing the results in our paper titled ``Regression and Classification of Compositional Data via a novel Supervised Log Ratio Method,'' which are organized in three folders:
* Functions: Functions for various methods, e.g., our supervised log-ratio method, and utility functions for their comparison.
* Experiments: Simulated data and results from comparing the methods on them.
* Data: Results from comparing the methods on three real data sets.

If you are interested in fitting our supervised log-ratios method for your own analyses, please see the GitHub repository [drjingma/slr](https://github.com/drjingma/slr) for instructions on how to install the most updated version of our R package.

## Functions

This folder contains:
* `slrs.R`: An implementation of our supervised log-ratio method.
* `compositionallasso.R`: An implementation of compositional lasso for regression problems.
* `codalasso.R`: An implementation of codalasso (i.e., compositional lasso for classification problems), borrowed from this [GitHub repository](https://github.com/cunningham-lab/codacore) which incorporates cross-validation. We edit this code to perform cross-validation on AUC.
* `logratiolasso.R`: An implementation of log-ratio lasso borrowed from this [GitHub repository](https://github.com/stephenbates19/logratiolasso) which we have edited to be able to use the 1-standard error rule in model selection via cross-validation.
* `util.R`: Helper functions that are used to facilitate data simulation and the comparison of the methods on metrics of interest.

## Experiments

This folder contains:
* R scripts for generating the simulated data and comparing the methods on each of them by computing various metrics. These metrics results are saved in the outputs subfolder. 
    * `sims_latent.R`: simulated data with continuous response, generated from independent covariates.
    * `sims_latent_binary.R`: simulated data with binary response, generated from independent covariates.
    * `sims_latent_correlated.R`: simulated data with continuous response, generated from correlated covariates using an error term with exponential decay.
    * `sims_latent_binary_correlated.R`: simulated data with binary response, generated from correlated covariates using an error term with exponential decay.
* R scripts for plotting the results that compare the methods on the four simulation settings above: `plot_metrics_latent.R`, `plot_metrics_latent_binary.R`, `plot_metrics_latent_correlated.R`, and `plot_metrics_latent_binary_correlated.R`.

## Data

This folder contains:
* R scripts for fitting the methods on a Crohn's disease data set, an HIV data set, and an HIV inflammation data set: `Crohn_fit.R`, `HIV_fit.R`, `sCD14_fit.R`.
* R scripts for randomly splitting the three data sets into training/test sets and comparing the methods on each of them by computing various metrics. These metrics results are saved in the outputs subfolder: `Crohn_metrics.R`, `HIV_metrics.R`, `sCD14_metrics.R`.
* R scripts for plotting the results from each method and comparisons of the methods:
    * `plot_mse_oneselectionbarplot.R` creates the bar graphs that compare the selections made by balance regression models for each of the data sets.
    * `plot_mse_metricboxplots_all.R` creates the boxplots that compare the methods on metrics of interest on each of the data sets.

