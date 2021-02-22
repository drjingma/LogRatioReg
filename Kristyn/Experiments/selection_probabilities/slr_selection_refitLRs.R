# Method: Supervised Log-Ratios
# Purpose: Compute selection probabilities via bootstrap
# Date: 12/2020
# Notes: This code is parallelized. It refits the coefficients on the balances
#        using a linear model (no constraint, since over balances instead of 
#        compositions)

getwd()

# libraries
library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
library(limSolve) # for constrained lm, lsei()
library(stats) # for hclust()
library(balance) # for sbp.fromHclust()

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Dr. Ma sources
library(Matrix)
library(glmnet)
library(compositions)
library(stats)
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))

# settings
tol = 1e-4
cv.n_lambda = 200
cv.K = 10
intercept = TRUE

# Bootstrap
bs.n = 100

# data
# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
load(paste0("Data/", "BMI.rda"))
# dim(raw_data) # 98 x 89
# dim(X) # 98 x 87
# dim(X.prop) # 98 x 87
log.X.prop = log(X.prop)
n = dim(X)[1]
num.genera = dim(X)[2]

################################################################################
# Actual work: Selecting the genera #
################################################################################

# They generate 100 bootstrap samples and use the same CV procedure to select 
#   the genera (for stable selection results)

bs.selected_variables = foreach(
  b = 1:bs.n, 
  .combine = cbind
) %dorng% {
  source("RCode/func_libs.R")
  library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
  library(limSolve) # for constrained lm, lsei()
  library(stats) # for hclust()
  library(balance) # for sbp.fromHclust()
  
  # resample the data
  bs.resample = sample(1:n, n, replace = TRUE)
  X.prop.bs = X.prop[bs.resample, ]
  y.bs = y[bs.resample]
  
  # refitted CV
  # Split the data into K folds
  shuffle = sample(1:n)
  idfold = (shuffle %% cv.K) + 1
  n_fold = as.vector(table(idfold))
  
  # Do cross-validation
  # calculate squared error (prediction error?) for each fold, 
  #   needed for CV(lambda) calculation
  cvm = rep(NA, cv.n_lambda) # want to have CV(lambda)
  cvm_sqerror = matrix(NA, cv.K, cv.n_lambda)
  # Fit Lasso for each fold removed
  for (j in 1:cv.K){
    # Training data
    Xtrain = X.prop.bs[idfold != j, ]
    Ytrain = y.bs[idfold != j]
    # Test data
    Xtest = X.prop.bs[idfold == j, ]
    Ytest = y.bs[idfold == j]

    # Fit SLR (with noise, because bootstrap leads to replicated observations)
    slr_j = fitSLR(Ytrain, Xtrain, nlam = cv.n_lambda, intercept = intercept)
    # transform the coefficients for our log-ratios model (slr) on balances to 
    #   coefficients for the linear log-contrasts model on log-contrasts
    slr_betas = as.matrix(slr_j$bet)
    btree = slr_j$btree
    non0.betas = slr_betas != 0 # diff lambda = diff col
    for(m in 1:cv.n_lambda){
      selected_variables = non0.betas[, m]
      # get refitted coefficients, after model selection and w/o penalization
      if(all(!selected_variables)){ # if none selected
        refit = lm(Ytrain ~ 1)
        predictCLM = function(X) coefficients(refit)
      } else{ # otherwise, fit on selected variables
        # fit to the subsetted data
        Xbtrain.sub = computeBalances(Xtrain, btree)[, selected_variables, drop = FALSE]
        lm.fit = lm(Ytrain ~ Xbtrain.sub)
        lm.coeffs = coefficients(lm.fit)
        predictCLM = function(X){
          lm.coeffs[1] + 
            computeBalances(X, btree)[, selected_variables, drop = FALSE] %*% 
            lm.coeffs[-1]
        }
      }
      # calculate squared error (prediction error?)
      Ypred = predictCLM(Xtest)
      cvm_sqerror[j, m] = sum(crossprod(Ytest - Ypred))
    }
  }
  
  # Calculate CV(lambda) for each value of lambda
  cvm = colMeans(cvm_sqerror)
  
  # Find lambda_min = argmin{CV(lambda)}
  lambda_min_index = which.min(cvm)
  lambda_min = slr_j$lambda[lambda_min_index]
  
  # final fit
  slr_select.bs = fitSLR(y.bs, X.prop.bs, lambda = lambda_min, nlam = 1, 
                         intercept = intercept)
  XYdata = data.frame(X.prop.bs, y = y.bs)
  selected_variables = LRtoLC(as.matrix(slr_select.bs$bet), slr_select.bs$btree) != 0 # diff lambda = diff col
  # record which variables were selected
  selected_variables
}
rownames(bs.selected_variables) = colnames(X)

bs.selected_variables_numeric = apply(bs.selected_variables, 2, as.numeric)
bs.selection_percentages = apply(bs.selected_variables_numeric, 1, FUN = 
                                   function(x) sum(x, na.rm = TRUE) / bs.n)
names(bs.selection_percentages) = rownames(bs.selected_variables)
bs.results = list(
  seed = rng.seed,  
  selected_variables = bs.selected_variables, 
  selection_percentages = bs.selection_percentages
)

saveRDS(bs.results,
        file = paste0("Kristyn/Experiments/output",
                      "/slr_selection", 
                      "_refitLRs",
                      "_int", intercept, 
                      "_B", bs.n, 
                      "_seed", rng.seed,
                      ".rds"))

sort(bs.selection_percentages)
