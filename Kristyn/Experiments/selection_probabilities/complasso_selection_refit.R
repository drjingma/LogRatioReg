# Method: Compositional Lasso
# Purpose: Compute selection probabilities via bootstrap
# Date: 12/2020
# Notes: This code is parallelized.

getwd()

# libraries
library(limSolve) # for constrained lm

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

# settings
tol = 1e-4

# Cross-validation
cv.n_lambda = 200
cv.K = 10

# Bootstrap
bs.n = 500

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

# bs.finalfits = list()

bs.selected_variables = foreach(
  b = 1:bs.n, 
  .combine = cbind, 
  .noexport = c("ConstrLassoC0")
) %dorng% {
  source("RCode/func_libs.R")
  library(limSolve)
  
  # resample the data
  bs.resample = sample(1:n, n, replace = TRUE)
  log.X.prop.bs = log.X.prop[bs.resample, ]
  y.bs = y[bs.resample]
  
  # refitted CV
  # Split the data into 10 folds
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
    Xtrain = log.X.prop.bs[idfold != j, ]
    Ytrain = y.bs[idfold != j]
    # Test data
    Xtest = log.X.prop.bs[idfold == j, ]
    Ytest = y.bs[idfold == j]
    
    # Fit LASSO on that fold using fitLASSOcompositional
    Lasso_j = ConstrLasso(
      Ytrain, Xtrain, Cmat = matrix(1, dim(Xtrain)[2], 1), 
      nlam = cv.n_lambda, tol = tol)
    non0.betas = Lasso_j$bet != 0 # diff lambda = diff col
    for(m in 1:cv.n_lambda){
      selected_variables = non0.betas[, m]
      # get refitted coefficients, after model selection and w/o penalization
      if(all(!selected_variables)){ # if none selected
        refit = lm(Ytrain ~ 1)
        predictCLM = function(X) coefficients(refit)
      } else{ # otherwise, fit on selected variables
        # fit to the subsetted data
        Xtrain.sub = Xtrain[, selected_variables, drop = FALSE]
        Xtrain.sub2 = cbind(1, Xtrain.sub)
        Q = as.matrix(rep(1, sum(selected_variables)))
        Q2 = rbind(0, Q)
        colnames(Xtrain.sub2)[1] = "Intercept"
        rownames(Q2) = colnames(Xtrain.sub2)
        lsei.fit = lsei(A = Xtrain.sub2, B = Ytrain, E = t(Q2), F = 0)
        predictCLM = function(X){
          lsei.fit$X[1] +
            X[, selected_variables, drop = FALSE] %*% lsei.fit$X[-1]
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
  lambda_min = Lasso_j$lambda[lambda_min_index]
  
  # final fit
  Lasso_select.bs = ConstrLasso(
    y.bs, log.X.prop.bs, Cmat = matrix(1, dim(log.X.prop)[2], 1), 
    lambda = lambda_min, nlam = 1, tol=tol)
  XYdata = data.frame(log.X.prop.bs, y = y.bs)
  non0.betas = Lasso_select.bs$bet != 0 # diff lambda = diff col
  selected_variables = non0.betas
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
                      "/complasso_selection", 
                      "_refit",
                      "_B", bs.n, 
                      "_seed", rng.seed,
                      ".rds"))

sort(bs.selection_percentages)