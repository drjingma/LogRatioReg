# workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
# setwd(workdir)
getwd()

# libraries
library(mvtnorm)
library(balance)
library(selbal)
library(microbenchmark)
library(ggplot2)
library(logratiolasso) # bates & tibshirani 2019

# Dr. Ma sources
source("RCode/func_libs.R")
source("COAT-master/coat.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "coat.R"))
source(paste0(functions_path, "principlebalances.R"))
source(paste0(functions_path, "selbal.R"))
source(paste0(functions_path, "constrainedlm.R"))

# helper functions
MSEyhat <- function(y, yhat){
  sqrt(sum((y - yhat)^2)) / length(y)
}
L2norm <- function(x){
  sqrt(crossprod(x))
}
L2distance <- function(x1, x2){
  sqrt(sum((x1 - x2)^2))
}

# settings
tol = 1e-4

# Cross-validation
cv.seed = 1234
cv.n_lambda = 100

# Bootstrap
bs.seed = 1997
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

# bs.finalfits = list()
bs.selected_variables = matrix(NA, num.genera, bs.n)
rownames(bs.selected_variables) = colnames(X)
for(b in 1:bs.n){
  set.seed(bs.seed + b)
  print(paste0("starting bootstrap ", b))
  
  # resample the data
  bs.resample = sample(1:n, n, replace = TRUE)
  log.X.prop.bs = log.X.prop[bs.resample, ]
  y.bs = y[bs.resample]
  
  # refitted CV
  # Split the data into 10 folds
  cv.K = 10
  set.seed(cv.seed)
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
    XYdata = data.frame(Xtrain, y = Ytrain)
    Lasso_j = ConstrLasso(
      Ytrain, Xtrain, Cmat = matrix(1, dim(Xtrain)[2], 1), nlam = cv.n_lambda,
      intercept = TRUE, scaling = TRUE, tol = tol)
    non0.betas = Lasso_j$bet != 0 # diff lambda = diff col
    for(m in 1:cv.n_lambda){
      # get refitted coefficients, after model selection and w/o penalization
      selected_variables = non0.betas[, m]
      if(all(!selected_variables)){ # if none selected
        refit = lm(Ytrain ~ 1)
        predictCLM = function(X) coefficients(refit)
      } else{ # otherwise, fit on selected variables
        # fit to the subsetted data
        Xtrain.sub = Xtrain[, selected_variables, drop = FALSE]
        Q = as.matrix(rep(1, sum(selected_variables)))
        # fit betabar
        stdXY = standardize(Xtrain.sub, Ytrain, center = TRUE, scale = TRUE)
        betabar = clm(stdXY$Xtilde, stdXY$Ytilde, Q)
        betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
        predictCLM = function(X){
          betabar.bstd$betahat0 +
            X[, rownames(betabar), drop = FALSE] %*% betabar.bstd$betahat
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
  Lasso_select = ConstrLasso(
    y.bs, log.X.prop.bs, Cmat = matrix(1, dim(log.X.prop.bs)[2], 1),
    lambda = lambda_min, nlam = 1,
    intercept = TRUE, scaling = TRUE, tol=tol)
  selected_variables = Lasso_select$bet != 0 # diff lambda = diff col
  # record which variables were selected in this fit
  bs.selected_variables[, b] = selected_variables
  # note: unless we're gonna estimate MSE(yhat), there is no need to get the
  #   actual constrained linear model.
  print(colnames(X)[selected_variables])
}
bs.selected_variables_numeric = apply(bs.selected_variables, 2, as.numeric)
bs.selection_percentages = apply(bs.selected_variables_numeric, 1, FUN = 
                                   function(x) sum(x, na.rm = TRUE))
names(bs.selection_percentages) = rownames(bs.selected_variables)
bs.results = list(
  seed = bs.seed,  
  selected_variables = bs.selected_variables, 
  selection_percentages = bs.selection_percentages
)
saveRDS(bs.results,
        file = paste0("Kristyn/ReproducingLin2014", 
                      "/bootstrap_simulations_121420",
                      "/bootstraps_cTsT.rds"))

sort(bs.selection_percentages)
