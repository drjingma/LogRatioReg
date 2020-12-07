# workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
# setwd(workdir)

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
source(paste0(functions_path, "classic_lasso.R"))
source(paste0(functions_path, "compositional_lasso.R"))
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "coat.R"))
source(paste0(functions_path, "principlebalances.R"))
source(paste0(functions_path, "selbal.R"))

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
# Choosing the tuning parameter -- an example of what goes on in the bootstrap #
################################################################################

# Lin et. al. 2014 applied the proposed method using a refitted version of 
#   ten-fold cross-validation to choose the tuning parameter, where
#   the prediction of each sample split is computed with the refitted coeffs
#   (obtained after model selection and without penalization)

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
  Xtrain = log.X.prop[idfold != j, ]
  Ytrain = y[idfold != j]
  # Test data
  Xtest = log.X.prop[idfold == j, ]
  Ytest = y[idfold == j]
  
  # Fit LASSO on that fold using fitLASSOcompositional
  # first, take out columns that have all 0.5's, because they shouldn't be selected anyway (and lead to problems)
  # cols.0.5 = apply(Xtrain, 2, FUN = function(vec) all(vec == 0.5))
  # Xtrain = Xtrain[, !cols.0.5]
  # Xtest = Xtest[, !cols.0.5]
  XYdata = data.frame(Xtrain, y = Ytrain)
  Lasso_j = ConstrLasso(
    Ytrain, Xtrain, Cmat = matrix(1, dim(Xtrain)[2], 1), nlam = cv.n_lambda, 
    intercept = TRUE, scaling = TRUE, tol = tol)
  non0.betas = Lasso_j$bet != 0 # diff lambda = diff col
  for(m in 1:cv.n_lambda){
    # get refitted coefficients, after model selection and w/o penalization
    selected_variables = non0.betas[, m]
    if(all(!selected_variables)){ # if none selected
      refit = lm(y ~ 1, data = XYdata)
    } else{ # otherwise, fit on selected variables
      refit = lm(
        as.formula(paste("y", "~-1+",
                         paste(colnames(XYdata)[which(selected_variables)], 
                               collapse = "+"),
                         sep = "")),
        data=XYdata)
    }
    # calculate squared error (prediction error?)
    # if(!all.equal(colnames(Xtest), colnames(Xtrain))) stop("Xtrain and Xtest don't have the same variable names")
    newx = Xtest[, selected_variables, drop = FALSE]
    Ypred = predict(refit, newdata = data.frame(newx), type = "response")
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
  y, log.X.prop, Cmat = matrix(1, dim(log.X.prop)[2], 1), 
  lambda = lambda_min, nlam = 1, 
  intercept=TRUE, scaling=TRUE, tol=tol)
XYdata = data.frame(log.X.prop, y = y)
non0.betas = Lasso_select$bet != 0 # diff lambda = diff col
selected_variables = non0.betas
if(all(!selected_variables)){ # if none selected
  finalfit = lm(y ~ 1, data = XYdata)
} else{ # otherwise, fit on selected variables
  finalfit = lm(
    as.formula(paste("y", "~-1+",
                     paste(colnames(XYdata)[which(selected_variables)], 
                           collapse = "+"),
                     sep = "")),
    data=XYdata)
}
finalfit
# ... now what?
# do this a bunch of times (bootstrap)

################################################################################
# Actual work: Selecting the genera #
################################################################################

# They generate 100 bootstrap samples and use the same CV procedure to select 
#   the genera (for stable selection results)

bs.finalfits = list()
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
    # first, take out columns that have all 0.5's, because they shouldn't be selected anyway (and lead to problems)
    # cols.0.5 = apply(Xtrain, 2, FUN = function(vec) all(vec == 0.5))
    # Xtrain = Xtrain[, !cols.0.5]
    # Xtest = Xtest[, !cols.0.5]
    XYdata = data.frame(Xtrain, y = Ytrain)
    Lasso_j = ConstrLasso(
      Ytrain, Xtrain, Cmat = matrix(1, dim(Xtrain)[2], 1), nlam = cv.n_lambda, 
      intercept = TRUE, scaling = TRUE, tol = tol)
    non0.betas = Lasso_j$bet != 0 # diff lambda = diff col
    for(m in 1:cv.n_lambda){
      # get refitted coefficients, after model selection and w/o penalization
      selected_variables = non0.betas[, m]
      if(all(!selected_variables)){ # if none selected
        refit = lm(y ~ 1, data = XYdata)
      } else{ # otherwise, fit on selected variables
        refit = lm(
          as.formula(paste("y", "~-1+",
                           paste(colnames(XYdata)[which(selected_variables)], 
                                 collapse = "+"),
                           sep = "")),
          data=XYdata)
      }
      # calculate squared error (prediction error?)
      # if(!all.equal(colnames(Xtest), colnames(Xtrain))) stop("Xtrain and Xtest don't have the same variable names")
      newx = Xtest[, selected_variables, drop = FALSE]
      Ypred = predict(refit, newdata = data.frame(newx), type = "response")
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
    lambda = lambda_min, nlam = 1, 
    intercept=TRUE, scaling=TRUE, tol=tol)
  XYdata = data.frame(log.X.prop.bs, y = y.bs)
  non0.betas = Lasso_select.bs$bet != 0 # diff lambda = diff col
  selected_variables = non0.betas
  if(all(!selected_variables)){ # if none selected
    finalfit.bs = lm(y ~ 1, data = XYdata)
  } else{ # otherwise, fit on selected variables
    finalfit.bs = lm(
      as.formula(paste("y", "~-1+",
                       paste(colnames(XYdata)[which(selected_variables)], 
                             collapse = "+"),
                       sep = "")),
      data=XYdata)
  }
  
  # save the final fit for this bootstrap
  bs.finalfits[[b]] = finalfit.bs
  
  # record which variables were selected in this fit
  bs.selected_variables[, b] = selected_variables
}
bs.selected_variables_numeric = apply(bs.selected_variables, 2, as.numeric)
bs.selection_percentages = apply(bs.selected_variables_numeric, 1, sum)
names(bs.selection_percentages) = rownames(bs.selected_variables)
bs.results = list(
  seed = bs.seed, 
  finalfits = bs.finalfits, 
  selected_variables = bs.selected_variables, 
  selection_percentages = bs.selection_percentages
)
saveRDS(bs.results, 
        file = paste0("Kristyn/ReproducingLin2014/lin2014_bootstrap_results_lmnoint.rds"))

which(bs.selection_percentages > 50)
sort(bs.selection_percentages)
bs.selection_percentages["Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella"]

