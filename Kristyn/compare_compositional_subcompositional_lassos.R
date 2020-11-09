workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
setwd(workdir)

# libraries
library(mvtnorm)
library(balance)
library(selbal)
library(propr)
library(microbenchmark)
library(ggplot2)
library(logratiolasso) # bates & tibshirani 2019
image_path = "/home/kristyn/Pictures"

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
source(paste0(functions_path, "propr.R"))
source(paste0(functions_path, "selbal.R"))

# settings

# Cross-validation
cv.seed = 1234
cv.n_lambda = 100

# Bootstrap
bs.seed = 1
bs.n = 100

# data
# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
DataFolder <- "/Data/"
load(paste0(workdir, DataFolder, "BMI.rda"))
# dim(raw_data) # 98 x 89
# dim(X) # 98 x 87
# dim(X.prop) # 98 x 87
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

################################################################################
# Lin et al. 2014 : Compositional Lasso
################################################################################

# Do cross-validation
# calculate squared error (prediction error?) for each fold, 
#   needed for CV(lambda) calculation
cvm = rep(NA, cv.n_lambda) # want to have CV(lambda)
cvm_sqerror = matrix(NA, cv.K, cv.n_lambda)
# Fit Lasso for each fold removed
for (j in 1:cv.K){
  # Training data
  Xtrain = X.prop[idfold != j, ]
  Ytrain = y[idfold != j]
  # Test data
  Xtest = X.prop[idfold == j, ]
  Ytest = y[idfold == j]
  
  # Fit LASSO on that fold using fitLASSOcompositional
  # first, take out columns that have all 0.5's, because they shouldn't be selected anyway (and lead to problems)
  cols.0.5 = apply(Xtrain, 2, FUN = function(vec) all(vec == 0.5))
  Xtrain = Xtrain[, !cols.0.5]
  Xtest = Xtest[, !cols.0.5]
  XYdata = data.frame(Xtrain, y = Ytrain)
  Lasso_j = fitCompositionalLASSO(Xtrain ,Ytrain, n_lambda = cv.n_lambda) # a problem in centering and scaling X cols with all 0.5's
  non0.betas = Lasso_j$beta_mat != 0 # diff lambda = diff col
  for(m in 1:cv.n_lambda){
    # get refitted coefficients, after model selection and w/o penalization
    selected_variables = non0.betas[, m]
    if(all(!selected_variables)){ # if none selected
      refit = lm(y ~ 1, data = XYdata)
    } else{ # otherwise, fit on selected variables
      refit = lm(
        as.formula(paste("y", "~",
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
cvm = (1 / n) * colSums(cvm_sqerror)

# Find lambda_min = argmin{CV(lambda)}
lambda_min_index = which.min(cvm)
lambda_min = Lasso_j$lambda_seq[lambda_min_index] # all Lasso_j generate the same lambda_seq
lambda_seq2 = Lasso_j$lambda_seq

# final fit
Lasso_select = fitCompositionalLASSO(X.prop, y, lambda_min)
XYdata = data.frame(X.prop, y = y)
non0.betas = Lasso_select$beta_mat != 0 # diff lambda = diff col
selected_variables = non0.betas
if(all(!selected_variables)){ # if none selected
  finalfit = lm(y ~ 1, data = XYdata)
} else{ # otherwise, fit on selected variables
  finalfit = lm(
    as.formula(paste("y", "~",
                     paste(colnames(XYdata)[which(selected_variables)], 
                           collapse = "+"),
                     sep = "")),
    data=XYdata)
}
finalfit


################################################################################
# Shi et al. 2016 : Sub-compositional Lasso
################################################################################

# Do cross-validation
# calculate squared error (prediction error?) for each fold, 
#   needed for CV(lambda) calculation
cvm2 = rep(NA, cv.n_lambda) # want to have CV(lambda)
cvm_sqerror2 = matrix(NA, cv.K, cv.n_lambda)
# Fit Lasso for each fold removed
for (j in 1:cv.K){
  # Training data
  Xtrain = X.prop[idfold != j, ]
  Ytrain = y[idfold != j]
  # Test data
  Xtest = X.prop[idfold == j, ]
  Ytest = y[idfold == j]
  
  # Fit LASSO on that fold using fitLASSOcompositional
  # first, take out columns that have all 0.5's, because they shouldn't be selected anyway (and lead to problems)
  cols.0.5 = apply(Xtrain, 2, FUN = function(vec) all(vec == 0.5))
  Xtrain = Xtrain[, !cols.0.5]
  Xtest = Xtest[, !cols.0.5]
  XYdata = data.frame(Xtrain, y = Ytrain)
  ##############################################################################
  Lasso_j = ConstrLasso(
    Ytrain, Xtrain, Cmat = rep(1, num.genera), lambda=NULL, nlam=20, intercept=TRUE, scaling=TRUE, maxiter=1000, tol=1e-8)
  # Lasso_j = fitCompositionalLASSO(Xtrain ,Ytrain, n_lambda = cv.n_lambda) # a problem in centering and scaling X cols with all 0.5's
  ##############################################################################
  non0.betas = Lasso_j$beta_mat != 0 # diff lambda = diff col
  for(m in 1:cv.n_lambda){
    # get refitted coefficients, after model selection and w/o penalization
    selected_variables = non0.betas[, m]
    if(all(!selected_variables)){ # if none selected
      refit = lm(y ~ 1, data = XYdata)
    } else{ # otherwise, fit on selected variables
      refit = lm(
        as.formula(paste("y", "~",
                         paste(colnames(XYdata)[which(selected_variables)], 
                               collapse = "+"),
                         sep = "")),
        data=XYdata)
    }
    # calculate squared error (prediction error?)
    # if(!all.equal(colnames(Xtest), colnames(Xtrain))) stop("Xtrain and Xtest don't have the same variable names")
    newx = Xtest[, selected_variables, drop = FALSE]
    Ypred = predict(refit, newdata = data.frame(newx), type = "response")
    cvm_sqerror2[j, m] = sum(crossprod(Ytest - Ypred))
  }
}

# Calculate CV(lambda) for each value of lambda
cvm2 = (1 / n) * colSums(cvm_sqerror2)

# Find lambda_min = argmin{CV(lambda)}
lambda_min_index2 = which.min(cvm2)
lambda_min2 = Lasso_j$lambda_seq[lambda_min_index2]

# final fit
################################################################################
Lasso_select2 #= fitCompositionalLASSO(X.prop, y, lambda_min)
################################################################################
XYdata = data.frame(X.prop, y = y)
non0.betas2 = Lasso_select2$beta_mat != 0 # diff lambda = diff col
selected_variables2 = non0.betas2
if(all(!selected_variables2)){ # if none selected
  finalfit2 = lm(y ~ 1, data = XYdata)
} else{ # otherwise, fit on selected variables
  finalfit2 = lm(
    as.formula(paste("y", "~",
                     paste(colnames(XYdata)[which(selected_variables2)], 
                           collapse = "+"),
                     sep = "")),
    data=XYdata)
}
finalfit2