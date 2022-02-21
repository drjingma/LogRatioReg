# uses limSolve::lsei() instead of my own clm()

# workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
# setwd(workdir)
getwd()

# ``To compare the prediction performance of the two methods, we randomly divided the data into a
# training set of 70 subjects and a test set of 28 subjects, and used the fitted model chosen by cross-
#   validation based on the training set to evaluate the prediction error on the test set. The prediction
# error averaged over 100 replicates was 30·30 for the proposed method and 30·55 for lasso (ii),
# with standard errors of 0·97 and 1·04, respectively, suggesting that the prediction performance
# of the proposed method is similar to or better than that of lasso (ii).''

# libraries
library(mvtnorm)
library(balance)
library(selbal)
library(microbenchmark)
library(ggplot2)
library(logratiolasso) # bates & tibshirani 2019
library(limSolve) # for constrained lm

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

# settings
std.center = TRUE
std.scale = FALSE
tol = 1e-4

# Cross-validation
cv.seed = 1234
cv.n_lambda = 100
cv.K = 10

# Repetitions
rep.seed = 1997
rep.n = 100
n.train = 70
n.test = 28

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

pred.err = matrix(NA, n.test, rep.n)
for(b in 1:rep.n){
  set.seed(rep.seed + b)
  print(paste0("starting bootstrap ", b))
  
  # split into train and test sets
  train.idx = sample(1:n, n.train)
  Xtrain = log.X.prop[train.idx, ]
  Ytrain = y[train.idx]
  Xtest = log.X.prop[-train.idx, ]
  Ytest = y[-train.idx]
  
  # refitted CV
  # Split the data into 10 folds
  
  # Fit Lasso on training set
  cv.fits = cv.func(
    method="ConstrLasso", y = Ytrain, x = Xtrain, 
    Cmat = matrix(1, dim(Xtrain)[2], 1), nlam = cv.n_lambda, 
    nfolds=5, tol = tol, seed = cv.seed + b)
  lambdamin.idx = which.min(cv.fits$cvm)
  intercept = cv.fits$int[lambdamin.idx]
  betabar = cv.fits$bet[, lambdamin.idx]
  names(betabar) = colnames(betabar)
  # refit?
  # selected_variables = (betahat != 0)
  # # get refitted coefficients, after model selection and w/o penalization
  # if(all(!selected_variables)){ # if none selected
  #   refit = lm(Ytrain ~ 1)
  #   predictCLM = function(X) coefficients(refit)
  # } else{ # otherwise, fit on selected variables
  #   # fit to the subsetted data
  #   Xtrain.sub = Xtrain[, selected_variables, drop = FALSE]
  #   Xtrain.sub2 = cbind(1, Xtrain.sub)
  #   Q = as.matrix(rep(1, sum(selected_variables)))
  #   Q2 = rbind(0, Q)
  #   colnames(Xtrain.sub2)[1] = "Intercept"
  #   rownames(Q2) = colnames(Xtrain.sub2)
  #   lsei.fit = lsei(A = Xtrain.sub2, B = Ytrain, E = t(Q2), F = 0)
  #   predictCLM = function(X){
  #     lsei.fit$X[1] +
  #       X[, selected_variables, drop = FALSE] %*% lsei.fit$X[-1]
  #   }
  predictCLM = function(X){
    intercept + X %*% betabar
  }
  # calculate squared error (prediction error?)
  Ypred = predictCLM(Xtest)
  pred.err[, b] = (Ytest - Ypred)^2
}

dim(pred.err)

mean.pred.err = colMeans(pred.err)
mean(mean.pred.err)

saveRDS(pred.err,
        file = paste0("Kristyn/ReproducingLin2014", 
                      "/bootstrap_simulations_121420",
                      "/prediction_performance_notparallel.rds"))
