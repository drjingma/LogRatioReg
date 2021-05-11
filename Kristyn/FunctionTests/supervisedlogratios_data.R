getwd()

# libraries
library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
library(limSolve) # for constrained lm, lsei()
library(stats) # for hclust()
library(balance) # for sbp.fromHclust()

# Dr. Ma sources
library(Matrix)
library(glmnet)
library(compositions)
library(stats)
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))

# Method Settings
nlam = 200
cv.n_lambda = 200
cv.K = 10

# Data Split
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
# Experiments
################################################################################

# set.seed(20)
# split into train and test sets
train.idx = sample(1:n, n.train)
Xtrain = X.prop[train.idx, ]
Ytrain = y[train.idx]
Xtest = X.prop[-train.idx, ]
Ytest = y[-train.idx]

# refitted CV
# Split the data into 10 folds
seed = 1
# Fit Lasso on training set
set.seed(seed)
slr = cvSLR.old(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K)
lambdamin.idx = which.min(slr$cvm)
# par(mfrow = c(2, 1))
plot(slr$lambda, slr$cvm, xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
     col = rgb(0,0,0, alpha = 0.25))
points(slr$lambda[lambdamin.idx], min(slr$cvm), col = 2, pch = 15)
text(slr$lambda[lambdamin.idx], min(slr$cvm), 
     paste0("(", round(slr$lambda[lambdamin.idx], 2), ", ", 
            round(min(slr$cvm), 2), ")"), pos = 3)
set.seed(seed)
slr2 = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K, 
             foldid = slr$foldid)
lambdamin.idx2 = which.min(slr2$cvm)
plot(slr2$lambda, slr2$cvm, xlim = c(range(slr2$lambda) + c(-0.1, 0.2)), 
     col = rgb(0,0,0, alpha = 0.25))
points(slr2$lambda[lambdamin.idx2], min(slr2$cvm), col = 2, pch = 15)
text(slr2$lambda[lambdamin.idx2], min(slr2$cvm), 
     paste0("(", round(slr2$lambda[lambdamin.idx2], 2), ", ", 
            round(min(slr2$cvm), 2), ")"), pos = 3)
all.equal(slr$lambda, slr2$lambda)
# get fitted model
a0 = slr$int[lambdamin.idx]
a0
# > a0
# s17 
# 24.28477  
betahat = slr$bet[, lambdamin.idx]
as.numeric(betahat)
# > as.numeric(betahat)
# [1]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [7]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [13]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [19]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [25]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [31]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [37]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [43]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [49]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [55]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [61]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [67]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [73]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [79]  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
# [85] -0.05039666 -1.09195961
btree = slr$btree
predictSLR = function(X){
  a0 + computeBalances(X, btree) %*% betahat
}
# calculate squared error (prediction error?)
Ypred = predictSLR(Xtest)
as.numeric(Ypred)
# > as.numeric(Ypred)
# [1] 24.24984 24.24984 24.31970 24.51785 26.33160 25.79855 24.28477 24.24984 24.40562
# [10] 24.28477 24.51584 24.13911 24.17997 24.37292 25.79855 24.21491 24.28477 24.28477
# [19] 24.28477 24.29396 24.14504 24.21491 24.28477 24.61150 22.21220 23.37691 24.28477
# [28] 24.28477
PE = crossprod(Ypred - Ytest) / n
PE
# > PE # seed = 15
# [,1]
# [1,] 12.17391
# > PE # seed = 1
# [,1]
# [1,] 7.979851
# > PE # seed = 20
# [,1]
# [1,] 13.70292

# look at the transformation
betahat.tr = LRtoLC(betahat, btree)
a0.tr = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat.tr)
predictLLC = function(Z){
  a0.tr + Z %*% betahat.tr
}
Ypred.tr = predictLLC(log(Xtest))
as.numeric(Ypred.tr)
# > as.numeric(Ypred.tr)
# [1] 24.24984 24.24984 24.31970 24.51785 26.33160 25.79855 24.28477 24.24984 24.40562
# [10] 24.28477 24.51584 24.13911 24.17997 24.37292 25.79855 24.21491 24.28477 24.28477
# [19] 24.28477 24.29396 24.14504 24.21491 24.28477 24.61150 22.21220 23.37691 24.28477
# [28] 24.28477
all.equal(Ypred, Ypred.tr)
# > all.equal(Ypred, Ypred.tr)
# [1] TRUE



################################################################################
# adding intercept argument ####################################################
# when intercept = TRUE, results should be the same
# Fit Lasso on training set
slr = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K, intercept = TRUE)
lambdamin.idx = which.min(slr$cvm)
# get fitted model
a0 = slr$int[lambdamin.idx]
# a0 = 
betahat = slr$bet[, lambdamin.idx]
# betahat = 
btree = slr$btree
predictSLR = function(X){
  a0 + computeBalances(X, btree) %*% betahat
}
# calculate squared error (prediction error?)
Ypred = predictSLR(Xtest)
# Ypred = 
PE = crossprod(Ypred - Ytest) / n
PE
# PE = 

# look at the transformation
betahat.tr = LRtoLC(betahat, btree)
a0.tr = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat.tr)
predictLLC = function(Z){
  a0.tr + Z %*% betahat.tr
}
Ypred.tr = predictLLC(log(Xtest))
# Ypred.tr = 
all.equal(Ypred, Ypred.tr)





################################################################################
# when intercept = FALSE
