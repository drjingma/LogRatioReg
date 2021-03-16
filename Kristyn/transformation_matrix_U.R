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
cv.n_lambda = 50
cv.K = 3

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

set.seed(20)
# split into train and test sets
train.idx = sample(1:n, n.train)
Xtrain = X.prop[train.idx, ]
Ytrain = y[train.idx]
Xtest = X.prop[-train.idx, ]
Ytest = y[-train.idx]

set.seed(1)
slr = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K)
lambdamin.idx = which.min(slr$cvm)
plot(slr$lambda, slr$cvm, xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
     col = rgb(0,0,0, alpha = 0.25))
points(slr$lambda[lambdamin.idx2], min(slr$cvm), col = 2, pch = 15)
text(slr$lambda[lambdamin.idx2], min(slr$cvm), 
     paste0("(", round(slr$lambda[lambdamin.idx2], 2), ", ", 
            round(min(slr$cvm), 2), ")"), pos = 3)
# get fitted model
a0 = slr$int[lambdamin.idx]
thetahat = slr$bet[, lambdamin.idx]
btree = slr$btree

# calculate predicted response -- version 1: ilr(X) = log(X) U
U = sbp.fromHclust(btree)
# unname(U) # what U looks like
ilrX1 = log(Xtest) %*% U
Ypred1 = a0 + ilrX1 %*% thetahat
PE1 = crossprod(Ypred1 - Ytest) / n

# look at the transformation
betahat = U %*% thetahat
a0.tr = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat)
Ypred.tr = a0.tr + log(Xtest) %*% betahat
all.equal(Ypred1, Ypred.tr)


# calculate predicted response -- version 2: ilr(X) using balance.fromSBP(X, U)
ilrX2 = balance.fromSBP(Xtest, U)
Ypred2 = a0 + ilrX2 %*% thetahat
PE2 = crossprod(Ypred2 - Ytest) / n


