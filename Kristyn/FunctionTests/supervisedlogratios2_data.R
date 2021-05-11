getwd()

# libraries
library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
library(limSolve) # for constrained lm, lsei()
library(stats) # for hclust()
library(balance) # for sbp.fromHclust()
library(rare)

# Dr. Ma sources
library(Matrix)
library(glmnet)
library(compositions)
library(stats)
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "supervisedlogratios2.R"))
sourceCpp(paste0(functions_path, "RcppExports.cpp"))
sourceCpp(paste0(functions_path, "RcppExports.cpp"))

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
# seed = 1
# Fit Lasso on training set
# set.seed(seed)
# slr0 = cvSLR.old(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K)
# lambdamin.idx = which.min(slr0$cvm)
# # par(mfrow = c(2, 1))
# plot(slr0$lambda, slr0$cvm, xlim = c(range(slr0$lambda) + c(-0.1, 0.2)), 
#      col = rgb(0,0,0, alpha = 0.25))
# points(slr0$lambda[lambdamin.idx], min(slr0$cvm), col = 2, pch = 15)
# text(slr0$lambda[lambdamin.idx], min(slr0$cvm), 
#      paste0("(", round(slr0$lambda[lambdamin.idx], 2), ", ", 
#             round(min(slr0$cvm), 2), ")"), pos = 3)
set.seed(1)
slr = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K)
lambdamin.idx = which.min(slr$cvm)
lambdamin = slr$lambda[lambdamin.idx]
plot(slr$lambda, slr$cvm, #xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
     col = rgb(0,0,0, alpha = 0.25))
# points(slr$lambda[lambdamin.idx], min(slr$cvm), col = 2, pch = 15)
# text(slr$lambda[lambdamin.idx], min(slr$cvm), 
#      paste0("(", round(slr$lambda[lambdamin.idx], 2), ", ", 
#             round(min(slr$cvm), 2), ")"), pos = 3)


set.seed(1)
btree = getSupervisedTree(Ytrain, Xtrain)
U = getU(btree)
A = Matrix(U, sparse = TRUE)
# slr2 = cvSLR2(y = Ytrain, X = Xtrain, A = A, Q = rare:::svdA(A), 
#               nlam = cv.n_lambda, nfolds = cv.K)
slr2.1 = cvSLR2(y = Ytrain, X = Xtrain, A = A, 
                nlam = cv.n_lambda, nfolds = cv.K, alpha = 1)
plot(slr2.1$lambda, slr2.1$cvm, col = rgb(0,0,0, alpha = 0.25))

lambdamin.idx2 = which.min(slr2.1$cvm)
lambdamin2 = slr2.1$lambda[lambdamin.idx2]
