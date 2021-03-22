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

################################################################################
# Using un-normalized SBP matrix, V, as transformation matrix, instead of U
################################################################################
# U.type = "sbp"
# 
# set.seed(1)
# slr = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K, 
#             U.type = U.type)
# lambdamin.idx = which.min(slr$cvm)
# plot(slr$lambda, slr$cvm, xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
#      col = rgb(0,0,0, alpha = 0.25))
# points(slr$lambda[lambdamin.idx], min(slr$cvm), col = 2, pch = 15)
# text(slr$lambda[lambdamin.idx], min(slr$cvm), 
#      paste0("(", round(slr$lambda[lambdamin.idx], 2), ", ", 
#             round(min(slr$cvm), 2), ")"), pos = 3)
# # get fitted model
# a0 = slr$int[lambdamin.idx]
# thetahat = slr$bet[, lambdamin.idx]
# btree = slr$btree
# 
# # calculate predicted response -- version 1: ilr(X) = log(X) U
# ilrX1 = unname(computeBalances(Xtest, btree, U.type = "sbp"))
# Ypred1.ilr = a0 + ilrX1 %*% thetahat
# # look at the transformation from theta to beta
# U1 = getU(btree, U.type = "sbp") # same as sbp.fromHclust(btree)
# betahat1 = U1 %*% thetahat
# betahat01 = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat1)
# a0.calc = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% thetahat)
# all.equal(unname(a0), unname(betahat01)) # the intercept doesn't change!!!
# # look at the transformed calculation for Yhat1
# Ypred1.llc = betahat01 + log(Xtest) %*% betahat1
# all.equal(Ypred1.ilr, Ypred1.llc) # they are equal!
# # calculated prediction error
# PE1 = crossprod(Ypred1.ilr - Ytest) / n
# 
# # calculate predicted response -- version 2: ilr(X) using balance.fromSBP(X, U)
# ilrX2 = unname(computeBalances(Xtest, btree, U.type = "balance.fromSBP"))
# Ypred2.ilr = a0 + ilrX2 %*% thetahat
# PE2 = crossprod(Ypred2.ilr - Ytest) / n
# # look at the transformation from theta to beta
# U2 = getU(btree, U.type = "balance.fromSBP")
# betahat2 = U2 %*% thetahat
# betahat02 = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat2) # BUT IT WORKS WHEN I USE betahat1
# # look at the transformed calculation for Yhat2
# Ypred2.llc = betahat02 + log(Xtest) %*% betaha2
# all.equal(Ypred2.ilr, Ypred2.llc) # not exactly the same...... why not?? # EQUAL WHEN I USE betahat1 in betahat01 calculation -- could be due to the slr calculation
# all.equal(ilrX2 %*% thetahat, log(Xtest) %*% betahat2)
# 
# # summary(Ypred2.ilr - Ypred2.llc)
# a0 - betahat02
# betahat01 - betahat02
# # data.frame(ilr = Ypred2.ilr, llc = Ypred2.llc, diff = Ypred2.ilr - Ypred2.llc)
# # summary(Ypred1.ilr - Ypred2.ilr)
# # hist(ilrX1 - ilrX2, 10)
# # data.frame(ilr1 = Ypred1.ilr, ilr2 = Ypred2.ilr, diff12 = Ypred1.ilr - Ypred2.ilr)
# data.frame(ilr1 = c(betahat01, betahat1),
#            ilr2 = c(betahat02, betahat2),
#            diff12 = c(betahat01, betahat1) - c(betahat02, betahat2))
# 
# # calculate ilrX3 -- version 3: orthonormalization of SBP by hand
# ilrX3 = unname(computeBalances(Xtest, btree, U.type = "orthonormal"))
# # compare ilrX calculations using SBP (contrasts), balance.fromSBP, and by hand!
# all.equal(ilrX1, ilrX2)
# all.equal(ilrX1, ilrX3)
# all.equal(ilrX2, ilrX3) # EQUAL! (balance.fromSBP and by hand)
# U4 = getU(btree, U.type = 3)
# ilrX4 = unname(log(Xtest) %*% U4)
# all.equal(ilrX3, ilrX4) # computing by hand works.


################################################################################
# Using normalized SBP matrix, U, as transformation matrix
################################################################################
U.type = "balance.fromSBP"

set.seed(1)
slr = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K, 
            U.type = U.type)
lambdamin.idx = which.min(slr$cvm)
plot(slr$lambda, slr$cvm, xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
     col = rgb(0,0,0, alpha = 0.25))
points(slr$lambda[lambdamin.idx], min(slr$cvm), col = 2, pch = 15)
text(slr$lambda[lambdamin.idx], min(slr$cvm), 
     paste0("(", round(slr$lambda[lambdamin.idx], 2), ", ", 
            round(min(slr$cvm), 2), ")"), pos = 3)
# get fitted model
thetahat0.U = slr$int[lambdamin.idx]
thetahat.U = slr$bet[, lambdamin.idx]
btree = slr$btree # doesnt depend on transformation matrix

# calculate predicted response -- version 1: ilr(X) = log(X) U
ilrX.V = unname(computeBalances(Xtest, btree, U.type = "sbp"))
Ypred.ilr.V = thetahat0.U + ilrX.V %*% thetahat.U
# look at the transformation from theta to beta
V = getU(btree, U.type = "sbp") # same as sbp.fromHclust(btree)
betahat.V = V %*% thetahat.U
betahat0.V = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat.V)
# look at the transformed calculation for Yhat1
Ypred.llc.V = betahat0.V + log(Xtest) %*% betahat.V
all.equal(Ypred.ilr.V, Ypred.llc.V) # they are NOT equal here!! possibly bc using wrong sbp type 
data.frame(Ypred.ilr.V = Ypred.ilr.V, Ypred.llc.V = Ypred.llc.V, diff = Ypred.ilr.V - Ypred.llc.V)
all.equal(ilrX.V %*% thetahat.U, log(Xtest) %*% betahat.V) # it's at the coefficients that they differ

# calculate predicted response -- version 2: ilr(X) using balance.fromSBP(X, U)
ilrX.U = unname(computeBalances(Xtest, btree, U.type = "balance.fromSBP"))
Ypred.ilr.U = thetahat0.U + ilrX.U %*% thetahat.U
PE2 = crossprod(Ypred.ilr.U - Ytest) / n
# look at the transformation from theta to beta
U = getU(btree, U.type = "balance.fromSBP")
betahat.U = U %*% thetahat.U
betahat0.U = as.numeric(mean(Ytrain) - colMeans(log(Xtrain)) %*% betahat.U) # BUT IT WORKS WHEN I USE betahat1
# look at the transformed calculation for Yhat2
Ypred.llc.U = betahat0.U + log(Xtest) %*% betahat.U
all.equal(Ypred.ilr.U, Ypred.llc.U) # equal here!!! matches U.type!!!!!!!!!!!!!!!!!
# data.frame(ilr.V = c(betahat0.V, betahat.V), 
#            ilr.U = c(betahat0.U, betahat.U), 
#            diff.VU = c(betahat0.V, betahat.V) - c(betahat0.U, betahat.U))
# calculate prediction error
PE2 = crossprod(Ypred2.ilr - Ytest) / n
PE2.2 = crossprod(Ypred2.llc - Ytest) / n
a0 - betahat01
betahat02 - betahat01
all.equal(unname(a0), unname(betahat02))

# calculate ilrX3 -- version 3: orthonormalization of SBP by hand
ilrX3 = unname(computeBalances(Xtest, btree, U.type = "orthonormal"))
# compare ilrX calculations using SBP (contrasts), balance.fromSBP, and by hand!
all.equal(ilrX1, ilrX2)
all.equal(ilrX1, ilrX3)
all.equal(ilrX2, ilrX3) # EQUAL! (balance.fromSBP and by hand)
U4 = getU(btree, U.type = 3)
ilrX4 = unname(log(Xtest) %*% U4)
all.equal(ilrX3, ilrX4) # computing by hand works.



