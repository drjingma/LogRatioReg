rm(list=ls())
# Purpose: observe different R^2 values for specified rho and sigma_eps values
#   for the one-balance model
#   this time, more generalizeable to different balances
# Date: 2/15/2022

################################################################################
# libraries and settings

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(propr)

source("RCode/func_libs.R")
source("RCode/slr-main.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/supervisedlogratioseta.R")
source("Kristyn/Functions/HSClust.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/simulatedata.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

# expdecay Sigma #############################################################

# Settings to toggle with
sigma.settings = "expdecaySigma"
theta.value = 1
intercept = TRUE
K = 10
n = 100
p = 30
scaling = TRUE
linkage = "average"
tol = 1e-4
nlam = 100
neta = p
#################
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
#################
rho = 0 #
desired_Rsquared = 1 #
sigma_eps = get_sigma_eps(
  sbp = SBP.true, ilr.trans.constant = ilrtrans.true$const, theta = theta.value, 
  Rsq = desired_Rsquared, rho = rho)
#################

# Population parameters
SigmaW = rgExpDecay(p, rho)$Sigma
muW = c(rep(log(p), 5), rep(0, p - 5))
names(muW) = paste0('s', 1:p)

##############################################################################
# generate data
# set.seed(1947)
set.seed(1234)
# generate X
logW.all <- mvrnorm(n = 2 * n, mu = muW, Sigma = SigmaW)
W.all <- exp(logW.all)
X.all <- sweep(W.all, 1, rowSums(W.all), FUN='/')
colnames(X.all) = paste0('s', 1:p)
# create the ilr(X.all) covariate by hand to
#   generate y
beta = as.numeric(ilrtrans.true$ilr.trans * theta.value)
y.all = as.numeric(log(X.all) %*% beta) + rnorm(n) * sigma_eps

# subset out training and test sets
X = X.all[1:n, ]
X.test = X.all[-(1:n), ]
Y <- y.all[1:n]
Y.test <- y.all[-(1:n)]

# about beta
names(beta) <- paste0('s', 1:p)
non0.beta = (beta != 0)
is0.beta = abs(beta) <= 10e-8
bspars = sum(non0.beta)

##############################################################################
# estimate Rsquared, given the true model
SSres = sum((y.all - as.numeric(log(X.all) %*% beta))^2)
SStot = sum((y.all - mean(y.all))^2)
Rsq = 1 - SSres/SStot

##############################################################################
# Dr. Ma's new slr
##############################################################################

# apply supervised log-ratios, using CV to select threshold
slrnew = slr(x = X, y = Y)
slrnew_activevars = names(slrnew$index)
slrnew_SBP = matrix(slrnew$index)
rownames(slrnew_SBP) = slrnew_activevars

slrnew.coefs = coefficients(slrnew$model)
slrnew.a0 = slrnew.coefs[1]
slrnew.thetahat = slrnew.coefs[-1]

slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew.betahat) = names(beta)
slrnew.betahat[slrnew_activevars, ] =
  as.numeric(slrnew.betahat.nonzero)

# compute metrics on the selected model #
slrnew.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test,
  ilrX.train = getIlrX(X[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
  ilrX.test = getIlrX(X.test[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
  n.train = n, n.test = n,
  thetahat0 = slrnew.a0, thetahat = slrnew.thetahat, betahat = slrnew.betahat,
  true.sbp = SBP.true,
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

##############################################################################
# compositional lasso (a linear log contrast method)
##############################################################################
classo = cv.func(
  method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
  nfolds = K, tol = tol, intercept = intercept, scaling = scaling)

cl.lam.min.idx = which.min(classo$cvm)
cl.a0 = classo$int[cl.lam.min.idx]
cl.betahat = classo$bet[, cl.lam.min.idx]

# compute metrics on the selected model #
cl.metrics = getMetricsLLC(
  y.train = Y, y.test = Y.test,
  logX.train = log(X),
  logX.test = log(X.test),
  n.train = n, n.test = n,
  betahat0 = cl.a0, betahat = cl.betahat,
  true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta, 
  metrics = c("prediction", "selection"))





