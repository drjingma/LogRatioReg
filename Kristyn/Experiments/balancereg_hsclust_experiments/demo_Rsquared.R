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
U.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
#################
get_sigma_eps = function(sbp, ilr.trans.constant, theta, Rsq, rho){
  # calculate Var[ilr(X)]
  pos.idx = which(sbp == 1)
  num.pos = length(pos.idx)
  neg.idx = which(sbp == -1)
  num.neg = length(neg.idx)
  pairs.pos <- pairs.neg <- c()
  if(num.pos > 1){
    pairs.pos = t(combn(pos.idx, 2))
  }
  if(num.neg > 1){
    pairs.neg = t(combn(neg.idx, 2))
  }
  Sigma.ij.pos = rho^abs(pairs.pos[, 1] - pairs.pos[, 2])
  Sigma.ij.neg = rho^abs(pairs.neg[, 1] - pairs.neg[, 2])
  varIlrX = ilr.trans.constant^2 * (
    1 / num.pos + 1 / num.neg + 
      2 * sum(Sigma.ij.pos) / num.pos^2 + 2 * sum(Sigma.ij.neg) / num.neg^2
    )
  # get sigma.eps.sq
  sigma.eps.sq = theta^2 * varIlrX * (1 / Rsq - 1)
  # return sigma.eps
  return(sqrt(sigma.eps.sq))
}
rho = 0.2 #
desired_Rsquared = 0.6 #
sigma_eps = get_sigma_eps(
  sbp = SBP.true, ilr.trans.constant = U.true$const, theta = theta.value, 
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
beta = as.numeric(U.true$ilr.trans * theta.value)
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

# ##############################################################################
# # supervised log-ratios (a balance regression method) with eta
# #   -- hierarchical spectral clustering
# ##############################################################################
# # apply hierarchical spectral clustering to the SLR similarity matrix
# slrSimMat = getSlrMatrix(
#   y = Y, X = X, type = "similarity")
# 
# # kmeans
# slrhsc_btree = HSClust(
#   W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
# slrhsc_SBP = sbp.fromHSClust(
#   levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
# # plotSBP(slrhsc_SBP)
# fields::image.plot(slrSimMat)
# 
# # # shi-malik
# # slrhsc_btree2 = HSClust(
# #   W = slrSimMat, force_levelMax = TRUE, method = "shimalik")
# # slrhsc_SBP2 = sbp.fromHSClust(
# #   levels_matrix = slrhsc_btree2$allLevels, row_names = names(beta))
# # plotSBP(slrhsc_SBP2)
# 
# # apply supervised log-ratios, using CV to select threshold
# slrhsc2 = cvBMLassoThresh(
#   y = Y, X = X,
#   W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
#   hsc_method = "kmeans", # "shimalik", "kmeans"
#   # thresh_method = "original",
#   # multiple_balances = TRUE,
#   force_levelMax = TRUE,
#   stopping_rule = NULL,
#   sbp = slrhsc_SBP,
#   neta = neta, nlam = nlam,
#   nfolds = K, foldid = NULL,
#   intercept = intercept,
#   standardize = scaling,
#   prints = TRUE,
#   seed = 123
# )
# 
# slrhsc2.eta.min.idx = slrhsc2$min.idx[2]
# slrhsc2.lam.min.idx = slrhsc2$min.idx[1]
# slrhsc2.a0 = slrhsc2$theta0[[slrhsc2.eta.min.idx]][slrhsc2.lam.min.idx]
# slrhsc2.thetahat = slrhsc2$theta[[slrhsc2.eta.min.idx]][, slrhsc2.lam.min.idx]
# slrhsc2.SBP = slrhsc2$sbp_thresh[[slrhsc2.eta.min.idx]]
# slrhsc2.betahat.nonzero = getBetaFromTheta(slrhsc2.thetahat, sbp = slrhsc2.SBP)
# slrhsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
# rownames(slrhsc2.betahat) = names(beta)
# slrhsc2.betahat[slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], ] =
#   as.numeric(slrhsc2.betahat.nonzero)
# slrhsc2.betahat
# 
# # compute metrics on the selected model #
# slrhsc2.metrics = getMetricsBalanceReg(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(
#     X[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
#     sbp = slrhsc2.SBP),
#   ilrX.test = getIlrX(
#     X.test[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
#     sbp = slrhsc2.SBP),
#   n.train = n, n.test = n,
#   thetahat0 = slrhsc2.a0, thetahat = slrhsc2.thetahat,
#   betahat = slrhsc2.betahat,
#   sbp = slrhsc2.SBP,
#   true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
# 
# # plot the tree given by slr-hsc, indicating significant covariates
# slrhsc2_leaf_types = rep("covariate", nrow(slrhsc2.SBP))
# slrhsc2_balance_types = rep("balance", ncol(slrhsc2.SBP))
# slrhsc2_nodes_types = data.frame(
#   name = c(colnames(slrhsc2.SBP), rownames(slrhsc2.SBP)),
#   type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
# )
# plotSBP(slrhsc2.SBP, title = "slr-hsc-eta (lasso)", nodes_types = slrhsc2_nodes_types)
# 
# # note: centering X and y didn't change the slrSimMat
# # note 2: using clrXj - clrXk instead of logXj - logXk didn't change slrSimMat
# 
# # heat map of cross-validated mse's
# # fields::image.plot(slrhsc2$cvm)

##############################################################################
# Dr. Ma's new slr
##############################################################################

# apply supervised log-ratios, using CV to select threshold
slr.main = slr(
  x = X, y = Y
)

slrhsc2.eta.min.idx = slrhsc2$min.idx[2]
slrhsc2.lam.min.idx = slrhsc2$min.idx[1]
slrhsc2.a0 = slrhsc2$theta0[[slrhsc2.eta.min.idx]][slrhsc2.lam.min.idx]
slrhsc2.thetahat = slrhsc2$theta[[slrhsc2.eta.min.idx]][, slrhsc2.lam.min.idx]
slrhsc2.SBP = slrhsc2$sbp_thresh[[slrhsc2.eta.min.idx]]
slrhsc2.betahat.nonzero = getBetaFromTheta(slrhsc2.thetahat, sbp = slrhsc2.SBP)
slrhsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrhsc2.betahat) = names(beta)
slrhsc2.betahat[slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], ] =
  as.numeric(slrhsc2.betahat.nonzero)
slrhsc2.betahat

# compute metrics on the selected model #
slrhsc2.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test,
  ilrX.train = getIlrX(
    X[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
    sbp = slrhsc2.SBP),
  ilrX.test = getIlrX(
    X.test[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
    sbp = slrhsc2.SBP),
  n.train = n, n.test = n,
  thetahat0 = slrhsc2.a0, thetahat = slrhsc2.thetahat,
  betahat = slrhsc2.betahat,
  sbp = slrhsc2.SBP,
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)