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
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/supervisedlogratioseta.R")
source("Kristyn/Functions/HSClust.R")
source("Kristyn/Functions/slrnew.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/simulatedata.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

# expdecay Sigma #############################################################

# Settings to toggle with
sigma.settings = "latentVarModel"
intercept = TRUE
K = 10
n = 100
p = 30
scaling = TRUE
tol = 1e-4
nlam = 100
neta = p
rho = 0.2 
sigma_eps1 = 0
sigma_eps2 = 0
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)

##############################################################################
# generate data
b0 = 0
b1 = 1
a0 = 0
a1 = 1
theta.value = 1
# set.seed(1947)
set.seed(1234)
# get latent variable
U.all = matrix(runif(2 * n), ncol = 1)
# simulate y from latent variable
y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
# simulate X: 
epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
a1 = theta.value * ilrtrans.true$ilr.trans[-p] 
#   alpha1j = {
#     c1=theta*ilr.const/k+   if j \in I+
#     -c2=-theta*ilr.const/k-  if j \in I-
#     0                       o/w
#   }
alrXj.all = a0 + U.all %*% t(a1) + epsj.all #log(Xj/Xp) =alpha0j+alpha1j*U+epsj
X.all <- alrinv(alrXj.all)
colnames(X.all) = paste0('s', 1:p)

# subset out training and test sets
X = X.all[1:n, ]
X.test = X.all[-(1:n), ]
Y <- y.all[1:n]
Y.test <- y.all[-(1:n)]

# about beta
non0.beta = as.vector(SBP.true != 0)
is0.beta = !non0.beta
# solve for beta
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
beta = (b1 / c1plusc2) * theta.value * as.vector(ilrtrans.true$ilr.trans)

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

# ##############################################################################
# # Dr. Ma's new slr
# ##############################################################################
# 
# # apply supervised log-ratios, using CV to select threshold
# slrnew = slr(x = X, y = Y)
# slrnew_activevars = names(slrnew$index)
# slrnew_SBP = matrix(slrnew$index)
# rownames(slrnew_SBP) = slrnew_activevars
# 
# slrnew.coefs = coefficients(slrnew$model)
# slrnew.a0 = slrnew.coefs[1]
# slrnew.thetahat = slrnew.coefs[-1]
# 
# slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
# slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
# rownames(slrnew.betahat) = colnames(X)
# slrnew.betahat[slrnew_activevars, ] =
#   as.numeric(slrnew.betahat.nonzero)
# 
# # compute metrics on the selected model #
# slrnew.metrics = getMetricsBalanceReg(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
#   ilrX.test = getIlrX(X.test[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
#   n.train = n, n.test = n,
#   thetahat0 = slrnew.a0, thetahat = slrnew.thetahat, betahat = slrnew.betahat,
#   true.sbp = SBP.true,
#   is0.true.beta = is0.beta, non0.true.beta = non0.beta,
#   metrics = c("prediction", "selection"))

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

