rm(list=ls())
# Purpose: observe different R^2 values for specified rho and sigma_eps values
# Date: 1/11/2022

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
values.theta = 1
linkage = "average"
tol = 1e-4
nlam = 100
neta = 50
intercept = TRUE
K = 10
n = 100
p = 30
scaling = TRUE
#################
# if rho = 0, 
#   sigma_eps = sqrt(2/3) => R^2 = 0.6
#   sigma_eps = sqrt(1/4) => R^2 = 0.8
# if rho = 0.2, 
#   sigma_eps = sqrt(0.7125333) => R^2 = 0.6
#   sigma_eps = sqrt(0.2672) => R^2 = 0.8
# if rho = 0.5, 
#   sigma_eps = sqrt(0.808333) => R^2 = 0.6
#   sigma_eps = sqrt(0.303125) => R^2 = 0.8
# rho = 0.5
# sigma_eps = sqrt(0.303125)
get_sigma_eps = function(theta_val, Rsq_val, rho_val){
  sigma_eps_sq.tmp = theta_val^2 * (1 - Rsq_val) / Rsq_val + 
    theta_val^2 * (1 - Rsq_val) * (rho_val^3 + 2 * rho_val^2 + 3 * rho_val) / 
    (10 * Rsq_val)
  return(sqrt(sigma_eps_sq.tmp))
}
rho = 0.2 #
desired_Rsquared = 0.6 #
sigma_eps = get_sigma_eps(
  theta_val = values.theta, Rsq_val = desired_Rsquared, rho_val = rho)
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
SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
# SBP.true = matrix(c(1, 1, -1, -1, -1, rep(0, p - 5)))
U.true.details = getUdetailed(sbp = SBP.true)
U.true = U.true.details$U
# note: there is no theta
# also note: beta is U.true * values.theta
beta = as.numeric(U.true * values.theta)
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
# supervised log-ratios (a balance regression method) with eta
#   -- hierarchical spectral clustering
##############################################################################
start.time = Sys.time()
# apply hierarchical spectral clustering to the SLR similarity matrix
slrSimMat = getSupervisedMatrix(
  y = Y, X = X, type = "similarity")

# kmeans
slrhsc_btree = HSClust(
  W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
slrhsc_SBP = sbp.fromHSClust(
  levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
plotSBP(slrhsc_SBP)
fields::image.plot(slrSimMat)

# # shi-malik
# slrhsc_btree2 = HSClust(
#   W = slrSimMat, force_levelMax = TRUE, method = "shimalik")
# slrhsc_SBP2 = sbp.fromHSClust(
#   levels_matrix = slrhsc_btree2$allLevels, row_names = names(beta))
# plotSBP(slrhsc_SBP2)

# apply supervised log-ratios, using CV to select threshold and also lambda
slrhsc2 = cvILReta(
  y = Y, X = X,
  W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
  clustering_method = "hsc",
  hsc_method = "kmeans", # "shimalik", "kmeans"
  force_levelMax = TRUE,
  sbp = slrhsc_SBP,
  lambda = NULL, nlam = nlam,
  eta = NULL, neta = neta,
  nfolds = K, foldid = NULL,
  intercept = intercept,
  standardize = scaling
)
end.time = Sys.time()
slrhsc2.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

slrhsc2.eta.min.idx = slrhsc2$min.idx[2]
slrhsc2.lam.min.idx = slrhsc2$min.idx[1]
slrhsc2.a0 = slrhsc2$theta0[[slrhsc2.eta.min.idx]][slrhsc2.lam.min.idx]
slrhsc2.thetahat = slrhsc2$theta[[slrhsc2.eta.min.idx]][, slrhsc2.lam.min.idx]
slrhsc2.SBP = slrhsc2$sbp_thresh[[slrhsc2.eta.min.idx]]
slrhsc2.betahat.nonzero = getBeta(slrhsc2.thetahat, sbp = slrhsc2.SBP)
slrhsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrhsc2.betahat) = names(beta)
slrhsc2.betahat[slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], ] =
  as.numeric(slrhsc2.betahat.nonzero)

# compute metrics on the selected model #
slrhsc2.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test,
  ilrX.train = computeBalances(
    X[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
    sbp = slrhsc2.SBP),
  ilrX.test = computeBalances(
    X.test[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
    sbp = slrhsc2.SBP),
  n.train = n, n.test = n,
  thetahat0 = slrhsc2.a0, thetahat = slrhsc2.thetahat,
  betahat = slrhsc2.betahat,
  sbp = slrhsc2.SBP,
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

# plot the tree given by slr-hsc, indicating significant covariates
slrhsc2_leaf_types = rep("covariate", nrow(slrhsc2.SBP))
slrhsc2_balance_types = rep("balance", ncol(slrhsc2.SBP))
slrhsc2_nodes_types = data.frame(
  name = c(colnames(slrhsc2.SBP), rownames(slrhsc2.SBP)),
  type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
)
plotSBP(slrhsc2.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_nodes_types)
fields::image.plot(slrSimMat)
# note: centering X and y didn't change the slrSimMat
# note 2: using clrXj - clrXk instead of logXj - logXk didn't change slrSimMat

# what does Cor(y, Xj) look like? 
#   -- looks like it would work just as well for thresholding
cor_logx_y = apply(log(X), 2, function(x) (stats::cor(x, Y))^2)
plot(1:30, cor_logx_y)

# heat map of cross-validated mse's
fields::image.plot(slrhsc2$cvm)
# points(slrhsc2$cvm[slrhsc2$min.idx], col = 1)
title(
  # columns
  ylab = paste0( # eta: there is at least one value less than eta
    "eta ", "[", 
    round(slrhsc2$eta[1], 3), ",", 
    round(slrhsc2$eta[neta], 3), "]"),
  # rows
  xlab = paste0(
    "lambda ", "[", 
    round(slrhsc2$lambda[1], 3), ",", 
    round(slrhsc2$lambda[nlam], 3), "]")
)
# the more variables you have (small lambda, large eta), the more 

# try incorporating cor(logXj, y)^2?
getSupervisedMatrix2 = function(
  y, X, rho.type = "square", type = "similarity"
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  cor_logx_y = apply(log(X), 2, function(x) (stats::cor(x, Y))^2)
  
  # checks
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(0, p, p) # diagonal == 1
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      if(rho.type == "square" | rho.type == "squared" | rho.type == "s" | 
         rho.type == 2){
        val = (stats::cor(Zjk, y))^2 * 
          min((cor_logx_y[j])^2, (cor_logx_y[k])^2)
      } else{
        val = abs(stats::cor(Zjk, y)) * 
          min(abs(cor_logx_y[j]), abs(cor_logx_y[k]))
      }
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # give the rows and columns the names of taxa in X, for sbp.fromHclust()
  rownames(cormat) = colnames(X)
  colnames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  if(type != "similarity"){
    cormat = 1 - cormat
  } 
  return(cormat)
}
slrSimMat2 = getSupervisedMatrix2(
  y = Y, X = X, type = "similarity")
fields::image.plot(slrSimMat)
fields::image.plot(slrSimMat2)
# apply(slrSimMat2, 1, function(row) !all(row < 8.182383e-05))

# try incorporating balances in the same part of the ratio?
getSupervisedMatrix3 = function(
  y, X, rho.type = "square", type = "similarity"
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  cor_logx_y = apply(log(X), 2, function(x) (stats::cor(x, Y))^2)
  
  # checks
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(0, p, p) # diagonal == 1
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      Zjk2 = log(X[, j]) + log(X[, k])
      if(rho.type == "square" | rho.type == "squared" | rho.type == "s" | 
         rho.type == 2){
        val = (stats::cor(Zjk, y))^2 + (stats::cor(Zjk2, y))^2
      } else{
        val = abs(stats::cor(Zjk, y)) + abs(stats::cor(Zjk2, y))
      }
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # give the rows and columns the names of taxa in X, for sbp.fromHclust()
  rownames(cormat) = colnames(X)
  colnames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  if(type != "similarity"){
    cormat = 1 - cormat
  } 
  return(cormat)
}
slrSimMat3 = getSupervisedMatrix3(
  y = Y, X = X, type = "similarity")
fields::image.plot(slrSimMat)
fields::image.plot(slrSimMat3)

# try incorporating balances in the same part of the ratio AND cor(Xj, y)?
getSupervisedMatrix4 = function(
  y, X, rho.type = "square", type = "similarity"
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  cor_logx_y = apply(log(X), 2, function(x) (stats::cor(x, Y))^2)
  
  # checks
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(0, p, p) # diagonal == 1
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      Zjk2 = log(X[, j]) + log(X[, k])
      if(rho.type == "square" | rho.type == "squared" | rho.type == "s" | 
         rho.type == 2){
        val = max((stats::cor(Zjk, y))^2 + (stats::cor(Zjk2, y))^2 - 
          max((cor_logx_y[j])^2, (cor_logx_y[k])^2), 0)
      } else{
        val = max(abs(stats::cor(Zjk, y)) + abs(stats::cor(Zjk2, y)) - 
          max(abs(cor_logx_y[j]), abs(cor_logx_y[k])), 0)
      }
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # give the rows and columns the names of taxa in X, for sbp.fromHclust()
  rownames(cormat) = colnames(X)
  colnames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  if(type != "similarity"){
    cormat = 1 - cormat
  } 
  return(cormat)
}
slrSimMat4 = getSupervisedMatrix4(
  y = Y, X = X, type = "similarity")
fields::image.plot(slrSimMat)
fields::image.plot(slrSimMat2) # incorporate cor(logXj, y)
fields::image.plot(slrSimMat3) # incorporate balances in same side of ratio
fields::image.plot(slrSimMat4) # incorporate both

##############################################################################
# supervised log-ratios (a balance regression method) with eta
#   -- hierarchical spectral clustering
##############################################################################
start.time = Sys.time()
# apply hierarchical spectral clustering to the SLR similarity matrix
slr2SimMat4 = getSupervisedMatrix4(
  y = Y, X = X, type = "similarity")
slr2hsc_btree = HSClust(
  W = slr2SimMat4, force_levelMax = TRUE, method = "kmeans")
slr2hsc_SBP = sbp.fromHSClust(
  levels_matrix = slr2hsc_btree$allLevels, row_names = names(beta))
# apply supervised log-ratios, using CV to select threshold and also lambda
slr2hsc2 = cvILReta(
  y = Y, X = X,
  W = slr2SimMat4, # normalized similarity matrix (all values between 0 & 1)
  clustering_method = "hsc",
  hsc_method = "kmeans", # "shimalik", "kmeans"
  force_levelMax = TRUE,
  sbp = slr2hsc_SBP,
  lambda = NULL, nlam = nlam,
  eta = NULL, neta = neta,
  nfolds = K, foldid = NULL,
  intercept = intercept,
  standardize = scaling
)
end.time = Sys.time()
slr2hsc2.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

slr2hsc2.eta.min.idx = slr2hsc2$min.idx[2]
slr2hsc2.lam.min.idx = slr2hsc2$min.idx[1]
slr2hsc2.a0 = slr2hsc2$theta0[[slr2hsc2.eta.min.idx]][slr2hsc2.lam.min.idx]
slr2hsc2.thetahat = slr2hsc2$theta[[slr2hsc2.eta.min.idx]][, slr2hsc2.lam.min.idx]
slr2hsc2.SBP = slr2hsc2$sbp_thresh[[slr2hsc2.eta.min.idx]]
slr2hsc2.betahat.nonzero = getBeta(slr2hsc2.thetahat, sbp = slr2hsc2.SBP)
slr2hsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slr2hsc2.betahat) = names(beta)
slr2hsc2.betahat[slr2hsc2$meets_threshold[[slr2hsc2.eta.min.idx]], ] =
  as.numeric(slr2hsc2.betahat.nonzero)

# compute metrics on the selected model #
slr2hsc2.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test,
  ilrX.train = computeBalances(
    X[, slr2hsc2$meets_threshold[[slr2hsc2.eta.min.idx]], drop = FALSE],
    sbp = slr2hsc2.SBP),
  ilrX.test = computeBalances(
    X.test[, slr2hsc2$meets_threshold[[slr2hsc2.eta.min.idx]], drop = FALSE],
    sbp = slr2hsc2.SBP),
  n.train = n, n.test = n,
  thetahat0 = slr2hsc2.a0, thetahat = slr2hsc2.thetahat,
  betahat = slr2hsc2.betahat,
  sbp = slr2hsc2.SBP,
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

# plot the tree given by slr-hsc, indicating significant covariates
slr2hsc2_leaf_types = rep("covariate", nrow(slr2hsc2.SBP))
slr2hsc2_balance_types = rep("balance", ncol(slr2hsc2.SBP))
slr2hsc2_nodes_types = data.frame(
  name = c(colnames(slr2hsc2.SBP), rownames(slr2hsc2.SBP)),
  type = c(slr2hsc2_balance_types, slr2hsc2_leaf_types)
)
plotSBP(slr2hsc2.SBP, title = "slr2-hsc-eta", nodes_types = slr2hsc2_nodes_types)
