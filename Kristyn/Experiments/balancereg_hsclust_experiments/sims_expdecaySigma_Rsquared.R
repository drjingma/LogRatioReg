# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 1/2/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_hsclust_experiments/outputs"

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Other simulation settings
numSims = 100

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSims
) %dorng% {
  # rm(list=ls())
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
  
  # Sigma with exponential decay ###############################################
  
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
  get_sigma_eps = function(theta_val, Rsq_val, rho_val){
    sigma_eps_sq.tmp = theta_val^2 * (1 - Rsq_val) / Rsq_val + 
      theta_val^2 * (1 - Rsq_val) * (rho_val^3 + 2 * rho_val^2 + 3 * rho_val) / 
      (10 * Rsq_val)
    return(sqrt(sigma_eps_sq.tmp))
  }
  rho = 0 #
  desired_Rsquared = 0.6 #
  if(rho == 0){
    if(desired_Rsquared == 0.6){
      sigma_eps = sqrt(2/3)
    } else if(desired_Rsquared == 0.8){
      sigma_eps = sqrt(1/4)
    }
  } else if(rho == 0.2){
    if(desired_Rsquared == 0.6){
      sigma_eps = sqrt(0.7125333)
    } else if(desired_Rsquared == 0.8){
      sigma_eps = sqrt(0.2672)
    }
  } else if(rho == 0.5){
    if(desired_Rsquared == 0.6){
      sigma_eps = sqrt(0.808333)
    } else if(desired_Rsquared == 0.8){
      sigma_eps = sqrt(0.303125)
    }
  }
  
  # Population parameters
  SigmaW = rgExpDecay(p, rho)$Sigma
  muW = c(rep(log(p), 5), rep(0, p - 5))
  names(muW) = paste0('s', 1:p)
  
  file.end = paste0(
    "_", sigma.settings,
    "_dim", n, "x", p, 
    "_Rsq", desired_Rsquared,
    "_rho", rho, 
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # set.seed(1947)
  # generate X
  logW.all <- mvrnorm(n = 2 * n, mu = muW, Sigma = SigmaW) 
  W.all <- exp(logW.all)
  X.all <- sweep(W.all, 1, rowSums(W.all), FUN='/')
  colnames(X.all) = paste0('s', 1:p)
  # create the ilr(X.all) covariate by hand to
  #   generate y
  SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
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
  # supervised log-ratios (a balance regression method) withOUT eta
  #   -- hierarchical clustering
  ##############################################################################
  # apply hierarchical clustering to the SLR distance matrix
  slrDistMat = getSupervisedMatrix(
    y = Y, X = X, type = "distance")
  slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
  slrhc_SBP = sbp.fromHclust(slrhc_btree)
  
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhc = cvILR(y = Y, X = X, sbp = slrhc_SBP, nlam = nlam, 
             nfolds = K, intercept = intercept, standardize = scaling)
  slrhc.lam.min.idx = which.min(slrhc$cvm)
  slrhc.a0 = slrhc$int[slrhc.lam.min.idx]
  slrhc.thetahat = slrhc$bet[, slrhc.lam.min.idx]
  slrhc.betahat = getBeta(slrhc.thetahat, sbp = slrhc$sbp)
  
  # compute metrics on the selected model #
  slrhc.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test, 
    ilrX.train = computeBalances(X, sbp = slrhc$sbp), 
    ilrX.test = computeBalances(X.test, sbp = slrhc$sbp), 
    n.train = n, n.test = n, 
    thetahat0 = slrhc.a0, thetahat = slrhc.thetahat, betahat = slrhc.betahat, 
    sbp = slrhc$sbp, 
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  # # plot the tree given by slr-hc, indicating significant covariates
  # slrhc_leaf_types = rep("covariate", nrow(slrhc$sbp))
  # slrhc_balance_types = rep("balance", ncol(slrhc$sbp))
  # slrhc_nodes_types = data.frame(
  #   name = c(colnames(slrhc$sbp), rownames(slrhc$sbp)),
  #   type = c(slrhc_balance_types, slrhc_leaf_types)
  # )
  # plotSBP(slrhc$sbp, title = "slr-hc", nodes_types = slrhc_nodes_types)
  # fields::image.plot(slrDistMat)
  
  saveRDS(c(
    slrhc.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/slr_hc_metrics", file.end))
  
  ##############################################################################
  # supervised log-ratios (a balance regression method) with eta
  #   -- hierarchical clustering
  ##############################################################################
  # apply hierarchical clustering to the SLR distance matrix, using thresholding
  #   parameter eta to select covariates to include in the binary tree
  
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhc2 = cvILReta(
    y = Y, X = X,
    W = slrDistMat,
    clustering_method = "hc",
    force_levelMax = TRUE,
    sbp = slrhc_SBP,
    lambda = NULL, nlam = nlam,
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL,
    intercept = intercept,
    standardize = scaling
  )
  
  slrhc2.eta.min.idx = slrhc2$min.idx[2]
  slrhc2.lam.min.idx = slrhc2$min.idx[1]
  slrhc2.a0 = slrhc2$theta0[[slrhc2.eta.min.idx]][slrhc2.lam.min.idx]
  slrhc2.thetahat = slrhc2$theta[[slrhc2.eta.min.idx]][, slrhc2.lam.min.idx]
  slrhc2.SBP = slrhc2$sbp_thresh[[slrhc2.eta.min.idx]]
  slrhc2.betahat.nonzero = getBeta(slrhc2.thetahat, sbp = slrhc2.SBP)
  slrhc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhc2.betahat) = names(beta)
  slrhc2.betahat[slrhc2$meets_threshold[[slrhc2.eta.min.idx]], ] = 
    as.numeric(slrhc2.betahat.nonzero)
  
  # compute metrics on the selected model #
  slrhc2.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test, 
    ilrX.train = computeBalances(
      X[, slrhc2$meets_threshold[[slrhc2.eta.min.idx]], drop = FALSE], 
      sbp = slrhc2.SBP), 
    ilrX.test = computeBalances(
      X.test[, slrhc2$meets_threshold[[slrhc2.eta.min.idx]], drop = FALSE], 
      sbp = slrhc2.SBP), 
    n.train = n, n.test = n, 
    thetahat0 = slrhc2.a0, thetahat = slrhc2.thetahat, betahat = slrhc2.betahat, 
    sbp = slrhc2.SBP, 
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  # plot the tree given by slr-hc-eta, indicating significant covariates
  slrhc2_leaf_types = rep("covariate", nrow(slrhc2.SBP))
  slrhc2_balance_types = rep("balance", ncol(slrhc2.SBP))
  slrhc2_nodes_types = data.frame(
    name = c(colnames(slrhc2.SBP), rownames(slrhc2.SBP)),
    type = c(slrhc2_balance_types, slrhc2_leaf_types)
  )
  plotSBP(slrhc2.SBP, title = "slr-hc-eta", nodes_types = slrhc2_nodes_types)
  # fields::image.plot(slrDistMat)
  
  saveRDS(c(
    slrhc2.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/slr_hc_eta_metrics", file.end))
  
  ##############################################################################
  # supervised log-ratios (a balance regression method) withOUT eta
  #   -- hierarchical spectral clustering
  ##############################################################################
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSupervisedMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
  
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhsc = cvILR(
    y = Y, X = X, sbp = slrhsc_SBP, nlam = nlam, nfolds = K, 
    intercept = intercept, standardize = scaling)
  slrhsc.lam.min.idx = which.min(slrhsc$cvm)
  slrhsc.a0 = slrhsc$int[slrhsc.lam.min.idx]
  slrhsc.thetahat = slrhsc$bet[, slrhsc.lam.min.idx]
  slrhsc.betahat = getBeta(slrhsc.thetahat, sbp = slrhsc$sbp)
  
  # compute metrics on the selected model #
  slrhsc.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test, 
    ilrX.train = computeBalances(X, sbp = slrhsc$sbp), 
    ilrX.test = computeBalances(X.test, sbp = slrhsc$sbp), 
    n.train = n, n.test = n, 
    thetahat0 = slrhsc.a0, thetahat = slrhsc.thetahat, betahat = slrhsc.betahat, 
    sbp = slrhsc$sbp, 
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  # # plot the tree given by slr-hsc, indicating significant covariates
  # slrhsc_leaf_types = rep("covariate", nrow(slrhsc$sbp))
  # slrhsc_balance_types = rep("balance", ncol(slrhsc$sbp))
  # slrhsc_nodes_types = data.frame(
  #   name = c(colnames(slrhsc$sbp), rownames(slrhsc$sbp)),
  #   type = c(slrhsc_balance_types, slrhsc_leaf_types)
  # )
  # plotSBP(slrhsc$sbp, title = "slr-hsc", nodes_types = slrhsc_nodes_types)
  # fields::image.plot(slrSimMat)
  
  saveRDS(c(
    slrhsc.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/slr_hsc_metrics", file.end))
  
  ##############################################################################
  # supervised log-ratios (a balance regression method) with eta
  #   -- hierarchical spectral clustering
  ##############################################################################
  
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
  
  # # plot the tree given by slr-hsc, indicating significant covariates
  # slrhsc2_leaf_types = rep("covariate", nrow(slrhsc2.SBP))
  # slrhsc2_balance_types = rep("balance", ncol(slrhsc2.SBP))
  # slrhsc2_nodes_types = data.frame(
  #   name = c(colnames(slrhsc2.SBP), rownames(slrhsc2.SBP)),
  #   type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
  # )
  # plotSBP(slrhsc2.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_nodes_types)
  # # fields::image.plot(slrSimMat)
  
  saveRDS(c(
    slrhsc2.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/slr_hsc_eta_metrics", file.end))
  
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
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  saveRDS(c(
    cl.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/classo_metrics", file.end))
  
  ##############################################################################
  # propr method (a balance regression method)
  ##############################################################################
  
  pr_res <- suppressMessages(propr(X, metric = "phs"))
  pr_btree = hclust(as.dist(pr_res@matrix),method = linkage)
  pr_SBP = sbp.fromHclust(pr_btree)
  pr = cvILR(y = Y, X = X, sbp = pr_SBP, nlam = nlam, 
             nfolds = K, intercept = intercept, standardize = scaling)
  
  pr.lam.min.idx = which.min(pr$cvm)
  pr.a0 = pr$int[pr.lam.min.idx]
  pr.thetahat = pr$bet[, pr.lam.min.idx]
  pr.betahat = getBeta(pr.thetahat, sbp = pr$sbp)
  
  # compute metrics on the selected model #
  pr.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test, 
    ilrX.train = computeBalances(X, sbp = pr$sbp), 
    ilrX.test = computeBalances(X.test, sbp = pr$sbp), 
    n.train = n, n.test = n, 
    thetahat0 = pr.a0, thetahat = pr.thetahat, betahat = pr.betahat, 
    sbp = pr$sbp, 
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  # # plot the tree given by slr-hsc, indicating significant covariates
  # pr_leaf_types = rep("covariate", nrow(pr$sbp))
  # pr_balance_types = rep("balance", ncol(pr$sbp))
  # pr_nodes_types = data.frame(
  #   name = c(colnames(pr$sbp), rownames(pr$sbp)),
  #   type = c(pr_balance_types, pr_leaf_types)
  # )
  # plotSBP(pr$sbp, title = "propr", nodes_types = pr_nodes_types)
  # fields::image.plot(pr_res@matrix)
  
  saveRDS(c(
    pr.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/propr_metrics", file.end))
  
  ### fin ###
}


