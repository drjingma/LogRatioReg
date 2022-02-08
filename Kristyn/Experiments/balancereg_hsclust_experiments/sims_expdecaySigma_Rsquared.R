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
nworkers = detectCores() / 2
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
  rho = 0.2 #
  desired_Rsquared = 0.8 #
  sigma_eps = get_sigma_eps(
    theta_val = values.theta, Rsq_val = desired_Rsquared, rho_val = rho)
  
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
  U.true = getIlrTrans(sbp = SBP.true)
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
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical clustering
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical clustering to the SLR distance matrix
  slrDistMat = getSlrMatrix(
    y = Y, X = X, type = "distance")
  slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
  slrhc_SBP = sbp.fromHclust(slrhc_btree)
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhc = cvBMLasso(y = Y, X = X, sbp = slrhc_SBP, nlam = nlam,
             nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  slrhc.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

  slrhc.lam.min.idx = which.min(slrhc$cvm)
  slrhc.a0 = slrhc$theta0[slrhc.lam.min.idx]
  slrhc.thetahat = slrhc$theta[, slrhc.lam.min.idx]
  slrhc.betahat = getBetaFromTheta(slrhc.thetahat, sbp = slrhc$sbp)

  # compute metrics on the selected model #
  slrhc.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrhc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrhc$sbp),
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
    "Rsq" = Rsq,
    "time" = slrhc.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hc_metrics", file.end))

  ##############################################################################
  # supervised log-ratios (a balance regression method), using
  #   distal balances
  #   -- hierarchical clustering + distal balances
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical clustering to the SLR distance matrix
  slrDistMat = getSlrMatrix(
    y = Y, X = X, type = "distance")
  slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
  slrhc_distal_SBP = sbp.fromHclust(slrhc_btree)
  slrhc_distal_SBP = sbp.subset(slrhc_distal_SBP)
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhc_distal = cvBMLasso(y = Y, X = X, sbp = slrhc_distal_SBP, nlam = nlam,
             nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  slrhc_distal.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

  slrhc_distal.lam.min.idx = which.min(slrhc_distal$cvm)
  slrhc_distal.a0 = slrhc_distal$theta0[slrhc_distal.lam.min.idx]
  slrhc_distal.thetahat = slrhc_distal$theta[, slrhc_distal.lam.min.idx]
  slrhc_distal.betahat = getBetaFromTheta(slrhc_distal.thetahat, sbp = slrhc_distal$sbp)

  # compute metrics on the selected model #
  slrhc_distal.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrhc_distal$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrhc_distal$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrhc_distal.a0, thetahat = slrhc_distal.thetahat, betahat = slrhc_distal.betahat,
    sbp = slrhc_distal$sbp,
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

  # # plot the tree given by slr-hc, indicating significant covariates
  # slrhc_distal_leaf_types = rep("covariate", nrow(slrhc_distal$sbp))
  # slrhc_distal_balance_types = rep("balance", ncol(slrhc_distal$sbp))
  # slrhc_distal_nodes_types = data.frame(
  #   name = c(colnames(slrhc_distal$sbp), rownames(slrhc_distal$sbp)),
  #   type = c(slrhc_distal_balance_types, slrhc_distal_leaf_types)
  # )
  # plotSBP(slrhc_distal$sbp, title = "slr-hc", nodes_types = slrhc_distal_nodes_types)
  # fields::image.plot(slrDistMat)

  saveRDS(c(
    slrhc_distal.metrics,
    "betaSparsity" = bspars,
    "Rsq" = Rsq,
    "time" = slrhc_distal.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hc_distal_metrics", file.end))

  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSlrMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhsc = cvBMLasso(
    y = Y, X = X, sbp = slrhsc_SBP, nlam = nlam, nfolds = K,
    intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  slrhsc.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

  slrhsc.lam.min.idx = which.min(slrhsc$cvm)
  slrhsc.a0 = slrhsc$theta0[slrhsc.lam.min.idx]
  slrhsc.thetahat = slrhsc$theta[, slrhsc.lam.min.idx]
  slrhsc.betahat = getBetaFromTheta(slrhsc.thetahat, sbp = slrhsc$sbp)

  # compute metrics on the selected model #
  slrhsc.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrhsc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrhsc$sbp),
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
    "Rsq" = Rsq,
    "time" = slrhsc.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_metrics", file.end))

  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering + thresholding with lasso
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSlrMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhsc2 = cvBMLassoThresh(
    y = Y, X = X,
    W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
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
  slrhsc2.betahat.nonzero = getBetaFromTheta(slrhsc2.thetahat, sbp = slrhsc2.SBP)
  slrhsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhsc2.betahat) = names(beta)
  slrhsc2.betahat[slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], ] =
    as.numeric(slrhsc2.betahat.nonzero)

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
    "Rsq" = Rsq,
    "time" = slrhsc2.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_thresh_lasso_metrics", file.end))
  
  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering + thresholding with mult. lm
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSlrMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhsc3 = cvBMThresh(
    y = Y, X = X,
    W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
    hsc_method = "kmeans", # "shimalik", "kmeans"
    multiple_balances = TRUE,
    force_levelMax = TRUE,
    sbp = slrhsc_SBP,
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL,
    intercept = intercept,
    standardize = scaling
  )
  end.time = Sys.time()
  slrhsc3.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrhsc3.eta.min.idx = slrhsc3$min.idx
  slrhsc3.a0 = slrhsc3$theta0[[slrhsc3.eta.min.idx]]
  slrhsc3.thetahat = slrhsc3$theta[[slrhsc3.eta.min.idx]]
  slrhsc3.SBP = slrhsc3$sbp_thresh[[slrhsc3.eta.min.idx]]
  slrhsc3.betahat.nonzero = getBetaFromTheta(slrhsc3.thetahat, sbp = slrhsc3.SBP)
  slrhsc3.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhsc3.betahat) = names(beta)
  slrhsc3.betahat[slrhsc3$meets_threshold[[slrhsc3.eta.min.idx]], ] =
    as.numeric(slrhsc3.betahat.nonzero)
  
  # compute metrics on the selected model #
  slrhsc3.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(
      X[, slrhsc3$meets_threshold[[slrhsc3.eta.min.idx]], drop = FALSE],
      sbp = slrhsc3.SBP),
    ilrX.test = getIlrX(
      X.test[, slrhsc3$meets_threshold[[slrhsc3.eta.min.idx]], drop = FALSE],
      sbp = slrhsc3.SBP),
    n.train = n, n.test = n,
    thetahat0 = slrhsc3.a0, thetahat = slrhsc3.thetahat,
    betahat = slrhsc3.betahat,
    sbp = slrhsc3.SBP,
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  # # plot the tree given by slr-hsc, indicating significant covariates
  # slrhsc3_leaf_types = rep("covariate", nrow(slrhsc3.SBP))
  # slrhsc3_balance_types = rep("balance", ncol(slrhsc3.SBP))
  # slrhsc3_nodes_types = data.frame(
  #   name = c(colnames(slrhsc3.SBP), rownames(slrhsc3.SBP)),
  #   type = c(slrhsc3_balance_types, slrhsc3_leaf_types)
  # )
  # plotSBP(slrhsc3.SBP, title = "slr-hsc-eta", nodes_types = slrhsc3_nodes_types)
  # # fields::image.plot(slrSimMat)
  
  saveRDS(c(
    slrhsc3.metrics,
    "betaSparsity" = bspars,
    "Rsq" = Rsq,
    "time" = slrhsc3.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_thresh_mlm_metrics", file.end))
  
  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering + thresholding with single lm
  #       (i.e. one balance)
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSlrMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhsc4 = cvBMThresh(
    y = Y, X = X,
    W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
    hsc_method = "kmeans", # "shimalik", "kmeans"
    multiple_balances = TRUE,
    force_levelMax = TRUE,
    sbp = slrhsc_SBP,
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL,
    intercept = intercept,
    standardize = scaling
  )
  end.time = Sys.time()
  slrhsc4.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrhsc4.eta.min.idx = slrhsc4$min.idx
  slrhsc4.a0 = slrhsc4$theta0[[slrhsc4.eta.min.idx]]
  slrhsc4.thetahat = slrhsc4$theta[[slrhsc4.eta.min.idx]]
  slrhsc4.SBP = slrhsc4$sbp_thresh[[slrhsc4.eta.min.idx]]
  slrhsc4.betahat.nonzero = getBetaFromTheta(slrhsc4.thetahat, sbp = slrhsc4.SBP)
  slrhsc4.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhsc4.betahat) = names(beta)
  slrhsc4.betahat[slrhsc4$meets_threshold[[slrhsc4.eta.min.idx]], ] =
    as.numeric(slrhsc4.betahat.nonzero)
  
  # compute metrics on the selected model #
  slrhsc4.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(
      X[, slrhsc4$meets_threshold[[slrhsc4.eta.min.idx]], drop = FALSE],
      sbp = slrhsc4.SBP),
    ilrX.test = getIlrX(
      X.test[, slrhsc4$meets_threshold[[slrhsc4.eta.min.idx]], drop = FALSE],
      sbp = slrhsc4.SBP),
    n.train = n, n.test = n,
    thetahat0 = slrhsc4.a0, thetahat = slrhsc4.thetahat,
    betahat = slrhsc4.betahat,
    sbp = slrhsc4.SBP,
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  # # plot the tree given by slr-hsc, indicating significant covariates
  # slrhsc4_leaf_types = rep("covariate", nrow(slrhsc4.SBP))
  # slrhsc4_balance_types = rep("balance", ncol(slrhsc4.SBP))
  # slrhsc4_nodes_types = data.frame(
  #   name = c(colnames(slrhsc4.SBP), rownames(slrhsc4.SBP)),
  #   type = c(slrhsc4_balance_types, slrhsc4_leaf_types)
  # )
  # plotSBP(slrhsc4.SBP, title = "slr-hsc-eta", nodes_types = slrhsc4_nodes_types)
  # # fields::image.plot(slrSimMat)
  
  saveRDS(c(
    slrhsc4.metrics,
    "betaSparsity" = bspars,
    "Rsq" = Rsq,
    "time" = slrhsc4.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_thresh_1lm_metrics", file.end))

  ##############################################################################
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  start.time = Sys.time()
  classo = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

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
    "Rsq" = Rsq,
    "time" = cl.timing
  ),
  paste0(output_dir, "/metrics", "/classo_metrics", file.end))

  ##############################################################################
  # propr method (a balance regression method)
  ##############################################################################
  start.time = Sys.time()
  pr_res <- suppressMessages(propr(X, metric = "phs"))
  pr_btree = hclust(as.dist(pr_res@matrix),method = linkage)
  pr_SBP = sbp.fromHclust(pr_btree)
  pr = cvBMLasso(y = Y, X = X, sbp = pr_SBP, nlam = nlam,
             nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  pr.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  pr.lam.min.idx = which.min(pr$cvm)
  pr.a0 = pr$theta0[pr.lam.min.idx]
  pr.thetahat = pr$theta[, pr.lam.min.idx]
  pr.betahat = getBetaFromTheta(pr.thetahat, sbp = pr$sbp)

  # compute metrics on the selected model #
  pr.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = pr$sbp),
    ilrX.test = getIlrX(X.test, sbp = pr$sbp),
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
    "Rsq" = Rsq,
    "time" = pr.timing
  ),
  paste0(output_dir, "/metrics", "/propr_metrics", file.end))
  
  ##############################################################################
  # oracle method (a balance regression method)
  #   -- hierarchical clustering
  ##############################################################################
  start.time = Sys.time()
  or_SBP = SBP.true
  or.ilrX = getIlrX(X, sbp = or_SBP)
  oracle = lm(Y ~ or.ilrX)
  end.time = Sys.time()
  or.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  or.coefs = coefficients(oracle)
  or.a0 = or.coefs[1]
  or.thetahat = or.coefs[-1]
  or.betahat = getBetaFromTheta(or.thetahat, sbp = or_SBP)

  # compute metrics on the selected model #
  or.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = or_SBP),
    ilrX.test = getIlrX(X.test, sbp = or_SBP),
    n.train = n, n.test = n,
    thetahat0 = or.a0, thetahat = or.thetahat, betahat = or.betahat,
    sbp = or_SBP,
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

  # # plot the tree given by slr-hsc, indicating significant covariates
  # or_leaf_types = rep("covariate", nrow(or_SBP))
  # or_balance_types = rep("balance", ncol(or_SBP))
  # or_nodes_types = data.frame(
  #   name = c(colnames(or_SBP), rownames(or_SBP)),
  #   type = c(or_balance_types, or_leaf_types)
  # )
  # plotSBP(or_SBP, title = "oracle", nodes_types = or_nodes_types)

  saveRDS(c(
    or.metrics,
    "betaSparsity" = bspars,
    "Rsq" = Rsq,
    "time" = or.timing
  ),
  paste0(output_dir, "/metrics", "/oracle_metrics", file.end))
  
  # ##############################################################################
  # # selbal method (a balance regression method)
  # ##############################################################################
  # library(selbal) # masks stats::cor()
  # 
  # start.time = Sys.time()
  # X.slbl = X
  # rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
  # colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
  # slbl = selbal.cv(x = X.slbl, y = as.vector(Y), n.fold = K)
  # end.time = Sys.time()
  # slbl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # U (transformation) matrix
  # U.slbl = rep(0, p)
  # names(U.slbl) = colnames(X.slbl)
  # pba.pos = unlist(subset(
  #   slbl$global.balance, subset = Group == "NUM", select = Taxa))
  # num.pos = length(pba.pos)
  # pba.neg = unlist(subset(
  #   slbl$global.balance, subset = Group == "DEN", select = Taxa))
  # num.neg = length(pba.neg)
  # U.slbl[pba.pos] = 1 / num.pos
  # U.slbl[pba.neg] = -1 / num.neg
  # norm.const = sqrt((num.pos * num.neg) / (num.pos + num.neg))
  # U.slbl = norm.const * U.slbl
  # # check: these are equal
  # # lm(as.vector(Y) ~ log(X) %*% as.matrix(U.slbl))
  # # slbl$glm
  # slbl.thetahat = coefficients(slbl$glm)[2]
  # slbl.betahat = U.slbl %*% as.matrix(slbl.thetahat)
  # 
  # # compute metrics on the selected model #
  # # prediction errors
  # # get prediction error on training set
  # slbl.Yhat.train = predict.glm(
  #   slbl$glm, newdata = data.frame(X.slbl), type = "response")
  # slbl.PE.train = crossprod(Y - slbl.Yhat.train) / n
  # # get prediction error on test set
  # X.slbl.test = X.test
  # rownames(X.slbl.test) = paste("Sample", 1:nrow(X.slbl.test), sep = "_")
  # colnames(X.slbl.test) = paste("V", 1:ncol(X.slbl.test), sep = "")
  # slbl.Yhat.test = predict.glm(
  #   slbl$glm, newdata = data.frame(X.slbl.test), type = "response")
  # slbl.PE.test = crossprod(Y.test - slbl.Yhat.test) / n
  # # estimation accuracy
  # slbl.EA = getEstimationAccuracy(true.beta = beta, betahat = slbl.betahat)
  # slbl.EA.active = getEstimationAccuracy(
  #   true.beta = beta[non0.beta], betahat = slbl.betahat[non0.beta])
  # slbl.EA.inactive = getEstimationAccuracy(
  #   true.beta = beta[is0.beta], betahat = slbl.betahat[is0.beta])
  # # selection accuracy
  # slbl.non0.betahat = abs(slbl.betahat) > 10e-8
  # slbl.SA = getSelectionAccuracy(
  #   is0.true.beta = is0.beta, non0.true.beta = non0.beta,
  #   non0.betahat = slbl.non0.betahat)
  # 
  # slbl.metrics = c(
  #   "PEtr" = slbl.PE.train,
  #   "PEte" = slbl.PE.test,
  #   "EA1" = slbl.EA$EA1,
  #   "EA2" = slbl.EA$EA2,
  #   "EAInfty" = slbl.EA$EAInfty,
  #   "EA1Active" = slbl.EA.active$EA1,
  #   "EA2Active" = slbl.EA.active$EA2,
  #   "EAInftyActive" = slbl.EA.active$EAInfty,
  #   "EA1Inactive" = slbl.EA.inactive$EA1,
  #   "EA2Inactive" = slbl.EA.inactive$EA2,
  #   "EAInftyInactive" = slbl.EA.inactive$EAInfty,
  #   "FP" = slbl.SA$FP,
  #   "FN" = slbl.SA$FN,
  #   "TPR" = slbl.SA$TPR,
  #   "precision" = slbl.SA$precision,
  #   "Fscore" = slbl.SA$Fscore
  # )
  # saveRDS(c(
  #   slbl.metrics,
  #   "betaSparsity" = bspars,
  #   "Rsq" = Rsq,
  #   "time" = slbl.timing
  # ),
  # paste0(output_dir, "/metrics", "/selbal_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


