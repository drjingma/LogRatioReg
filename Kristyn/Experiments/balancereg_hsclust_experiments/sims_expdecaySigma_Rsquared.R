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
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "pminus4"
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
  rho = 0.5 #
  desired_Rsquared = 0.8 #
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
  #################
  
  SigmaW = rgExpDecay(p, rho)$Sigma
  # fields::image.plot(SigmaW)
  
  # theta settings
  SigmaW_hsclust = HSClust(
    W = getSimilarityMatrix(unnormalized_similarity_matrix = SigmaW),
    levelMax = p, force_levelMax = TRUE)
  SBP = sbp.fromHSClust(levels_matrix = SigmaW_hsclust$allLevels)
  
  # a preliminary plot of the tree given by covariance matrix SigmaW
  # nodes_types = data.frame(
  #   name = c(colnames(SBP), rownames(SBP)),
  #   type = c(rep("balance", ncol(SBP)), rep("covariate", nrow(SBP)))
  # )
  # plotSBP(SBP, title = "Sigma", nodes_types = nodes_types) 
  
  # for each column (contrast), find which variables are included (1 or -1)
  indices.theta = p - 4
  
  # print(indices.theta)
  # error checking indices.theta found based on theta.settings argument
  if(is.null(indices.theta)){
    stop("invalid indices.theta")
  }
  
  # get theta
  theta = rep(0, p - 1)
  if(is.null(values.theta)){
    # assume values.theta = 1
    values.theta = 1
  }
  if(length(values.theta) == 1){
    # if values.theta = 1, assume it's the same value for all nonzero indices
    values.theta = rep(values.theta, length(indices.theta))
  } else if(length(values.theta) != length(indices.theta)){
    # when 1 < length(values.theta) < total # of nonzero values
    stop("indices.theta does not have same length as values.theta")
  }
  theta[indices.theta] = values.theta
  theta = as.matrix(theta)
  
  # about beta
  beta = as.vector(getBeta(theta, sbp = SBP))
  names(beta) <- paste0('s', 1:p)
  non0.beta = (beta != 0)
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
  # Population parameters, continued
  muW = c(rep(log(p), 5), rep(0, p - 5))
  names(muW) = names(beta)
  
  # plot the tree given by covariance matrix SigmaW, indicating 
  #   significant covariates and balances (theta's)
  leaf_types = rep("insignif cov", nrow(SBP))
  leaf_types[non0.beta] = "signif cov"
  balance_types = rep("insignif bal", ncol(SBP))
  balance_types[theta[, 1] != 0] = "signif bal"
  nodes_types = data.frame(
    name = c(colnames(SBP), rownames(SBP)),
    type = c(balance_types, leaf_types)
  )
  # plotSBP(SBP, title = "Sigma", nodes_types = nodes_types)
  
  file.end = paste0(
    "_", sigma.settings,
    "_", theta.settings, 
    "_val", values.theta[1],
    "_dim", n, "x", p, 
    "_Rsq", desired_Rsquared,
    "_rho", rho, 
    "_noise", sigma_eps,
    "_int", intercept,
    "_scale", scaling,
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # seed = 123
  # set.seed(seed)
  # generate X
  logW.all <- mvrnorm(n = 2 * n, mu = muW, Sigma = SigmaW) 
  W.all <- exp(logW.all)
  X.all <- sweep(W.all, 1, rowSums(W.all), FUN='/')
  colnames(X.all) = names(beta)
  # create the ilr(X.all) covariate by hand to
  #   generate y
  SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  U.true.details = getUdetailed(sbp = SBP.true)
  U.true = U.true.details$U
  y.all = as.numeric(log(X.all) %*% U.true * values.theta) + rnorm(n) * sigma_eps
  
  # subset out training and test sets
  X = X.all[1:n, ]
  X.test = X.all[-(1:n), ]
  Y <- y.all[1:n]
  Y.test <- y.all[-(1:n)]
  
  ##############################################################################
  # estimate Rsquared, given the true model
  SSres = sum((y.all - as.numeric(log(X.all) %*% U.true * values.theta))^2)
  SStot = sum((y.all - mean(y.all))^2)
  Rsq = 1 - SSres/SStot
  
  ##############################################################################
  # supervised log-ratios (a balance regression method) withOUT eta
  #   -- hierarchical clustering
  ##############################################################################
  # apply hierarchical clustering to the SLR distance matrix
  slrDistMat = getSupervisedMatrix(
    y = Y, X = X, rho.type = rho.type, type = "distance")
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
  
  # slrhc2_sbp = slrhc2$sbp_thresh[[slrhc2$min.idx[2]]]
  # slrhc2_thetahat = slrhc2$theta[[slrhc2$min.idx[2]]][, slrhc2$min.idx[1]]
  # slrhc2_betahat = getBeta(theta = slrhc2_thetahat, sbp = slrhc2_sbp)[, 1]
  # slrhc2_leaf_types = rep("not-selected cov", nrow(slrhc2_sbp))
  # slrhc2_leaf_types[slrhc2_betahat != 0] = "selected cov"
  # slrhc2_balance_types = rep("not-selected bal", ncol(slrhc2_sbp))
  # slrhc2_balance_types[slrhc2_thetahat != 0] = "selected bal"
  # slrhc2_nodes_types = data.frame(
  #   name = c(colnames(slrhc2_sbp), rownames(slrhc2_sbp)),
  #   type = c(slrhc2_balance_types, slrhc2_leaf_types)
  # )
  # plotSBP(slrhc2_sbp, title = "SLR hclust - eta tree", 
  #         nodes_types = slrhc2_nodes_types)
  
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
    y = Y, X = X, rho.type = rho.type, type = "similarity")
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
  
  # slrhsc2_sbp = slrhsc2$sbp_thresh[[slrhsc2$min.idx[2]]]
  # slrhsc2_thetahat = slrhsc2$theta[[slrhsc2$min.idx[2]]][, slrhsc2$min.idx[1]]
  # slrhsc2_betahat = getBeta(theta = slrhsc2_thetahat, sbp = slrhsc2_sbp)[, 1]
  # slrhsc2_leaf_types = rep("not-selected cov", nrow(slrhsc2_sbp))
  # slrhsc2_leaf_types[slrhsc2_betahat != 0] = "selected cov"
  # slrhsc2_balance_types = rep("not-selected bal", ncol(slrhsc2_sbp))
  # slrhsc2_balance_types[slrhsc2_thetahat != 0] = "selected bal"
  # slrhsc2_nodes_types = data.frame(
  #   name = c(colnames(slrhsc2_sbp), rownames(slrhsc2_sbp)),
  #   type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
  # )
  # plotSBP(slrhsc2_sbp, title = "SLR HSClust - eta tree", 
  #         nodes_types = slrhsc2_nodes_types)
  
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
  
  saveRDS(c(
    pr.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/propr_metrics", file.end))
  
  ### fin ###
}


