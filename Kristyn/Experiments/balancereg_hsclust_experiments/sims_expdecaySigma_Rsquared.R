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
  rho = 0
  desired_Rsquared = 0.6
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
  # supervised log-ratios (a balance regression method) 
  #   -- hierarchical spectral clustering
  ##############################################################################
  
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrMat = getSupervisedMatrix(y = Y, X = X, rho.type = rho.type)
  # fields::image.plot(slrMat)
  slr.hsclust = HSClust(
    W = slrMat, 
    force_levelMax = TRUE, method = "kmeans")
  slrhsc.SBP_full = sbp.fromHSClust(
    levels_matrix = slr.hsclust$allLevels, row_names = names(beta))
  slrhsc = cvILReta(
    y = Y, X = X, 
    W = slrMat, # normalized similarity matrix (all values between 0 & 1)
    hsc_method = "kmeans", # "shimalik", "kmeans"
    force_levelMax = TRUE, 
    sbp = slrhsc.SBP_full,
    lambda = NULL, nlam = nlam, 
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL, 
    intercept = intercept, 
    standardize = scaling
  )
  
  slrhsc_sbp = slrhsc$sbp_thresh[[slrhsc$min.idx[2]]]
  slrhsc_thetahat = slrhsc$theta[[slrhsc$min.idx[2]]][, slrhsc$min.idx[1]]
  slrhsc_betahat = getBeta(theta = slrhsc_thetahat, sbp = slrhsc_sbp)[, 1]
  slrhsc_leaf_types = rep("not-selected cov", nrow(slrhsc_sbp))
  slrhsc_leaf_types[slrhsc_betahat != 0] = "selected cov"
  slrhsc_balance_types = rep("not-selected bal", ncol(slrhsc_sbp))
  slrhsc_balance_types[slrhsc_thetahat != 0] = "selected bal"
  slrhsc_nodes_types = data.frame(
    name = c(colnames(slrhsc_sbp), rownames(slrhsc_sbp)),
    type = c(slrhsc_balance_types, slrhsc_leaf_types)
  )
  # plotSBP(slrhsc_sbp, title = "SLR HSClust tree", 
  #         nodes_types = slrhsc_nodes_types)
  
  slrhsc.eta.min.idx = slrhsc$min.idx[2]
  slrhsc.lam.min.idx = slrhsc$min.idx[1]
  slrhsc.a0 = slrhsc$theta0[[slrhsc.eta.min.idx]][slrhsc.lam.min.idx]
  slrhsc.thetahat = slrhsc$theta[[slrhsc.eta.min.idx]][, slrhsc.lam.min.idx]
  slrhsc.SBP = slrhsc$sbp_thresh[[slrhsc.eta.min.idx]]
  slrhsc.betahat.nonzero = getBeta(slrhsc.thetahat, sbp = slrhsc.SBP)
  slrhsc.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhsc.betahat) = names(beta)
  slrhsc.betahat[slrhsc$meets_threshold[[slrhsc.eta.min.idx]], ] = 
    as.numeric(slrhsc.betahat.nonzero)
  
  # compute metrics on the selected model #
  slrhsc.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test, 
    ilrX.train = computeBalances(
      X[, slrhsc$meets_threshold[[slrhsc.eta.min.idx]], drop = FALSE], 
      sbp = slrhsc.SBP), 
    ilrX.test = computeBalances(
      X.test[, slrhsc$meets_threshold[[slrhsc.eta.min.idx]], drop = FALSE], 
      sbp = slrhsc.SBP), 
    n.train = n, n.test = n, 
    thetahat0 = slrhsc.a0, thetahat = slrhsc.thetahat, betahat = slrhsc.betahat, 
    sbp = slrhsc.SBP, 
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

  saveRDS(c(
    slrhsc.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/slr_hsclust_metrics", file.end))
  
  ##############################################################################
  # supervised log-ratios (a balance regression method) 
  #   -- hierarchical clustering
  ##############################################################################
  
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhcDistMat = getSupervisedMatrix(
    y = Y, X = X, rho.type = rho.type, type = "distance")
  slrhc_btree = hclust(as.dist(slrhcDistMat), method = linkage)
  slrhc.SBP_full = sbp.fromHclust(slrhc_btree)
  slrhc = cvILReta(
    y = Y, X = X, 
    W = slrMat,
    hsc_method = "kmeans", 
    force_levelMax = TRUE, 
    sbp = slrhc.SBP_full,
    lambda = NULL, nlam = nlam, 
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL, 
    intercept = intercept, 
    standardize = scaling
  )

  # slrhc_sbp = slrhc$sbp_thresh[[slrhc$min.idx[2]]]
  # slrhc_thetahat = slrhc$theta[[slrhc$min.idx[2]]][, slrhc$min.idx[1]]
  # slrhc_betahat = getBeta(theta = slrhc_thetahat, sbp = slrhc_sbp)[, 1]
  # slrhc_leaf_types = rep("not-selected cov", nrow(slrhc_sbp))
  # slrhc_leaf_types[slrhc_betahat != 0] = "selected cov"
  # slrhc_balance_types = rep("not-selected bal", ncol(slrhc_sbp))
  # slrhc_balance_types[slrhc_thetahat != 0] = "selected bal"
  # slrhc_nodes_types = data.frame(
  #   name = c(colnames(slrhc_sbp), rownames(slrhc_sbp)),
  #   type = c(slrhc_balance_types, slrhc_leaf_types)
  # )
  # plotSBP(slrhc_sbp, title = "SLR hclust tree", 
  #         nodes_types = slrhc_nodes_types)
  
  slrhc.eta.min.idx = slrhc$min.idx[2]
  slrhc.lam.min.idx = slrhc$min.idx[1]
  slrhc.a0 = slrhc$theta0[[slrhc.eta.min.idx]][slrhc.lam.min.idx]
  slrhc.thetahat = slrhc$theta[[slrhc.eta.min.idx]][, slrhc.lam.min.idx]
  slrhc.SBP = slrhc$sbp_thresh[[slrhc.eta.min.idx]]
  slrhc.betahat.nonzero = getBeta(slrhc.thetahat, sbp = slrhc.SBP)
  slrhc.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhc.betahat) = names(beta)
  slrhc.betahat[slrhc$meets_threshold[[slrhc.eta.min.idx]], ] = 
    as.numeric(slrhc.betahat.nonzero)
  
  # compute metrics on the selected model #
  slrhc.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test, 
    ilrX.train = computeBalances(
      X[, slrhc$meets_threshold[[slrhc.eta.min.idx]], drop = FALSE], 
      sbp = slrhc.SBP), 
    ilrX.test = computeBalances(
      X.test[, slrhc$meets_threshold[[slrhc.eta.min.idx]], drop = FALSE], 
      sbp = slrhc.SBP), 
    n.train = n, n.test = n, 
    thetahat0 = slrhc.a0, thetahat = slrhc.thetahat, betahat = slrhc.betahat, 
    sbp = slrhc.SBP, 
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  saveRDS(c(
    slrhc.metrics, 
    "betaSparsity" = bspars,
    "Rsq" = Rsq
  ), 
  paste0(output_dir, "/metrics", "/slr_hclust_metrics", file.end))
}


