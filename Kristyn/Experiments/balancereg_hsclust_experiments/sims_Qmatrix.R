# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 1/27/2021

################################################################################
# libraries and settings

# (put in the loop)
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
# theta_overlapping_balance = FALSE # (put in the loop)

# ################################################################################
# # Defining theta for the model #
# ################################################################################
# 
# library(balance)
# source("Kristyn/Functions/supervisedlogratios.R")
# source("Kristyn/Functions/HSClust.R")
# 
# # for plots
# library(ggraph) # make dendrogram
# library(igraph) # transform dataframe to graph object: graph_from_data_frame()
# library(tidygraph)
# 
# # data
# Q = read.table(file = "Data/Q.txt")
# Q = as.matrix(Q)
# 
# # Settings to toggle with (make sure these are the same as in the loop)
# sigma.settings = "Qmatrix"
# values.theta = 1
# linkage = "average"
# tol = 1e-4
# nlam = 100
# neta = 50
# intercept = TRUE
# K = 10
# n = 100
# p = ncol(Q)
# scaling = TRUE
# #################
# rho = 0
# sigma_eps = 0
# 
# # Population parameters
# if(theta_overlapping_balance){
#   registerDoRNG(123) # for tree with an overlapped balance
# } else{
#   registerDoRNG(456) # for tree with no overlapping balances
# }
# 
# ### apply hierarchical spectral clustering to similarity matrix Q?
# # hsclust_Q = HSClust(
# #   W = Q, # assuming Q is a similarity matrix
# #   force_levelMax = TRUE, method = "kmeans")
# # sbp_Q = sbp.fromHSClust(levels_matrix = hsclust_Q$allLevels)
# ### apply hierarchical clustering to similarity matrix Q?
# hclust_Q = hclust(
#   as.dist(getSimilarityMatrix(unnormalized_similarity_matrix = Q)), 
#   method = linkage)
# sbp_Q = sbp.fromHclust(hclust_Q)
# 
# # define theta (directly related to sparsity)
# # randomly choose 5% of balances/internal nodes that at most 5 active variables
# #   i.e. excluse balances/internal nodes with more than 5 active variables
# contrast_vars = apply(sbp_Q, 2, FUN = function(col) which(col != 0))
# contrast_lens = sapply(contrast_vars, length)
# viable_contrasts = colnames(sbp_Q)[contrast_lens <= 5]
# num_select_viable_contrasts = signif(length(viable_contrasts) * 0.05, 0)
# # View(sbp_Q[, viable_contrasts])
# theta = rep(0, p - 1)
# names(theta) = colnames(sbp_Q)
# selected_viable_contrasts = sample(
#   x = viable_contrasts, size = num_select_viable_contrasts, 
#   replace = FALSE)
# theta[selected_viable_contrasts] = values.theta
# non0.theta = (theta != 0)
# which(non0.theta)
# 
# # define beta
# beta = as.vector(getBeta(theta = theta, sbp = sbp_Q))
# names(beta) <- rownames(sbp_Q)
# non0.beta = (beta != 0)
# is0.beta = abs(beta) <= 10e-8
# bspars = sum(non0.beta)
# 
# # plot the selected balances
# balance_types = rep("bal", ncol(sbp_Q))
# balance_types[non0.theta] = "selected bal"
# covariate_types = rep("cov", nrow(sbp_Q))
# covariate_types[non0.beta] = "selected cov"
# model_nodes_types = data.frame(
#   name = c(colnames(sbp_Q), rownames(sbp_Q)),
#   type = c(balance_types, covariate_types)
# )
# plotSBP(sbp = sbp_Q, nodes_types = model_nodes_types, text_size = 1.5) # 8x25
# 
# # save theta
# saveRDS(
#   theta, 
#   file = paste0(
#     output_dir, 
#     "/theta_overlapping_balance_", theta_overlapping_balance, ".rds")
# )

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSims
) %dorng% {
  # rm(list=ls())
  
  output_dir = "Kristyn/Experiments/balancereg_hsclust_experiments/outputs"
  theta_overlapping_balance = FALSE
  
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
  
  ################################################################################
  # data
  Q = read.table(file = "Data/Q.txt")
  Q = as.matrix(Q)
  
  # Sigma with exponential decay ###############################################
  
  # Settings to toggle with
  sigma.settings = "Qmatrix"
  values.theta = 1
  linkage = "average"
  tol = 1e-4
  nlam = 100
  neta = 50
  intercept = TRUE
  K = 10
  n = 100
  p = ncol(Q)
  scaling = TRUE
  #################
  rho = 0
  sigma_eps = 0
  
  # Population parameters
  
  ### apply hierarchical spectral clustering to similarity matrix Q?
  # hsclust_Q = HSClust(
  #   W = Q, # assuming Q is a similarity matrix
  #   force_levelMax = TRUE, method = "kmeans")
  # sbp_Q = sbp.fromHSClust(levels_matrix = hsclust_Q$allLevels)
  ### apply hierarchical clustering to similarity matrix Q?
  hclust_Q = hclust(
    as.dist(getSimilarityMatrix(unnormalized_similarity_matrix = Q)), 
    method = linkage)
  sbp_Q = sbp.fromHclust(hclust_Q)
  
  # define covariate parameters
  SigmaW = rgExpDecay(p, rho)$Sigma
  muW = c(rep(log(p), 5), rep(0, p - 5))
  names(muW) = paste0('s', 1:p)
  
  # import theta
  theta = readRDS(paste0(
    output_dir, 
    "/theta_overlapping_balance_", theta_overlapping_balance, ".rds")
  )
  
  # define beta
  beta = as.vector(getBeta(theta = theta, sbp = sbp_Q))
  names(beta) <- rownames(sbp_Q)
  non0.beta = (beta != 0)
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
  # saving stuff
  file.end = paste0(
    "_", sigma.settings,
    "_balanceoverlap", theta_overlapping_balance,
    "_dim", n, "x", p, 
    "_noise", sigma_eps, 
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  
  fake.data = simulateBalanceReg(
    mu = muW, Sigma = SigmaW, sbp = sbp_Q, n = 2 * n, theta = theta, 
    sigma.noise = sigma_eps)
  colnames(fake.data$X) = names(beta)
  
  # subset out training and test sets
  X = fake.data$X[1:n, ]
  X.test = fake.data$X[-(1:n), ]
  Y <- fake.data$y[1:n, , drop = TRUE]
  Y.test <- fake.data$y[-(1:n), , drop = TRUE]
  
  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical clustering
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical clustering to the SLR distance matrix
  slrDistMat = getSupervisedMatrix(
    y = Y, X = X, type = "distance")
  slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
  slrhc_SBP = sbp.fromHclust(slrhc_btree)
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhc = cvILR(y = Y, X = X, sbp = slrhc_SBP, nlam = nlam,
             nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  slrhc.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

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
  slrDistMat = getSupervisedMatrix(
    y = Y, X = X, type = "distance")
  slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
  slrhc_distal_SBP = sbp.fromHclust(slrhc_btree)
  slrhc_distal_SBP = sbp.subset(slrhc_distal_SBP)
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhc_distal = cvILR(y = Y, X = X, sbp = slrhc_distal_SBP, nlam = nlam,
             nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  slrhc_distal.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

  slrhc_distal.lam.min.idx = which.min(slrhc_distal$cvm)
  slrhc_distal.a0 = slrhc_distal$int[slrhc_distal.lam.min.idx]
  slrhc_distal.thetahat = slrhc_distal$bet[, slrhc_distal.lam.min.idx]
  slrhc_distal.betahat = getBeta(slrhc_distal.thetahat, sbp = slrhc_distal$sbp)

  # compute metrics on the selected model #
  slrhc_distal.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = computeBalances(X, sbp = slrhc_distal$sbp),
    ilrX.test = computeBalances(X.test, sbp = slrhc_distal$sbp),
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
    "time" = slrhc_distal.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hc_distal_metrics", file.end))


  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical clustering + thresholding
  ##############################################################################
  # apply hierarchical clustering to the SLR distance matrix, using thresholding
  #   parameter eta to select covariates to include in the binary tree
  start.time = Sys.time()
  # apply hierarchical clustering to the SLR distance matrix
  slrDistMat = getSupervisedMatrix(
    y = Y, X = X, type = "distance")
  slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
  slrhc_SBP = sbp.fromHclust(slrhc_btree)
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
  end.time = Sys.time()
  slrhc2.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

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

  # # plot the tree given by slr-hc-eta, indicating significant covariates
  # slrhc2_leaf_types = rep("covariate", nrow(slrhc2.SBP))
  # slrhc2_balance_types = rep("balance", ncol(slrhc2.SBP))
  # slrhc2_nodes_types = data.frame(
  #   name = c(colnames(slrhc2.SBP), rownames(slrhc2.SBP)),
  #   type = c(slrhc2_balance_types, slrhc2_leaf_types)
  # )
  # plotSBP(slrhc2.SBP, title = "slr-hc-eta", nodes_types = slrhc2_nodes_types)
  # # fields::image.plot(slrDistMat)

  saveRDS(c(
    slrhc2.metrics,
    "betaSparsity" = bspars,
    "time" = slrhc2.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hc_eta_metrics", file.end))

  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering
  ##############################################################################
  start.time = Sys.time()
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
  end.time = Sys.time()
  slrhsc.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

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
    "time" = slrhsc.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_metrics", file.end))

  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering + thresholding
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSupervisedMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
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
    "time" = slrhsc2.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_eta_metrics", file.end))

  # ##############################################################################
  # # supervised log-ratios (a balance regression method)
  # #   -- hierarchical spectral clustering + thresholding + natural stop
  # #     i.e. if spectral clustering halts, let it
  # ##############################################################################
  # start.time = Sys.time()
  # # apply hierarchical spectral clustering to the SLR similarity matrix
  # slrSimMat = getSupervisedMatrix(
  #   y = Y, X = X, type = "similarity")
  # slrhsc_btree_natstop = HSClust(
  #   W = slrSimMat, force_levelMax = FALSE, method = "kmeans")
  # slrhsc_SBP_natstop = sbp.fromHSClust(
  #   levels_matrix = slrhsc_btree_natstop$allLevels, row_names = names(beta))
  # # apply supervised log-ratios, using CV to select threshold and also lambda
  # slrhsc2_natstop = cvILReta(
  #   y = Y, X = X,
  #   W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
  #   clustering_method = "hsc",
  #   hsc_method = "kmeans", # "shimalik", "kmeans"
  #   force_levelMax = TRUE,
  #   sbp = slrhsc_SBP_natstop,
  #   lambda = NULL, nlam = nlam,
  #   eta = NULL, neta = neta,
  #   nfolds = K, foldid = NULL,
  #   intercept = intercept,
  #   standardize = scaling
  # )
  # end.time = Sys.time()
  # slrhsc2_natstop.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # slrhsc2_natstop.eta.min.idx = slrhsc2_natstop$min.idx[2]
  # slrhsc2_natstop.lam.min.idx = slrhsc2_natstop$min.idx[1]
  # slrhsc2_natstop.a0 = slrhsc2_natstop$theta0[[slrhsc2_natstop.eta.min.idx]][
  #   slrhsc2_natstop.lam.min.idx]
  # slrhsc2_natstop.thetahat = slrhsc2_natstop$theta[[
  #   slrhsc2_natstop.eta.min.idx]][, slrhsc2_natstop.lam.min.idx]
  # slrhsc2_natstop.SBP = slrhsc2_natstop$sbp_thresh[[
  #   slrhsc2_natstop.eta.min.idx]]
  # slrhsc2_natstop.betahat.nonzero = getBeta(
  #   slrhsc2_natstop.thetahat, sbp = slrhsc2_natstop.SBP)
  # slrhsc2_natstop.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  # rownames(slrhsc2_natstop.betahat) = names(beta)
  # slrhsc2_natstop.betahat[slrhsc2_natstop$meets_threshold[[
  #   slrhsc2_natstop.eta.min.idx]], ] =
  #   as.numeric(slrhsc2_natstop.betahat.nonzero)
  # 
  # # compute metrics on the selected model #
  # slrhsc2_natstop.metrics = getMetricsBalanceReg(
  #   y.train = Y, y.test = Y.test,
  #   ilrX.train = computeBalances(
  #     X[, slrhsc2_natstop$meets_threshold[[slrhsc2_natstop.eta.min.idx]],
  #       drop = FALSE],
  #     sbp = slrhsc2_natstop.SBP),
  #   ilrX.test = computeBalances(
  #     X.test[, slrhsc2_natstop$meets_threshold[[slrhsc2_natstop.eta.min.idx]],
  #            drop = FALSE],
  #     sbp = slrhsc2_natstop.SBP),
  #   n.train = n, n.test = n,
  #   thetahat0 = slrhsc2_natstop.a0, thetahat = slrhsc2_natstop.thetahat,
  #   betahat = slrhsc2_natstop.betahat,
  #   sbp = slrhsc2_natstop.SBP,
  #   true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  # 
  # # # plot the tree given by slr-hsc, indicating significant covariates
  # # slrhsc2_natstop_leaf_types = rep("covariate", nrow(slrhsc2_natstop.SBP))
  # # slrhsc2_natstop_balance_types = rep("balance", ncol(slrhsc2_natstop.SBP))
  # # slrhsc2_natstop_nodes_types = data.frame(
  # #   name = c(colnames(slrhsc2_natstop.SBP), rownames(slrhsc2_natstop.SBP)),
  # #   type = c(slrhsc2_natstop_balance_types, slrhsc2_natstop_leaf_types)
  # # )
  # # plotSBP(slrhsc2_natstop.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_natstop_nodes_types)
  # # # fields::image.plot(slrSimMat)
  # 
  # saveRDS(c(
  #   slrhsc2_natstop.metrics,
  #   "betaSparsity" = bspars,
  #   "time" = slrhsc2_natstop.timing
  # ),
  # paste0(output_dir, "/metrics", "/slr_hsc_eta_natstop_metrics", file.end))
  # 
  # ##############################################################################
  # # supervised log-ratios (a balance regression method)
  # #   -- hierarchical spectral clustering + thresholding + Newman-Girvan
  # #     modularity as stopping rule
  # ##############################################################################
  # start.time = Sys.time()
  # # apply hierarchical spectral clustering to the SLR similarity matrix
  # slrSimMat = getSupervisedMatrix(
  #   y = Y, X = X, type = "similarity")
  # slrhsc_btree_ngmstop = HSClust(
  #   W = slrSimMat, force_levelMax = FALSE, method = "kmeans",
  #   stopping_rule = "ngm")
  # slrhsc_SBP_ngmstop = sbp.fromHSClust(
  #   levels_matrix = slrhsc_btree_ngmstop$allLevels, row_names = names(beta))
  # # apply supervised log-ratios, using CV to select threshold and also lambda
  # slrhsc2_ngmstop = cvILReta(
  #   y = Y, X = X,
  #   W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
  #   clustering_method = "hsc",
  #   hsc_method = "kmeans", # "shimalik", "kmeans"
  #   force_levelMax = TRUE,
  #   sbp = slrhsc_SBP_ngmstop,
  #   lambda = NULL, nlam = nlam,
  #   eta = NULL, neta = neta,
  #   nfolds = K, foldid = NULL,
  #   intercept = intercept,
  #   standardize = scaling
  # )
  # end.time = Sys.time()
  # slrhsc2_ngmstop.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # slrhsc2_ngmstop.eta.min.idx = slrhsc2_ngmstop$min.idx[2]
  # slrhsc2_ngmstop.lam.min.idx = slrhsc2_ngmstop$min.idx[1]
  # slrhsc2_ngmstop.a0 = slrhsc2_ngmstop$theta0[[slrhsc2_ngmstop.eta.min.idx]][slrhsc2_ngmstop.lam.min.idx]
  # slrhsc2_ngmstop.thetahat = slrhsc2_ngmstop$theta[[slrhsc2_ngmstop.eta.min.idx]][, slrhsc2_ngmstop.lam.min.idx]
  # slrhsc2_ngmstop.SBP = slrhsc2_ngmstop$sbp_thresh[[slrhsc2_ngmstop.eta.min.idx]]
  # slrhsc2_ngmstop.betahat.nonzero = getBeta(slrhsc2_ngmstop.thetahat, sbp = slrhsc2_ngmstop.SBP)
  # slrhsc2_ngmstop.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  # rownames(slrhsc2_ngmstop.betahat) = names(beta)
  # slrhsc2_ngmstop.betahat[slrhsc2_ngmstop$meets_threshold[[slrhsc2_ngmstop.eta.min.idx]], ] =
  #   as.numeric(slrhsc2_ngmstop.betahat.nonzero)
  # 
  # # compute metrics on the selected model #
  # slrhsc2_ngmstop.metrics = getMetricsBalanceReg(
  #   y.train = Y, y.test = Y.test,
  #   ilrX.train = computeBalances(
  #     X[, slrhsc2_ngmstop$meets_threshold[[slrhsc2_ngmstop.eta.min.idx]], drop = FALSE],
  #     sbp = slrhsc2_ngmstop.SBP),
  #   ilrX.test = computeBalances(
  #     X.test[, slrhsc2_ngmstop$meets_threshold[[slrhsc2_ngmstop.eta.min.idx]], drop = FALSE],
  #     sbp = slrhsc2_ngmstop.SBP),
  #   n.train = n, n.test = n,
  #   thetahat0 = slrhsc2_ngmstop.a0, thetahat = slrhsc2_ngmstop.thetahat,
  #   betahat = slrhsc2_ngmstop.betahat,
  #   sbp = slrhsc2_ngmstop.SBP,
  #   true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  # 
  # # # plot the tree given by slr-hsc, indicating significant covariates
  # # slrhsc2_ngmstop_leaf_types = rep("covariate", nrow(slrhsc2_ngmstop.SBP))
  # # slrhsc2_ngmstop_balance_types = rep("balance", ncol(slrhsc2_ngmstop.SBP))
  # # slrhsc2_ngmstop_nodes_types = data.frame(
  # #   name = c(colnames(slrhsc2_ngmstop.SBP), rownames(slrhsc2_ngmstop.SBP)),
  # #   type = c(slrhsc2_ngmstop_balance_types, slrhsc2_ngmstop_leaf_types)
  # # )
  # # plotSBP(slrhsc2_ngmstop.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_ngmstop_nodes_types)
  # # # fields::image.plot(slrSimMat)
  # 
  # saveRDS(c(
  #   slrhsc2_ngmstop.metrics,
  #   "betaSparsity" = bspars,
  #   "time" = slrhsc2_ngmstop.timing
  # ),
  # paste0(output_dir, "/metrics", "/slr_hsc_eta_ngmstop_metrics", file.end))

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
  pr = cvILR(y = Y, X = X, sbp = pr_SBP, nlam = nlam,
             nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  pr.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

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
    "time" = pr.timing
  ),
  paste0(output_dir, "/metrics", "/propr_metrics", file.end))
  
  ##############################################################################
  # oracle method (a balance regression method)
  #   -- hierarchical clustering
  ##############################################################################
  start.time = Sys.time()
  or_SBP = sbp_Q
  or.ilrX = computeBalances(X, sbp = or_SBP)
  oracle = lm(Y ~ or.ilrX)
  end.time = Sys.time()
  or.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  or.coefs = coefficients(oracle)
  or.a0 = or.coefs[1]
  or.thetahat = or.coefs[2]
  or.betahat = getBeta(or.thetahat, sbp = or_SBP)

  # compute metrics on the selected model #
  or.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = computeBalances(X, sbp = or_SBP),
    ilrX.test = computeBalances(X.test, sbp = or_SBP),
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
    "time" = or.timing
  ),
  paste0(output_dir, "/metrics", "/oracle_metrics", file.end))
  
  ##############################################################################
  # selbal method (a balance regression method)
  ##############################################################################
  library(selbal) # masks stats::cor()

  start.time = Sys.time()
  X.slbl = X
  rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
  colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
  slbl = selbal.cv(x = X.slbl, y = as.vector(Y), n.fold = K)
  end.time = Sys.time()
  slbl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # U (transformation) matrix
  U.slbl = rep(0, p)
  names(U.slbl) = colnames(X.slbl)
  pba.pos = unlist(subset(
    slbl$global.balance, subset = Group == "NUM", select = Taxa))
  num.pos = length(pba.pos)
  pba.neg = unlist(subset(
    slbl$global.balance, subset = Group == "DEN", select = Taxa))
  num.neg = length(pba.neg)
  U.slbl[pba.pos] = 1 / num.pos
  U.slbl[pba.neg] = -1 / num.neg
  norm.const = sqrt((num.pos * num.neg) / (num.pos + num.neg))
  U.slbl = norm.const * U.slbl
  # check: these are equal
  # lm(as.vector(Y) ~ log(X) %*% as.matrix(U.slbl))
  # slbl$glm
  slbl.thetahat = coefficients(slbl$glm)[2]
  slbl.betahat = U.slbl %*% as.matrix(slbl.thetahat)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slbl.Yhat.train = predict.glm(
    slbl$glm, newdata = data.frame(X.slbl), type = "response")
  slbl.PE.train = crossprod(Y - slbl.Yhat.train) / n
  # get prediction error on test set
  X.slbl.test = X.test
  rownames(X.slbl.test) = paste("Sample", 1:nrow(X.slbl.test), sep = "_")
  colnames(X.slbl.test) = paste("V", 1:ncol(X.slbl.test), sep = "")
  slbl.Yhat.test = predict.glm(
    slbl$glm, newdata = data.frame(X.slbl.test), type = "response")
  slbl.PE.test = crossprod(Y.test - slbl.Yhat.test) / n
  # estimation accuracy
  slbl.EA = getEstimationAccuracy(true.beta = beta, betahat = slbl.betahat)
  slbl.EA.active = getEstimationAccuracy(
    true.beta = beta[non0.beta], betahat = slbl.betahat[non0.beta])
  slbl.EA.inactive = getEstimationAccuracy(
    true.beta = beta[is0.beta], betahat = slbl.betahat[is0.beta])
  # selection accuracy
  slbl.non0.betahat = abs(slbl.betahat) > 10e-8
  slbl.SA = getSelectionAccuracy(
    is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    non0.betahat = slbl.non0.betahat)

  slbl.metrics = c(
    "PEtr" = slbl.PE.train,
    "PEte" = slbl.PE.test,
    "EA1" = slbl.EA$EA1,
    "EA2" = slbl.EA$EA2,
    "EAInfty" = slbl.EA$EAInfty,
    "EA1Active" = slbl.EA.active$EA1,
    "EA2Active" = slbl.EA.active$EA2,
    "EAInftyActive" = slbl.EA.active$EAInfty,
    "EA1Inactive" = slbl.EA.inactive$EA1,
    "EA2Inactive" = slbl.EA.inactive$EA2,
    "EAInftyInactive" = slbl.EA.inactive$EAInfty,
    "FP" = slbl.SA$FP,
    "FN" = slbl.SA$FN,
    "TPR" = slbl.SA$TPR,
    "precision" = slbl.SA$precision,
    "Fscore" = slbl.SA$Fscore
  )
  saveRDS(c(
    slbl.metrics,
    "betaSparsity" = bspars,
    "time" = slbl.timing
  ),
  paste0(output_dir, "/metrics", "/selbal_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


