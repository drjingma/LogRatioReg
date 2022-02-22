# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 2/9/2022

################################################################################
# libraries and settings

# (put in the loop)
output_dir = "Kristyn/Experiments/current_experiments/outputs"

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
# single_balance = FALSE # (put in the loop)
# theta_overlapping_balance = TRUE # (put in the loop)

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
# theta.value = 1
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
# 
# # define theta (directly related to sparsity)
# if(single_balance){
#   registerDoRNG(456)
#   theta = rep(0, p - 1)
#   theta[5] = theta.value
# } else{
#   if(theta_overlapping_balance){
#     registerDoRNG(123) # for tree with an overlapped balance
#   } else{
#     registerDoRNG(456) # for tree with no overlapping balances
#   }
#   # randomly choose 5% of balances/internal nodes that at most 5 active variables
#   #   i.e. excluse balances/internal nodes with more than 5 active variables
#   contrast_vars = apply(sbp_Q, 2, FUN = function(col) which(col != 0))
#   contrast_lens = sapply(contrast_vars, length)
#   viable_contrasts = colnames(sbp_Q)[contrast_lens <= 5]
#   num_select_viable_contrasts = signif(length(viable_contrasts) * 0.05, 0)
#   # View(sbp_Q[, viable_contrasts])
#   theta = rep(0, p - 1)
#   names(theta) = colnames(sbp_Q)
#   selected_viable_contrasts = sample(
#     x = viable_contrasts, size = num_select_viable_contrasts,
#     replace = FALSE)
#   theta[selected_viable_contrasts] = theta.value
# }
# non0.theta = (theta != 0)
# which(non0.theta)
# 
# # define beta
# beta = as.vector(getBetaFromTheta(theta = theta, sbp = sbp_Q))
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
# if(single_balance){
#   saveRDS(
#     theta,
#     file = paste0(
#       output_dir,
#       "/Qmatrix_theta",
#       "_singlebalance", single_balance,
#       ".rds")
#   )
# } else{
#   saveRDS(
#     theta,
#     file = paste0(
#       output_dir,
#       "/Qmatrix_theta",
#       "_balanceoverlap", theta_overlapping_balance,
#       ".rds")
#   )
# }


################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSims
) %dorng% {
  # rm(list=ls())
  
  output_dir = "Kristyn/Experiments/balancereg_hsclust_experiments/outputs"
  single_balance = FALSE
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
  source("Kristyn/Functions/slrnew.R")
  
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
  theta.value = 1
  n = 100
  p = nrow(Q)
  intercept = TRUE
  scaling = TRUE
  K = 10
  linkage = "average"
  tol = 1e-4
  nlam = 100
  neta = p
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
  if(single_balance){
    theta = readRDS(paste0(
        output_dir,
        "/Qmatrix_theta",
        "_singlebalance", single_balance,
        ".rds")
    )
  } else{
    theta = readRDS(paste0(
        output_dir,
        "/Qmatrix_theta",
        "_balanceoverlap", theta_overlapping_balance,
        ".rds")
    )
  }
  
  # define beta
  beta = as.vector(getBetaFromTheta(theta = theta, sbp = sbp_Q))
  names(beta) <- rownames(sbp_Q)
  non0.beta = (beta != 0)
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
  # saving stuff
  if(single_balance){
    file.end = paste0(
      "_", sigma.settings,
      "_singlebalance", single_balance,
      "_dim", n, "x", p, 
      "_noise", sigma_eps, 
      "_sim", b,
      ".rds")
  } else{
    file.end = paste0(
      "_", sigma.settings,
      "_balanceoverlap", theta_overlapping_balance,
      "_dim", n, "x", p, 
      "_noise", sigma_eps, 
      "_sim", b,
      ".rds")
  }
  
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
  
  # ##############################################################################
  # # supervised log-ratios (a balance regression method)
  # #   -- hierarchical spectral clustering + thresholding with lasso
  # ##############################################################################
  # start.time = Sys.time()
  # # apply hierarchical spectral clustering to the SLR similarity matrix
  # slrSimMat = getSlrMatrix(
  #   y = Y, X = X, type = "similarity")
  # slrhsc_btree = HSClust(
  #   W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  # slrhsc_SBP = sbp.fromHSClust(
  #   levels_matrix = slrhsc_btree$allLevels, row_names = names(beta))
  # # apply supervised log-ratios, using CV to select threshold and also lambda
  # slrhsc2 = cvBMLassoThresh(
  #   y = Y, X = X,
  #   W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
  #   hsc_method = "kmeans", # "shimalik", "kmeans"
  #   force_levelMax = TRUE,
  #   sbp = slrhsc_SBP,
  #   lambda = NULL, nlam = nlam,
  #   eta = NULL, neta = neta,
  #   nfolds = K, foldid = NULL,
  #   intercept = intercept,
  #   standardize = scaling
  # )
  # end.time = Sys.time()
  # slrhsc2.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
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
  #   true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  # 
  # # # plot the tree given by slr-hsc, indicating significant covariates
  # # slrhsc2_leaf_types = rep("covariate", nrow(slrhsc2.SBP))
  # # slrhsc2_balance_types = rep("balance", ncol(slrhsc2.SBP))
  # # slrhsc2_nodes_types = data.frame(
  # #   name = c(colnames(slrhsc2.SBP), rownames(slrhsc2.SBP)),
  # #   type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
  # # )
  # # plotSBP(slrhsc2.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_nodes_types)
  # # # fields::image.plot(slrSimMat)
  # 
  # saveRDS(c(
  #   slrhsc2.metrics,
  #   "betaSparsity" = bspars,
  #   "Rsq" = Rsq,
  #   "time" = slrhsc2.timing
  # ),
  # paste0(output_dir, "/metrics", "/slr_hsc_thresh_lasso_metrics", file.end))
  # 
  # ##############################################################################
  # # compositional lasso (a linear log contrast method)
  # ##############################################################################
  # start.time = Sys.time()
  # classo = cv.func(
  #   method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
  #   nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  # end.time = Sys.time()
  # cl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # cl.lam.min.idx = which.min(classo$cvm)
  # cl.a0 = classo$int[cl.lam.min.idx]
  # cl.betahat = classo$bet[, cl.lam.min.idx]
  # 
  # # compute metrics on the selected model #
  # cl.metrics = getMetricsLLC(
  #   y.train = Y, y.test = Y.test,
  #   logX.train = log(X),
  #   logX.test = log(X.test),
  #   n.train = n, n.test = n,
  #   betahat0 = cl.a0, betahat = cl.betahat,
  #   true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  # 
  # saveRDS(c(
  #   cl.metrics,
  #   "betaSparsity" = bspars,
  #   "Rsq" = Rsq,
  #   "time" = cl.timing
  # ),
  # paste0(output_dir, "/metrics", "/classo_metrics", file.end))
  # 
  # ##############################################################################
  # # new slr method (a balance regression method)
  # #   -- hierarchical spectral clustering
  # ##############################################################################
  # start.time = Sys.time()
  # slrnew = slr(x = X, y = Y)
  # slrnew_activevars = names(slrnew$index)
  # slrnew_SBP = matrix(slrnew$index)
  # rownames(slrnew_SBP) = slrnew_activevars
  # end.time = Sys.time()
  # slrnew.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # slrnew.coefs = coefficients(slrnew$model)
  # slrnew.a0 = slrnew.coefs[1]
  # slrnew.thetahat = slrnew.coefs[-1]
  # 
  # slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
  # slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  # rownames(slrnew.betahat) = names(beta)
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
  #   true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  # 
  # saveRDS(c(
  #   slrnew.metrics,
  #   "betaSparsity" = bspars,
  #   "Rsq" = Rsq,
  #   "time" = slrnew.timing
  # ),
  # paste0(output_dir, "/metrics", "/slr_new_metrics", file.end))
  
  ##############################################################################
  # new slr method (a balance regression method)
  #   -- hierarchical spectral clustering (without rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slrnew = slr(x = X, y = Y, rank1approx = FALSE)
  slrnew_activevars = names(slrnew$index)
  slrnew_SBP = matrix(slrnew$index)
  rownames(slrnew_SBP) = slrnew_activevars
  end.time = Sys.time()
  slrnew.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
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
    true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
  
  saveRDS(c(
    slrnew.metrics,
    "betaSparsity" = bspars,
    "Rsq" = Rsq,
    "time" = slrnew.timing
  ),
  paste0(output_dir, "/metrics", "/slr_new_noapprox_metrics", file.end))
  
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
  # 
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


