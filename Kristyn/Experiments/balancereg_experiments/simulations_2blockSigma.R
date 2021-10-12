# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Note: Here, we simulate from a 2-block-diagonal Sigma and intend to:
#   - have one active contrast corresponding to the two blocks, so that 
#       slr correctly specifies the log-ratio formed by the 2 blocks
# Date: 10/11/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_experiments/outputs"

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
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  
  library(balance) # for sbp.fromHclust()
  library(selbal)
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  
  # Settings to toggle with
  sigma.settings = "2blockSigma"
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "dense" # block, dense
  # "block2" => choose j corresp. to block of correlated variables in block 2
  # "dense" => j = 1
  values.theta = 1
  linkage = "average"
  tol = 1e-4
  nlam = 200
  intercept = TRUE
  K = 10
  n = 100
  p = 200
  cor_ij = 0.2 # 0.2, 0.5
  scaling = TRUE
  
  # Population parameters
  sigma_eps = 0.1 # 0.1, 0.5
  num.blocks = 2
  SigmaWblock = matrix(cor_ij, p / num.blocks, p / num.blocks)
  for(i in 1:nrow(SigmaWblock)) SigmaWblock[i, i] = 1
  SigmaW = as.matrix(bdiag(SigmaWblock, SigmaWblock))
  SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
  U = getU(btree = SigmaWtree) # transformation matrix
  # plot(SigmaWtree)
  # sbp.fromHclust(SigmaWtree)
  
  # theta settings
  if(theta.settings == "dense"){
    indices.theta = 1
  } else if(theta.settings == "block2"){
    SBP = sbp.fromHclust(SigmaWtree)
    # for each column (contrast), find which variables are included (1 or -1)
    contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
    # get the contrasts with length p / 2 -- there are 2 of them
    indices.theta = unname(which(sapply(
      contrast.vars, length) == p / num.blocks)[1])
    ###
  } else{
    stop("invalid theta.settings")
  }
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
  beta = as.vector(getBeta(theta, U = U))
  names(beta) <- paste0('s', 1:p)
  non0.beta = (beta != 0)
  slr.non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
  # Population parameters, continued
  muW = c(rep(log(p), 5), rep(0, p - 5))
  
  file.end = paste0(
    "_dim", n, "x", p, 
    "_", sigma.settings,
    "_", theta.settings, 
    "_noise", sigma_eps,
    "_cor", cor_ij, 
    "_int", intercept,
    "_scale", scaling,
    "_K", K,
    "_seed", rng.seed,
    ".rds")
  
  ##############################################################################
  # simulate data
  logW <- mvrnorm(n = n * 2, mu = muW, Sigma = SigmaW) 
  W <- exp(logW)
  rownames(W) <- colnames(W) <- names(beta)
  XAll <- sweep(W, 1, rowSums(W), FUN='/')
  ilrXAll = computeBalances(XAll, U = U)
  
  # generate Y
  yAll = ilrXAll %*% theta + rnorm(n) * sigma_eps
  
  # subset out training and test sets
  X = XAll[1:n, ]
  X.test = XAll[-(1:n), ]
  Y <- yAll[1:n, , drop = TRUE]
  Y.test <- yAll[-(1:n), , drop = TRUE]
  
  ##############################################################################
  # supervised log-ratios
  ##############################################################################
  
  if(!file.exists(paste0(output_dir, "/models", "/slr_model", b, file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/slr_timing", b, file.end))){
    slr.model.already.existed = FALSE
    # apply supervised log-ratios, using CV to select lambda
    start.time = Sys.time()
    slr = cvSLR(
      y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
      rho.type = rho.type, linkage = linkage, standardize = scaling)
    end.time = Sys.time()
    saveRDS(slr, paste0(output_dir, "/models", "/slr_model", b, file.end))
    
    # timing metric
    slr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      slr.timing, 
      paste0(output_dir, "/timing", "/slr_timing", b, file.end))
  } else{
    slr.model.already.existed = TRUE
    slr = readRDS(paste0(output_dir, "/models", "/slr_model", b, file.end))
    slr.timing = readRDS(paste0(
      output_dir, "/timing", "/slr_timing", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/metrics", "/slr_metrics", b, file.end)) | 
     slr.model.already.existed == FALSE){
    
    # binary tree
    slr.btree = slr$btree
    
    # choose lambda
    slr.lam.min.idx = which.min(slr$cvm)
    slr.lam.min = slr$lambda[slr.lam.min.idx]
    slr.a0 = slr$int[slr.lam.min.idx]
    slr.thetahat = slr$bet[, slr.lam.min.idx]
    slr.Uhat = getU(btree = slr.btree)
    slr.betahat = getBeta(slr.thetahat, U = slr.Uhat)
    
    # evaluate model #
    
    # 1. prediction error #
    # 1a. on training set #
    slr.PE.train = getMSEyhat(
      Y, n, slr.a0, slr.thetahat, computeBalances(X, slr.btree))
    # 1b. on test set #
    slr.PE.test = getMSEyhat(
      Y.test, n, slr.a0, slr.thetahat, computeBalances(X.test, slr.btree))
    
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    slr.EA = getEstimationAccuracy(beta, slr.betahat)
    # 2b. estimation accuracy for active set
    slr.EA.active = getEstimationAccuracy(beta[non0.beta], slr.betahat[non0.beta])
    # 2c. estimation accuracy for inactive set
    slr.EA.inactive = getEstimationAccuracy(beta[is0.beta], slr.betahat[is0.beta])
    
    # 3. selection accuracy #
    # 3a. selection of beta #
    ### using SBP matrix
    slr.SBP = sbp.fromHclust(slr.btree)
    slr.non0.thetahat = (slr.thetahat != 0)
    slr.sel.cols.SBP = slr.SBP[, slr.non0.thetahat, drop = FALSE]
    slr.non0.betahat = apply(slr.sel.cols.SBP, 1, function(row) any(row != 0))
    slr.SA = getSelectionAccuracy(is0.beta, non0.beta, slr.non0.betahat)
    
    saveRDS(c(
      "PEtr" = slr.PE.train, 
      "PEte" = slr.PE.test, 
      "EA1" = slr.EA$EA1, 
      "EA2" = slr.EA$EA2, 
      "EAInfty" = slr.EA$EAInfty, 
      "EA1Active" = slr.EA.active$EA1, 
      "EA2Active" = slr.EA.active$EA2, 
      "EAInftyActive" = slr.EA.active$EAInfty, 
      "EA1Inactive" = slr.EA.inactive$EA1, 
      "EA2Inactive" = slr.EA.inactive$EA2, 
      "EAInftyInactive" = slr.EA.inactive$EAInfty, 
      "FP" = slr.SA$FP, 
      "FN" = slr.SA$FN, 
      "TPR" = slr.SA$TPR, 
      "precision" = slr.SA$precision, 
      "Fscore" = slr.SA$Fscore,
      "timing" = slr.timing,
      "betaSparsity" = bspars
    ), 
    paste0(output_dir, "/metrics", "/slr_metrics", b, file.end))
  } else{
    slr.btree = slr$btree
    slr.SBP = sbp.fromHclust(slr.btree)
  }
  
  if(!file.exists(paste0(output_dir, "/roccurves", "/slr_roc", b, file.end)) | 
     slr.model.already.existed == FALSE){
    # roc
    slr.roc <- apply(slr$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, slr.SBP))
    
    saveRDS(slr.roc, paste0(output_dir, "/roccurves", "/slr_roc", b, file.end))
  }
  
  ##############################################################################
  # compositional lasso
  ##############################################################################
  
  if(!file.exists(paste0(output_dir, "/models", "/classo_model", b, file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/classo_timing", b, file.end))){
    cl.model.already.existed = FALSE
    # apply compositional lasso, using CV to select lambda
    start.time = Sys.time()
    classo = cv.func(
      method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
      nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
    end.time = Sys.time()
    saveRDS(classo, paste0(output_dir, "/models", "/classo_model", b, file.end))
    
    # timing metric
    cl.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      cl.timing, 
      paste0(output_dir, "/timing", "/classo_timing", b, file.end))
  } else{
    cl.model.already.existed = TRUE
    classo = readRDS(paste0(output_dir, "/models", "/classo_model", b, file.end))
    cl.timing = readRDS(paste0(
      output_dir, "/timing", "/classo_timing", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/metrics", "/classo_metrics", b, file.end)) | 
     cl.model.already.existed == FALSE){
    
    # choose lambda
    cl.lam.min.idx = which.min(classo$cvm)
    cl.lam.min = classo$lambda[cl.lam.min.idx]
    cl.a0 = classo$int[cl.lam.min.idx]
    cl.betahat = classo$bet[, cl.lam.min.idx]
    
    # evaluate model #
    
    # 1. prediction error #
    # 1a. on training set #
    cl.PE.train = getMSEyhat(Y, n, cl.a0, cl.betahat, log(X))
    # 1b. on test set #
    cl.PE.test = getMSEyhat(Y.test, n, cl.a0, cl.betahat, log(X.test))
    
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    cl.EA = getEstimationAccuracy(beta, cl.betahat)
    # 2b. estimation accuracy for active set
    cl.EA.active = getEstimationAccuracy(beta[non0.beta], cl.betahat[non0.beta])
    # 2c. estimation accuracy for inactive set
    cl.EA.inactive = getEstimationAccuracy(beta[is0.beta], cl.betahat[is0.beta])
    
    # 3. selection accuracy (i.t.o. beta) #
    cl.non0.betahat = abs(cl.betahat) > 10e-8
    cl.SA = getSelectionAccuracy(is0.beta, non0.beta, cl.non0.betahat)
    
    saveRDS(c(
      "PEtr" = cl.PE.train, 
      "PEte" = cl.PE.test, 
      "EA1" = cl.EA$EA1, 
      "EA2" = cl.EA$EA2, 
      "EAInfty" = cl.EA$EAInfty, 
      "EA1Active" = cl.EA.active$EA1, 
      "EA2Active" = cl.EA.active$EA2, 
      "EAInftyActive" = cl.EA.active$EAInfty, 
      "EA1Inactive" = cl.EA.inactive$EA1, 
      "EA2Inactive" = cl.EA.inactive$EA2, 
      "EAInftyInactive" = cl.EA.inactive$EAInfty, 
      "FP" = cl.SA$FP, 
      "FN" = cl.SA$FN, 
      "TPR" = cl.SA$TPR, 
      "precision" = cl.SA$precision, 
      "Fscore" = cl.SA$Fscore,
      "timing" = cl.timing,
      "betaSparsity" = bspars
    ), 
    paste0(output_dir, "/metrics", "/classo_metrics", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/roccurves", "/classo_roc", b, file.end)) | 
     cl.model.already.existed == FALSE){
    # roc
    cl.roc <- apply(classo$bet, 2, function(a) 
      roc.for.coef(a, beta))
    
    saveRDS(cl.roc, paste0(output_dir, "/roccurves", "/classo_roc", b, file.end))
  }
  
  ##############################################################################
  # oracle method
  ##############################################################################
  
  if(!file.exists(paste0(output_dir, "/models", "/oracle_model", b, file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/oracle_timing", b, file.end))){
    or.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    start.time = Sys.time()
    oracle = cvILR(y = Y, X = X, btree = SigmaWtree, U = U, nlam = nlam, 
                   nfolds = K, intercept = intercept, standardize = scaling)
    end.time = Sys.time()
    saveRDS(oracle, paste0(output_dir, "/models", "/oracle_model", b, file.end))
    
    # timing metric
    or.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      or.timing, 
      paste0(output_dir, "/timing", "/oracle_timing", b, file.end))
  } else{
    or.model.already.existed = TRUE
    oracle = readRDS(paste0(output_dir, "/models", "/oracle_model", b, file.end))
    or.timing = readRDS(paste0(
      output_dir, "/timing", "/oracle_timing", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/metrics", "/oracle_metrics", b, file.end)) | 
     or.model.already.existed == FALSE){
    
    # binary tree
    or.btree = SigmaWtree
    
    # timing metric
    or.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    
    # choose lambda
    or.lam.min.idx = which.min(oracle$cvm)
    or.lam.min = oracle$lambda[or.lam.min.idx]
    or.a0 = oracle$int[or.lam.min.idx]
    or.thetahat = oracle$bet[, or.lam.min.idx]
    or.Uhat = getU(btree = or.btree)
    or.betahat = getBeta(or.thetahat, U = or.Uhat)
    
    # evaluate model #
    
    # 1. prediction error #
    # 1a. on training set #
    or.PE.train = getMSEyhat(
      Y, n, or.a0, or.thetahat, computeBalances(X, or.btree))
    # 1b. on test set #
    or.PE.test = getMSEyhat(
      Y.test, n, or.a0, or.thetahat, computeBalances(X.test, or.btree))
    
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    or.EA = getEstimationAccuracy(beta, or.betahat)
    # 2b. estimation accuracy for active set
    or.EA.active = getEstimationAccuracy(beta[non0.beta], or.betahat[non0.beta])
    # 2c. estimation accuracy for inactive set
    or.EA.inactive = getEstimationAccuracy(beta[is0.beta], or.betahat[is0.beta])
    
    # 3. selection accuracy #
    # 3a. selection of beta #
    ### using SBP matrix
    or.SBP = sbp.fromHclust(or.btree)
    row.names(or.SBP) = colnames(W)
    or.non0.thetahat = (or.thetahat != 0)
    or.sel.cols.SBP = or.SBP[, or.non0.thetahat, drop = FALSE]
    or.non0.betahat = apply(or.sel.cols.SBP, 1, function(row) any(row != 0))
    or.SA = getSelectionAccuracy(is0.beta, non0.beta, or.non0.betahat)
    
    saveRDS(c(
      "PEtr" = or.PE.train, 
      "PEte" = or.PE.test, 
      "EA1" = or.EA$EA1, 
      "EA2" = or.EA$EA2, 
      "EAInfty" = or.EA$EAInfty, 
      "EA1Active" = or.EA.active$EA1, 
      "EA2Active" = or.EA.active$EA2, 
      "EAInftyActive" = or.EA.active$EAInfty, 
      "EA1Inactive" = or.EA.inactive$EA1, 
      "EA2Inactive" = or.EA.inactive$EA2, 
      "EAInftyInactive" = or.EA.inactive$EAInfty, 
      "FP" = or.SA$FP, 
      "FN" = or.SA$FN, 
      "TPR" = or.SA$TPR, 
      "precision" = or.SA$precision, 
      "Fscore" = or.SA$Fscore,
      "timing" = or.timing,
      "betaSparsity" = bspars
    ), 
    paste0(output_dir, "/metrics", "/oracle_metrics", b, file.end))
  } else{
    or.btree = SigmaWtree
    or.SBP = sbp.fromHclust(or.btree)
    row.names(or.SBP) = colnames(W)
  }
  
  if(!file.exists(paste0(output_dir, "/roccurves", "/oracle_roc", b, file.end)) | 
     or.model.already.existed == FALSE){
    # roc
    or.roc <- apply(oracle$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, or.SBP))
    
    saveRDS(or.roc, paste0(output_dir, "/roccurves", "/oracle_roc", b, file.end))
  }
  
  ##############################################################################
  # propr method
  ##############################################################################
  
  if(!file.exists(paste0(output_dir, "/models", "/propr_model", b, file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/propr_timing", b, file.end))){
    pr.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    start.time = Sys.time()
    pr <- propr(X, metric = "phs")
    pr.tree = hclust(as.dist(pr@matrix),method = linkage)
    pr = cvILR(y = Y, X = X, btree = pr.tree, nlam = nlam, 
               nfolds = K, intercept = intercept, standardize = scaling)
    end.time = Sys.time()
    saveRDS(pr, paste0(output_dir, "/models", "/propr_model", b, file.end))
    
    # timing metric
    pr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      pr.timing, 
      paste0(output_dir, "/timing", "/propr_timing", b, file.end))
  } else{
    pr.model.already.existed = TRUE
    pr = readRDS(paste0(output_dir, "/models", "/propr_model", b, file.end))
    pr.timing = readRDS(paste0(
      output_dir, "/timing", "/propr_timing", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/metrics", "/propr_metrics", b, file.end)) | 
     pr.model.already.existed == FALSE){
    
    # binary tree
    pr.btree = pr$btree
    
    # timing metric
    pr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    
    # choose lambda
    pr.lam.min.idx = which.min(pr$cvm)
    pr.lam.min = pr$lambda[pr.lam.min.idx]
    pr.a0 = pr$int[pr.lam.min.idx]
    pr.thetahat = pr$bet[, pr.lam.min.idx]
    pr.Uhat = getU(btree = pr.btree)
    pr.betahat = getBeta(pr.thetahat, U = pr.Uhat)
    
    # evaluate model #
    
    # 1. prediction error #
    # 1a. on training set #
    pr.PE.train = getMSEyhat(
      Y, n, pr.a0, pr.thetahat, computeBalances(X, pr.btree))
    # 1b. on test set #
    pr.PE.test = getMSEyhat(
      Y.test, n, pr.a0, pr.thetahat, computeBalances(X.test, pr.btree))
    
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    pr.EA = getEstimationAccuracy(beta, pr.betahat)
    # 2b. estimation accuracy for active set
    pr.EA.active = getEstimationAccuracy(beta[non0.beta], pr.betahat[non0.beta])
    # 2c. estimation accuracy for inactive set
    pr.EA.inactive = getEstimationAccuracy(beta[is0.beta], pr.betahat[is0.beta])
    
    # 3. selection accuracy #
    # 3a. selection of beta #
    ### using SBP matrix
    pr.SBP = sbp.fromHclust(pr.btree)
    pr.non0.thetahat = (pr.thetahat != 0)
    pr.sel.cols.SBP = pr.SBP[, pr.non0.thetahat, drop = FALSE]
    pr.non0.betahat = apply(pr.sel.cols.SBP, 1, function(row) any(row != 0))
    pr.SA = getSelectionAccuracy(is0.beta, non0.beta, pr.non0.betahat)
    
    saveRDS(c(
      "PEtr" = pr.PE.train, 
      "PEte" = pr.PE.test, 
      "EA1" = pr.EA$EA1, 
      "EA2" = pr.EA$EA2, 
      "EAInfty" = pr.EA$EAInfty, 
      "EA1Active" = pr.EA.active$EA1, 
      "EA2Active" = pr.EA.active$EA2, 
      "EAInftyActive" = pr.EA.active$EAInfty, 
      "EA1Inactive" = pr.EA.inactive$EA1, 
      "EA2Inactive" = pr.EA.inactive$EA2, 
      "EAInftyInactive" = pr.EA.inactive$EAInfty, 
      "FP" = pr.SA$FP, 
      "FN" = pr.SA$FN, 
      "TPR" = pr.SA$TPR, 
      "precision" = pr.SA$precision, 
      "Fscore" = pr.SA$Fscore,
      "timing" = pr.timing,
      "betaSparsity" = bspars
    ), 
    paste0(output_dir, "/metrics", "/propr_metrics", b, file.end))
  } else{
    pr.btree = SigmaWtree
    pr.SBP = sbp.fromHclust(pr.btree)
  }
  
  if(!file.exists(paste0(output_dir, "/roccurves", "/propr_roc", b, file.end)) | 
     pr.model.already.existed == FALSE){
    # roc
    pr.roc <- apply(pr$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, pr.SBP))
    
    saveRDS(or.roc, paste0(output_dir, "/roccurves", "/propr_roc", b, file.end))
  }
}


