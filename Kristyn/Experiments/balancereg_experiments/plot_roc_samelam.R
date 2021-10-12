rm(list=ls())
# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Date: 10/11/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_experiments/outputs"

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100
rng.seed = 123

# Settings to toggle with
sigma.settings = "lin14Sigma" # 2blockSigma, 4blockSigma, 10blockSigma, lin14Sigma
rho.type = "square" # 1 = "absolute value", 2 = "square"
theta.settings = "dense" # "dense", "sparse", "both", "multsparse"
# if "2blockSigma" then "dense"
# if "4blockSigma", then "2blocks"
# if "10blockSigma", then "pairperblock" or "1blockpair4halves"
# if "lin14Sigma" then "dense" or "multsparse"
mu.settings = "" # matchbeta
linkage = "average"
tol = 1e-4
nlam = 200
intercept = TRUE
K = 10
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
cor_ij = 0.2 # 0.2, 0.5
scaling = TRUE
sigma_eps = 0.1  # 0.1, 0.5

if(sigma.settings == "lin14Sigma"){
  if(mu.settings == "matchbeta"){
    file.end = paste0( # for old simulations
      "_dim", n, "x", p,
      "_", sigma.settings,
      "_", theta.settings,
      "_", mu.settings,
      "_noise", sigma_eps,
      "_rho", rho,
      "_int", intercept,
      "_scale", scaling,
      "_K", K,
      "_seed", rng.seed,
      ".rds")
  } else{
    file.end = paste0( # for old simulations
      "_dim", n, "x", p,
      "_", sigma.settings,
      "_", theta.settings,
      "_noise", sigma_eps,
      "_rho", rho,
      "_int", intercept,
      "_scale", scaling,
      "_K", K,
      "_seed", rng.seed,
      ".rds")
  }
} else{ # for block-diagonal Sigma, either "2blockSigma" or "4blockSigma"
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
}

has.selbal = FALSE
has.coat = FALSE
has.oracle = TRUE
has.propr = TRUE
if(FALSE){
  has.selbal = TRUE
}
# if(sigma.settings == "10blockSigma"){
#   has.coat = FALSE
# }

################################################################################
# get lambda sequences

# import lambdas
cl.lams.list = list()
slr.lams.list = list()
if(has.selbal) selbal.lams.list = list()
if(has.oracle) or.lams.list = list()
if(has.coat) coat.lams.list = list()
if(has.propr) pr.lams.list = list()
for(i in 1:numSims){
  # classo
  cl.lam.tmp = readRDS(
    paste0(output_dir, "/models", "/classo_model", i, file.end
  ))$lambda
  rownames(cl.lam.tmp) = NULL
  cl.lams.list[[i]] = data.table(cl.lam.tmp)
  # slr
  slr.lam.tmp = readRDS(
    paste0(output_dir, "/models", "/slr_model", i, file.end
    ))$lambda
  rownames(slr.lam.tmp) = NULL
  slr.lams.list[[i]] = data.table(slr.lam.tmp)
  if(has.selbal){
    # selbal
    selbal.lam.tmp = readRDS(
      paste0(output_dir, "/models", "/selbal_model", i, file.end
      ))$lambda
    rownames(selbal.lam.tmp) = NULL
    selbal.lams.list[[i]] = data.table(selbal.lam.tmp)
  }
  if(has.oracle){
    # oracle
    or.lam.tmp = readRDS(
      paste0(output_dir, "/models", "/oracle_model", i, file.end
      ))$lambda
    rownames(or.lam.tmp) = NULL
    or.lams.list[[i]] = data.table(or.lam.tmp)
  }
  if(has.coat){
    # coat
    coat.lam.tmp = readRDS(
      paste0(output_dir, "/models", "/coat_model", i, file.end
      ))$lambda
    rownames(coat.lam.tmp) = NULL
    coat.lams.list[[i]] = data.table(coat.lam.tmp)
  }
  if(has.propr){
    # propr
    pr.lam.tmp = readRDS(
      paste0(output_dir, "/models", "/propr_model", i, file.end
      ))$lambda
    rownames(pr.lam.tmp) = NULL
    pr.lams.list[[i]] = data.table(pr.lam.tmp)
  }
}
cl.lams = as.matrix(do.call(cbind, cl.lams.list))
slr.lams = as.matrix(do.call(cbind, slr.lams.list))
if(has.selbal) selbal.lams = as.matrix(do.call(cbind, selbal.lams.list))
if(has.oracle) or.lams = as.matrix(do.call(cbind, or.lams.list))
if(has.coat) coat.lams = as.matrix(do.call(cbind, coat.lams.list))
if(has.propr) pr.lams = as.matrix(do.call(cbind, pr.lams.list))

# lambda min and max values
cl.bounds = c(min(cl.lams), max(cl.lams))
slr.bounds = c(min(slr.lams), max(slr.lams))
if(has.selbal) selbal.bounds = c(min(has.selbal), max(has.selbal))
if(has.oracle) or.bounds = c(min(has.oracle), max(has.oracle))
if(has.coat) coat.bounds = c(min(coat.lams), max(coat.lams))
if(has.propr) pr.bounds = c(min(pr.lams), max(pr.lams))

# new lambda sequences, one for each method
cl.lambda.seq = exp(seq(max(cl.bounds), min(cl.bounds),length.out = nlam))
slr.lambda.seq = exp(seq(max(slr.bounds), min(slr.bounds),length.out = nlam))
if(has.selbal){
  selbal.lambda.seq = exp(seq(max(selbal.bounds), min(selbal.bounds),length.out = nlam))
}
if(has.oracle){
  or.lambda.seq = exp(seq(max(or.bounds), min(or.bounds),length.out = nlam))
}
if(has.coat){
  coat.lambda.seq = exp(seq(max(coat.bounds), min(coat.bounds),length.out = nlam))
}
if(has.propr){
  pr.lambda.seq = exp(seq(max(pr.bounds), min(pr.bounds),length.out = nlam))
}

################################################################################
# solution path calculations

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
registerDoRNG(rng.seed)

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
  library(propr)
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  source("Kristyn/Functions/supervisedlogratiosalpha.R")
  source("COAT-master/coat.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  
  # Settings to toggle with
  sigma.settings = "10blockSigma"
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "1blockpair4halves"  
  # "pairperblock" => choose j corresp. to one pair of covariates for each block
  # "2blockpairs4halves" => 
  #   2 contrasts corresponding to 2 blocks each (accounts for 4 blocks so far), 
  #   4 contrasts, each corresponding to half (or approx. half) of the variables 
  #     in 4 different blocks (accounts for 8 blocks so far), and 
  #   the other two blocks with inactive variables (i.e. not in any of the 
  #     selected contrasts).
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
  num.blocks = 10
  SigmaWblock = matrix(cor_ij, p / num.blocks, p / num.blocks)
  for(i in 1:nrow(SigmaWblock)) SigmaWblock[i, i] = 1
  SigmaW = as.matrix(bdiag(
    SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, 
    SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock))
  SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
  U = getU(btree = SigmaWtree) # transformation matrix
  # plot(SigmaWtree)
  # sbp.fromHclust(SigmaWtree)
  
  # theta settings
  SBP = sbp.fromHclust(SigmaWtree)
  if(theta.settings == "pairperblock"){
    # for each column (contrast), find which variables are included (1 or -1)
    contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
    # get the contrasts with length p / 2 -- have 2 blocks of correlated vars
    #   not necessary, but may save on unnecessary computation in the next step
    block.contrasts = which(sapply(contrast.vars, length) == 2)
    # pick one such contrast
    if(length(block.contrasts) >= num.blocks){
      block.contrasts.pairs = do.call(rbind, contrast.vars[block.contrasts])
      block.labels = cut(1:p, num.blocks)
      indices.theta = rep(NA, num.blocks)
      for(i in 1:num.blocks){
        pairs.in.block.i = apply(
          block.contrasts.pairs, 1, FUN = function(x) all(
            as.numeric(x) %in% ((i - 1) * (p / num.blocks) + (1:(p / num.blocks)))))
        contrasts.block.i = block.contrasts[pairs.in.block.i]
        indices.theta[i] = sample(contrasts.block.i, 1)
      }
      # SBP[, indices.theta]
    } else{
      stop("there aren't 10 different contrasts corresponding to different pairs in each block!")
    }
  } else if(theta.settings == "1blockpair4halves"){
    # "1blockpair4halves" => 
    #   1 contrast corresponding to 2 blocks (accounts for 2 blocks so far), 
    #   4 contrasts, each corresponding to half (or approx. half) of the vars 
    #     in 4 different blocks (accounts for 8 blocks so far), and 
    #   the other 4 blocks with inactive vars (i.e. not in any of the 
    #     selected contrasts).
    contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
    # get the 1 contrast corresponding to 2 blocks
    block.contrasts.1blockpair = which(
      sapply(contrast.vars, length) == 2 * (p / num.blocks))
    indices.theta1 = unname(block.contrasts.1blockpair)
    # get the 4 contrasts, each corresponding to half (or approx. half) of the 
    #   vars in 4 different blocks
    block.contrasts.halves = which(
      sapply(contrast.vars, length) == 9) # 0.5 * (p / num.blocks))
    indices.theta2 = unname(sample(block.contrasts.halves, 4))
    indices.theta = c(indices.theta1, indices.theta2)
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
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam", "/slr_model", b, file.end))){
    slr.model.already.existed = FALSE
    # apply supervised log-ratios, using CV to select lambda
    slr = cvSLR(
      y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
      rho.type = rho.type, linkage = linkage, standardize = scaling)
    saveRDS(slr, paste0(output_dir, "/roc_samelam", "/slr_model", b, file.end))
  } else{
    slr.model.already.existed = TRUE
    slr = readRDS(paste0(output_dir, "/roc_samelam", "/slr_model", b, file.end))
  }
  
    slr.btree = slr$btree
    slr.SBP = sbp.fromHclust(slr.btree)
  
  if(!file.exists(paste0(output_dir, "/roc_samelam", "/slr_roc", b, file.end)) | 
     slr.model.already.existed == FALSE){
    # roc
    slr.roc <- apply(slr$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, slr.SBP))
    
    saveRDS(slr.roc, paste0(output_dir, "/roc_samelam", "/slr_roc", b, file.end))
  }
  
  ##############################################################################
  # compositional lasso
  ##############################################################################
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam", "/classo_model", b, file.end))){
    cl.model.already.existed = FALSE
    # apply compositional lasso, using CV to select lambda
    classo = cv.func(
      method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
      nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
    saveRDS(classo, paste0(output_dir, "/roc_samelam", "/classo_model", b, file.end))
  } else{
    cl.model.already.existed = TRUE
    classo = readRDS(paste0(output_dir, "/roc_samelam", "/classo_model", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/roc_samelam", "/classo_roc", b, file.end)) | 
     cl.model.already.existed == FALSE){
    # roc
    cl.roc <- apply(classo$bet, 2, function(a) 
      roc.for.coef(a, beta))
    
    saveRDS(cl.roc, paste0(output_dir, "/roc_samelam", "/classo_roc", b, file.end))
  }
  
  ##############################################################################
  # oracle method
  ##############################################################################
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam", "/oracle_model", b, file.end))){
    or.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    oracle = cvILR(y = Y, X = X, btree = SigmaWtree, U = U, nlam = nlam, 
                   nfolds = K, intercept = intercept, standardize = scaling)
    saveRDS(oracle, paste0(output_dir, "/roc_samelam", "/oracle_model", b, file.end))
  } else{
    or.model.already.existed = TRUE
    oracle = readRDS(paste0(output_dir, "/roc_samelam", "/oracle_model", b, file.end))
  }
  
    or.btree = SigmaWtree
    or.SBP = sbp.fromHclust(or.btree)
    row.names(or.SBP) = colnames(W)
  
  if(!file.exists(paste0(output_dir, "/roc_samelam", "/oracle_roc", b, file.end)) | 
     or.model.already.existed == FALSE){
    # roc
    or.roc <- apply(oracle$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, or.SBP))
    
    saveRDS(or.roc, paste0(output_dir, "/roc_samelam", "/oracle_roc", b, file.end))
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
    pr.btree = pr$btree
    pr.SBP = sbp.fromHclust(pr.btree)
  }
  
  if(!file.exists(paste0(output_dir, "/roccurves", "/propr_roc", b, file.end)) | 
     pr.model.already.existed == FALSE){
    # roc
    pr.roc <- apply(pr$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, pr.SBP))
    
    saveRDS(pr.roc, paste0(output_dir, "/roccurves", "/propr_roc", b, file.end))
  }
  
  ##############################################################################
  # supervised log-ratios alpha = 0.5
  ##############################################################################
  
  if(!file.exists(paste0(output_dir, "/models", "/slralpha0.5_model", b, file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/slralpha0.5_timing", b, file.end))){
    slr0.5.model.already.existed = FALSE
    # apply supervised log-ratios, using CV to select lambda
    alpha = 0.5
    start.time = Sys.time()
    slr0.5 = cvSLRalpha(
      y = Y, X = X, nlam = nlam, nfolds = K, alpha = alpha, 
      intercept = intercept, rho.type = rho.type, linkage = linkage, 
      scaling = scaling)
    end.time = Sys.time()
    saveRDS(slr0.5, paste0(output_dir, "/models", "/slralpha0.5_model", b, file.end))
    
    # timing metric
    slr0.5.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      slr0.5.timing, 
      paste0(output_dir, "/timing", "/slralpha0.5_timing", b, file.end))
  } else{
    slr0.5.model.already.existed = TRUE
    slr0.5 = readRDS(paste0(output_dir, "/models", "/slralpha0.5_model", b, file.end))
    slr0.5.timing = readRDS(paste0(
      output_dir, "/timing", "/slralpha0.5_timing", b, file.end))
  }
  
  if(!file.exists(paste0(output_dir, "/metrics", "/slralpha0.5_metrics", b, file.end)) | 
     slr0.5.model.already.existed == FALSE){
    
    # binary tree
    slr0.5.btree = slr0.5$btree
    
    # choose lambda
    slr0.5.lam.min.idx = which.min(slr0.5$cvm)
    slr0.5.lam.min = slr0.5$lambda[slr0.5.lam.min.idx]
    slr0.5.a0 = slr0.5$int[slr0.5.lam.min.idx]
    slr0.5.thetahat = slr0.5$bet[, slr0.5.lam.min.idx]
    slr0.5.Uhat = getU(btree = slr0.5.btree)
    slr0.5.betahat = getBeta(slr0.5.thetahat, U = slr0.5.Uhat)
    
    # evaluate model #
    
    # 1. prediction error #
    # 1a. on training set #
    slr0.5.PE.train = getMSEyhat(
      Y, n, slr0.5.a0, slr0.5.thetahat, computeBalances(X, slr0.5.btree))
    # 1b. on test set #
    slr0.5.PE.test = getMSEyhat(
      Y.test, n, slr0.5.a0, slr0.5.thetahat, computeBalances(X.test, slr0.5.btree))
    
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    slr0.5.EA = getEstimationAccuracy(beta, slr0.5.betahat)
    # 2b. estimation accuracy for active set
    slr0.5.EA.active = getEstimationAccuracy(beta[non0.beta], slr0.5.betahat[non0.beta])
    # 2c. estimation accuracy for inactive set
    slr0.5.EA.inactive = getEstimationAccuracy(beta[is0.beta], slr0.5.betahat[is0.beta])
    
    # 3. selection accuracy #
    # 3a. selection of beta #
    ### using SBP matrix
    slr0.5.SBP = sbp.fromHclust(slr0.5.btree)
    slr0.5.non0.thetahat = (slr0.5.thetahat != 0)
    slr0.5.sel.cols.SBP = slr0.5.SBP[, slr0.5.non0.thetahat, drop = FALSE]
    slr0.5.non0.betahat = apply(slr0.5.sel.cols.SBP, 1, function(row) any(row != 0))
    slr0.5.SA = getSelectionAccuracy(is0.beta, non0.beta, slr0.5.non0.betahat)
    
    saveRDS(c(
      "PEtr" = slr0.5.PE.train, 
      "PEte" = slr0.5.PE.test, 
      "EA1" = slr0.5.EA$EA1, 
      "EA2" = slr0.5.EA$EA2, 
      "EAInfty" = slr0.5.EA$EAInfty, 
      "EA1Active" = slr0.5.EA.active$EA1, 
      "EA2Active" = slr0.5.EA.active$EA2, 
      "EAInftyActive" = slr0.5.EA.active$EAInfty, 
      "EA1Inactive" = slr0.5.EA.inactive$EA1, 
      "EA2Inactive" = slr0.5.EA.inactive$EA2, 
      "EAInftyInactive" = slr0.5.EA.inactive$EAInfty, 
      "FP" = slr0.5.SA$FP, 
      "FN" = slr0.5.SA$FN, 
      "TPR" = slr0.5.SA$TPR, 
      "precision" = slr0.5.SA$precision, 
      "Fscore" = slr0.5.SA$Fscore,
      "timing" = slr0.5.timing,
      "betaSparsity" = bspars
    ), 
    paste0(output_dir, "/metrics", "/slralpha0.5_metrics", b, file.end))
  } else{
    slr0.5.btree = slr0.5$btree
    slr0.5.SBP = sbp.fromHclust(slr0.5.btree)
  }
  
  if(!file.exists(paste0(output_dir, "/roccurves", "/slralpha0.5_roc", b, file.end)) | 
     slr0.5.model.already.existed == FALSE){
    # roc
    slr0.5.roc <- apply(slr0.5$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, slr0.5.SBP))
    
    saveRDS(slr0.5.roc, paste0(output_dir, "/roccurves", "/slralpha0.5_roc", b, file.end))
  }
  
}



################################################################################
# plot roc curves

# import roc curves and organize TPR, S.hat, lambda information
# cl
cl.roc.list = list()
# each row corresponds to a lambda in the lambda sequence (different in ea. sim)
# each column corresponds to a different simulation
cl.TPR.mat = matrix(NA, nlam, numSims) 
cl.S.hat.mat = matrix(NA, nlam, numSims)
cl.TP.mat = matrix(NA, nlam, numSims)
# slr
slr.roc.list = list()
slr.TPR.mat = matrix(NA, nlam, numSims)
slr.S.hat.mat = matrix(NA, nlam, numSims)
slr.TP.mat = matrix(NA, nlam, numSims)
if(has.oracle){
  # oracle
  or.roc.list = list()
  or.TPR.mat = matrix(NA, nlam, numSims)
  or.S.hat.mat = matrix(NA, nlam, numSims)
  or.TP.mat = matrix(NA, nlam, numSims)
}
if(has.coat){
  # coat
  coat.roc.list = list()
  coat.TPR.mat = matrix(NA, nlam, numSims)
  coat.S.hat.mat = matrix(NA, nlam, numSims)
  coat.TP.mat = matrix(NA, nlam, numSims)
}
if(has.propr){
  # coat
  pr.roc.list = list()
  pr.TPR.mat = matrix(NA, nlam, numSims)
  pr.S.hat.mat = matrix(NA, nlam, numSims)
  pr.TP.mat = matrix(NA, nlam, numSims)
}
for(i in 1:numSims){
  # cl
  cl.sim.tmp = readRDS(paste0(
    output_dir, "/roccurves", "/classo_roc", i, file.end
  ))
  cl.roc.list[[i]] = cl.sim.tmp
  cl.TPR.mat[, i] = cl.sim.tmp["tpr", ]
  cl.S.hat.mat[, i] = cl.sim.tmp["S_hat", ]
  cl.TP.mat[, i] = cl.sim.tmp["TP", ]
  # slr
  slr.sim.tmp = readRDS(paste0(
    output_dir, "/roccurves", "/slr_roc", i, file.end
  ))
  slr.roc.list[[i]] = slr.sim.tmp
  slr.TPR.mat[, i] = slr.sim.tmp["tpr", ]
  slr.S.hat.mat[, i] = slr.sim.tmp["S_hat", ]
  slr.TP.mat[, i] = slr.sim.tmp["TP", ]
  if(has.oracle){
    # oracle
    or.sim.tmp = readRDS(paste0(
      output_dir, "/roccurves", "/oracle_roc", i, file.end
    ))
    or.roc.list[[i]] = or.sim.tmp
    or.TPR.mat[, i] = or.sim.tmp["tpr", ]
    or.S.hat.mat[, i] = or.sim.tmp["S_hat", ]
    or.TP.mat[, i] = or.sim.tmp["TP", ]
  }
  if(has.coat){
    # coat
    coat.sim.tmp = readRDS(paste0(
      output_dir, "/roccurves", "/coat_roc", i, file.end
    ))
    coat.roc.list[[i]] = coat.sim.tmp
    coat.TPR.mat[, i] = coat.sim.tmp["tpr", ]
    coat.S.hat.mat[, i] = coat.sim.tmp["S_hat", ]
    coat.TP.mat[, i] = coat.sim.tmp["TP", ]
  }
  if(has.propr){
    # oracle
    pr.sim.tmp = readRDS(paste0(
      output_dir, "/roccurves", "/propr_roc", i, file.end
    ))
    pr.roc.list[[i]] = pr.sim.tmp
    pr.TPR.mat[, i] = pr.sim.tmp["tpr", ]
    pr.S.hat.mat[, i] = pr.sim.tmp["S_hat", ]
    pr.TP.mat[, i] = pr.sim.tmp["TP", ]
  }
}

# average over each possible S.hat/TP value
# stack columns so which() is more interpretable
cl.TPR.vec = as.vector(cl.TPR.mat)
cl.S.hat.vec = as.vector(cl.S.hat.mat)
cl.TP.vec = as.vector(cl.TP.mat)
slr.TPR.vec = as.vector(slr.TPR.mat)
slr.S.hat.vec = as.vector(slr.S.hat.mat)
slr.TP.vec = as.vector(slr.TP.mat)
if(has.oracle){
  or.TPR.vec = as.vector(or.TPR.mat)
  or.S.hat.vec = as.vector(or.S.hat.mat)
  or.TP.vec = as.vector(or.TP.mat)
}
if(has.coat){
  coat.TPR.vec = as.vector(coat.TPR.mat)
  coat.S.hat.vec = as.vector(coat.S.hat.mat)
  coat.TP.vec = as.vector(coat.TP.mat)
}
if(has.propr){
  pr.TPR.vec = as.vector(pr.TPR.mat)
  pr.S.hat.vec = as.vector(pr.S.hat.mat)
  pr.TP.vec = as.vector(pr.TP.mat)
}

# get the averages
S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec)))
if(has.oracle & !has.coat & !has.propr){
  S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec)))
} else if(has.oracle & !has.coat & has.propr){
  S.hat.vals = sort(unique(c(
    cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec, pr.S.hat.vec)))
}
cl.TPR.avg = rep(NA, length(S.hat.vals))
cl.TP.avg = rep(NA, length(S.hat.vals))
slr.TPR.avg = rep(NA, length(S.hat.vals))
slr.TP.avg = rep(NA, length(S.hat.vals))
if(has.oracle){
  or.TPR.avg = rep(NA, length(S.hat.vals))
  or.TP.avg = rep(NA, length(S.hat.vals))
}
if(has.coat){
  coat.TPR.avg = rep(NA, length(S.hat.vals))
  coat.TP.avg = rep(NA, length(S.hat.vals))
}
if(has.propr){
  pr.TPR.avg = rep(NA, length(S.hat.vals))
  pr.TP.avg = rep(NA, length(S.hat.vals))
}
for(i in 1:length(S.hat.vals)){
  val.tmp = S.hat.vals[i]
  # classo
  cl.which.idx.tmp = which(cl.S.hat.vec == val.tmp)
  cl.TPR.avg[i] = mean(cl.TPR.vec[cl.which.idx.tmp])
  cl.TP.avg[i] = mean(cl.TP.vec[cl.which.idx.tmp])
  # slr
  slr.which.idx.tmp = which(slr.S.hat.vec == val.tmp)
  slr.TPR.avg[i] = mean(slr.TPR.vec[slr.which.idx.tmp])
  slr.TP.avg[i] = mean(slr.TP.vec[slr.which.idx.tmp])
  if(has.oracle){
    # oracle
    or.which.idx.tmp = which(or.S.hat.vec == val.tmp)
    or.TPR.avg[i] = mean(or.TPR.vec[or.which.idx.tmp])
    or.TP.avg[i] = mean(or.TP.vec[or.which.idx.tmp])
  }
  if(has.coat){
    # coat
    coat.which.idx.tmp = which(coat.S.hat.vec == val.tmp)
    coat.TPR.avg[i] = mean(coat.TPR.vec[coat.which.idx.tmp])
    coat.TP.avg[i] = mean(coat.TP.vec[coat.which.idx.tmp])
  }
  if(has.propr){
    # propr
    pr.which.idx.tmp = which(pr.S.hat.vec == val.tmp)
    pr.TPR.avg[i] = mean(pr.TPR.vec[pr.which.idx.tmp])
    pr.TP.avg[i] = mean(pr.TP.vec[pr.which.idx.tmp])
  }
}

# plot
data.gg = rbind(
  data.frame(
    S_hat = S.hat.vals, TPR = cl.TPR.avg, TP = cl.TP.avg, Method = "classo"), 
  data.frame(
    S_hat = S.hat.vals, TPR = slr.TPR.avg, TP = slr.TP.avg, Method = "slr")
)
if(has.oracle){
  data.gg = rbind(
    data.gg, 
    data.frame(
      S_hat = S.hat.vals, TPR = or.TPR.avg, TP = or.TP.avg, Method = "oracle"))
}
if(has.coat){
  data.gg = rbind(
    data.gg, 
    data.frame(
      S_hat = S.hat.vals, TPR = coat.TPR.avg, TP = coat.TP.avg, Method = "coat"))
}
if(has.propr){
  data.gg = rbind(
    data.gg, 
    data.frame(
      S_hat = S.hat.vals, TPR = pr.TPR.avg, TP = pr.TP.avg, Method = "propr"))
}
data.gg$Method = factor(data.gg$Method, levels = levels.gg)
tp_roc = ggplot(
  data.gg[!is.na(data.gg$TPR),], aes(x = S_hat, y = TP, color = Method)) + 
  geom_line(alpha = 0.5, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
tpr_roc = ggplot(
  data.gg[!is.na(data.gg$TPR),], aes(x = S_hat, y = TPR, color = Method)) + 
  geom_line(alpha = 0.5, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
ggarrange(tp_roc, tpr_roc)
if(sigma.settings == "lin14Sigma" & mu.settings == "matchbeta"){
  ggsave(
    filename = paste0(
      "20211011_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_", mu.settings, "_rocs.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
} else{
  ggsave(
    filename = paste0(
      "20211011_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_rocs.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
}


