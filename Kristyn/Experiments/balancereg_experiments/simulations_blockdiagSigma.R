# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Date: 09/13/2021

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
  roc.for.coef <- function(beta_hat, beta, eps = 1e-08){
    TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
    FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
    tpr <- TP/(TP + FN)
    S_hat <- sum((abs(beta_hat) > eps))
    out <- c(S_hat,tpr)
    names(out) <- c('S_hat','tpr')
    return(out)
  }
  
  roc.for.coef.LR <- function(beta_hat,beta,sbp,eps=1e-08){
    
    if (is.null(sbp)){
      stop('A sequential binary partition tree is needed for roc evaluation!')
    }
    
    # first identify the variable at the LR scale
    index <- which(abs(beta_hat) > eps)
    
    # map to original variable
    if (length(index)==0){
      S_hat <- NULL
    } else  if (length(index)==1){
      S_hat <- names(which(abs(sbp[,index])>0))
    } else {
      S_hat <- names(which(rowSums(abs(sbp[,index]))>0))
    }
    S0 <- names(which((abs(beta) > eps)))
    TP <- intersect(S_hat, S0)
    tpr <- length(TP)/length(S0)
    out <- c(length(S_hat),tpr)
    names(out) <- c('S_hat','tpr')
    return(out)
  }
  
  # Settings to toggle with
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "block" 
  # "dense" => j = 1 (of theta = theta_1, ..., theta_j, ..., theta_{p-1})
  # "sparse" => j = p - 1
  # "both" => theta = c(1, p - 1)
  # "multsparse" = c(p - 3:1)
  
  # "block" => choose the block of correlated variables in Gamma2
  # "pair+block" => choose one pair of correlated variable in Gamma1 +
  #   the block of correlated variables in Gamma2
  values.theta = 1
  linkage = "average"
  tol = 1e-4
  nlam = 100
  intercept = TRUE
  K = 10
  n = 100
  p = 200
  rho = 0.5 # 0.2, 0.5
  scaling = TRUE
  
  # theta settings
  if(theta.settings == "dense"){
    indices.theta = 1
  } else if(theta.settings == "sparse"){
    indices.theta = p - 1
  } else if(theta.settings == "both"){
    indices.theta = c(1, p - 1)
  } else if(theta.settings == "multsparse"){
    indices.theta = p - 3:1
  }
  
  # Population parameters
  sigma_eps = 0.5
  muW = c(rep(log(p), 5), rep(0, p - 5))
  SigmaW11 = matrix(0.9, p / 2, p / 2)
  for(i in 1:nrow(SigmaW11)) SigmaW11[i, i] = 1
  SigmaW12 = matrix(0, p / 2, p / 2)
  SigmaW = cbind(rbind(SigmaW11, SigmaW12), rbind(SigmaW12, SigmaW11))
  SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
  U = getU(btree = SigmaWtree) # transformation matrix
  
  # theta settings for block-diagonal Sigma
  SBP = sbp.fromHclust(SigmaWtree)
  block1vars = 1:(p / 2)
  block2vars = ((p / 2) + 1):p
  # for each column (contrast), find which variables are included (1 or -1)
  contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
  if(theta.settings == "block" | theta.settings == "pairblock"){
    # get the contrasts with length p / 2 -- there are 2 of them
    #   not necessary, but may save on unnecessary computation in the next step
    block.contrasts = which(sapply(contrast.vars, length) == p / 2)
    # find out which one contains block2vars
    is.block2vars.contrast = sapply(
      contrast.vars[block.contrasts], FUN = function(x) 
        isTRUE(all.equal(unname(sort(x)), block2vars)))
    block2vars.contrast = block.contrasts[is.block2vars.contrast]
    indices.theta = unname(block2vars.contrast)
    if(theta.settings == "pairblock"){
      # find pair, too
      #   again, not necessary, but saves computation
      pair.contrasts = which(sapply(contrast.vars, length) == 2)
      # find out which one is a pair of variables in block1
      is.block1pairvars.contrast = sapply(
        contrast.vars[pair.contrasts], FUN = function(x) all(x %in% block1vars))
      block1pairvars.contrast = pair.contrasts[is.block1pairvars.contrast]
      indices.theta = unname(c(block2vars.contrast, block1pairvars.contrast))
    }
  } else{
    stop("invalid theta.settings")
  }
  
  file.end = paste0(
    "_dim", n, "x", p, 
    "_", theta.settings, 
    "_rho", rho, 
    "_int", intercept,
    "_scale", scaling,
    "_K", K,
    "_seed", rng.seed,
    ".rds")
  
  ##############################################################################
  
  # simulate data
  logW <- mvrnorm(n = n * 2, mu = muW, Sigma = SigmaW) 
  W <- exp(logW)
  colnames(W) <- paste0('s', 1:p)
  XAll <- sweep(W, 1, rowSums(W), FUN='/')
  ZAll = log(XAll)
  ilrXAll = computeBalances(XAll, U = U)
  
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
  # get beta
  beta = as.vector(getBeta(theta, U = U))
  names(beta) <- colnames(W)
  # generate Y
  yAll = ilrXAll %*% theta + rnorm(n) * sigma_eps
  
  # subset out training and test sets
  X = XAll[1:n, ]
  X.test = XAll[-(1:n), ]
  Z = ZAll[1:n, ]
  Z.test = ZAll[-(1:n), ]
  ilrX = ilrXAll[1:n, ]
  ilrX.test = ilrXAll[-(1:n), ]
  Y <- yAll[1:n, , drop = TRUE]
  Y.test <- yAll[-(1:n), , drop = TRUE]
  
  # about beta
  non0.beta = (beta != 0)
  slr.non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  beta.active = beta[non0.beta]
  beta.inactive = beta[is0.beta]
  
  ##############################################################################
  # supervised log-ratios
  ##############################################################################
  
  # apply supervised log-ratios, using CV to select lambda
  start.time = Sys.time()
  slr = cvSLR(
    y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
    rho.type = rho.type, linkage = linkage, standardize = scaling)
  end.time = Sys.time()
  slr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
  slr.btree = slr$btree
  # plot(slr.btree)
  
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
  # get prediction error on training set
  slr.Yhat.train = slr.a0 + computeBalances(X, slr.btree) %*% slr.thetahat
  slr.PE.train = as.vector(crossprod(Y - slr.Yhat.train) / n)
  # 1b. on test set #
  # get prediction error on test set
  slr.Yhat.test = slr.a0 + computeBalances(X.test, slr.btree) %*% slr.thetahat
  slr.PE.test = as.vector(crossprod(Y.test - slr.Yhat.test) / n)
  
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  slr.EA1 = sum(abs(slr.betahat - beta))
  slr.EA2 = as.vector(sqrt(crossprod(slr.betahat - beta)))
  slr.EAInfty = max(abs(slr.betahat - beta))
  # 2b. estimation accuracy for active set
  slr.betahat.active = slr.betahat[non0.beta]
  slr.EA1.active = sum(abs(slr.betahat.active - beta.active))
  slr.EA2.active = as.vector(sqrt(crossprod(
    slr.betahat.active - beta.active)))
  slr.EAInfty.active = max(abs(slr.betahat.active - beta.active))
  # 2c. estimation accuracy for inactive set
  slr.betahat.inactive = slr.betahat[is0.beta]
  slr.EA1.inactive = sum(abs(slr.betahat.inactive - beta.inactive))
  slr.EA2.inactive = as.vector(sqrt(crossprod(
    slr.betahat.inactive - beta.inactive)))
  slr.EAInfty.inactive = max(abs(slr.betahat.inactive - beta.inactive))
  
  # 3. selection accuracy #
  # 3a. selection of beta #
  ### using SBP matrix
  slr.SBP = sbp.fromHclust(slr.btree)
  slr.non0.thetahat = (slr.thetahat != 0)
  slr.sel.cols.SBP = slr.SBP[, slr.non0.thetahat, drop = FALSE]
  slr.non0.betahat = apply(slr.sel.cols.SBP, 1, function(row) any(row != 0))
  #
  slr.FP = sum(is0.beta & slr.non0.betahat)
  slr.FN = sum((non0.beta != slr.non0.betahat) & non0.beta)
  slr.TPR = sum((non0.beta == slr.non0.betahat) & slr.non0.betahat) / 
    sum(non0.beta)
  
  saveRDS(c(
    "PEtr" = slr.PE.train, 
    "PEte" = slr.PE.test, 
    "EA1" = slr.EA1, 
    "EA2" = slr.EA2, 
    "EAInfty" = slr.EAInfty, 
    "EA1Active" = slr.EA1.active, 
    "EA2Active" = slr.EA2.active, 
    "EAInftyActive" = slr.EAInfty.active, 
    "EA1Inactive" = slr.EA1.inactive, 
    "EA2Inactive" = slr.EA2.inactive, 
    "EAInftyInactive" = slr.EAInfty.inactive, 
    "FP" = slr.FP, 
    "FN" = slr.FN, 
    "TPR" = slr.TPR, 
    "timing" = slr.timing,
    "betaSparsity" = bspars
  ), 
  file = paste0(output_dir, "/slr_metrics", b, file.end))
  
  # roc
  slr.roc <- apply(slr$bet, 2, function(a) 
    roc.for.coef.LR(a, beta, slr.SBP))
  
  saveRDS(slr.roc, file = paste0(output_dir, "/slr_roc", b, file.end))
  
  ##############################################################################
  # compositional lasso
  ##############################################################################
  
  # apply compositional lasso, using CV to select lambda
  start.time = Sys.time()
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = Z, Cmat = matrix(1, p, 1), nlam = nlam, 
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  end.time = Sys.time()
  classo.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
  
  # choose lambda
  cl.lam.min.idx = which.min(complasso$cvm)
  cl.lam.min = complasso$lambda[cl.lam.min.idx]
  cl.a0 = complasso$int[cl.lam.min.idx]
  cl.betahat = complasso$bet[, cl.lam.min.idx]
  
  # evaluate model #
  
  # 1. prediction error #
  # 1a. on training set #
  # get prediction error on training set
  cl.Yhat.train = cl.a0 + log(X) %*% cl.betahat
  cl.PE.train = as.vector(crossprod(Y - cl.Yhat.train) / n)
  # 1b. on test set #
  # get prediction error on test set
  cl.Yhat.test = cl.a0 + log(X.test) %*% cl.betahat
  cl.PE.test = as.vector(crossprod(Y.test - cl.Yhat.test) / n)
  
  # 2. estimation accuracy (i.t.o. beta) #
  cl.EA1 = sum(abs(cl.betahat - beta))
  cl.EA2 = as.vector(sqrt(crossprod(cl.betahat - beta)))
  cl.EAInfty = max(abs(cl.betahat - beta))
  # 2b. estimation accuracy for active set
  cl.betahat.active = cl.betahat[non0.beta]
  cl.EA1.active = sum(abs(cl.betahat.active - beta.active))
  cl.EA2.active = as.vector(sqrt(crossprod(
    cl.betahat.active - beta.active)))
  cl.EAInfty.active = max(abs(cl.betahat.active - beta.active))
  # 2c. estimation accuracy for inactive set
  cl.betahat.inactive = cl.betahat[is0.beta]
  cl.EA1.inactive = sum(abs(cl.betahat.inactive - beta.inactive))
  cl.EA2.inactive = as.vector(sqrt(crossprod(
    cl.betahat.inactive - beta.inactive)))
  cl.EAInfty.inactive = max(abs(cl.betahat.inactive - beta.inactive))
  
  # 3. selection accuracy (i.t.o. beta) #
  cl.non0.betahat = abs(cl.betahat) > 10e-8
  cl.FP = sum(is0.beta & cl.non0.betahat)
  cl.FN = sum((non0.beta != cl.non0.betahat) & non0.beta)
  cl.TPR = sum((non0.beta == cl.non0.betahat) & cl.non0.betahat) / 
    sum(non0.beta)
  
  saveRDS(c(
    "PEtr" = cl.PE.train, 
    "PEte" = cl.PE.test, 
    "EA1" = cl.EA1, 
    "EA2" = cl.EA2, 
    "EAInfty" = cl.EAInfty, 
    "EA1Active" = cl.EA1.active, 
    "EA2Active" = cl.EA2.active, 
    "EAInftyActive" = cl.EAInfty.active, 
    "EA1Inactive" = cl.EA1.inactive, 
    "EA2Inactive" = cl.EA2.inactive, 
    "EAInftyInactive" = cl.EAInfty.inactive, 
    "FP" = cl.FP, 
    "FN" = cl.FN, 
    "TPR" = cl.TPR, 
    "timing" = classo.timing,
    "betaSparsity" = bspars
  ), 
  file = paste0(output_dir, "/classo_metrics", b, file.end))
  
  # roc
  cl.roc <- apply(complasso$bet, 2, function(a) 
    roc.for.coef(a, beta))
  
  saveRDS(cl.roc, file = paste0(output_dir, "/classo_roc", b, file.end))
  
  ##############################################################################
  # oracle method
  ##############################################################################
  
  # apply oracle method, using CV to select lambda
  start.time = Sys.time()
  oracle = cvILR(y = Y, X = X, btree = SigmaWtree, U = U, nlam = nlam, 
                 nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  or.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
  or.btree = SigmaWtree
  # plot(or.btree)
  
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
  # get prediction error on training set
  or.Yhat.train = or.a0 + computeBalances(X, or.btree) %*% or.thetahat
  or.PE.train = as.vector(crossprod(Y - or.Yhat.train) / n)
  # 1b. on test set #
  # get prediction error on test set
  or.Yhat.test = or.a0 + computeBalances(X.test, or.btree) %*% or.thetahat
  or.PE.test = as.vector(crossprod(Y.test - or.Yhat.test) / n)
  
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  or.EA1 = sum(abs(or.betahat - beta))
  or.EA2 = as.vector(sqrt(crossprod(or.betahat - beta)))
  or.EAInfty = max(abs(or.betahat - beta))
  # 2b. estimation accuracy for active set
  or.betahat.active = or.betahat[non0.beta]
  or.EA1.active = sum(abs(or.betahat.active - beta.active))
  or.EA2.active = as.vector(sqrt(crossprod(
    or.betahat.active - beta.active)))
  or.EAInfty.active = max(abs(or.betahat.active - beta.active))
  # 2c. estimation accuracy for inactive set
  or.betahat.inactive = or.betahat[is0.beta]
  or.EA1.inactive = sum(abs(or.betahat.inactive - beta.inactive))
  or.EA2.inactive = as.vector(sqrt(crossprod(
    or.betahat.inactive - beta.inactive)))
  or.EAInfty.inactive = max(abs(or.betahat.inactive - beta.inactive))
  
  # 3. selection accuracy #
  # 3a. selection of beta #
  ### using SBP matrix
  or.SBP = sbp.fromHclust(or.btree)
  row.names(or.SBP) = colnames(W)
  or.non0.thetahat = (or.thetahat != 0)
  or.sel.cols.SBP = or.SBP[, or.non0.thetahat, drop = FALSE]
  or.non0.betahat = apply(or.sel.cols.SBP, 1, function(row) any(row != 0))
  #
  or.FP = sum(is0.beta & or.non0.betahat)
  or.FN = sum((non0.beta != or.non0.betahat) & non0.beta)
  or.TPR = sum((non0.beta == or.non0.betahat) & or.non0.betahat) / 
    sum(non0.beta)
  
  saveRDS(c(
    "PEtr" = or.PE.train, 
    "PEte" = or.PE.test, 
    "EA1" = or.EA1, 
    "EA2" = or.EA2, 
    "EAInfty" = or.EAInfty, 
    "EA1Active" = or.EA1.active, 
    "EA2Active" = or.EA2.active, 
    "EAInftyActive" = or.EAInfty.active, 
    "EA1Inactive" = or.EA1.inactive, 
    "EA2Inactive" = or.EA2.inactive, 
    "EAInftyInactive" = or.EAInfty.inactive, 
    "FP" = or.FP, 
    "FN" = or.FN, 
    "TPR" = or.TPR, 
    "timing" = or.timing,
    "betaSparsity" = bspars
  ), 
  file = paste0(output_dir, "/oracle_metrics", b, file.end))
  
  # roc
  or.roc <- apply(oracle$bet, 2, function(a) 
    roc.for.coef.LR(a, beta, or.SBP))
  
  saveRDS(or.roc, file = paste0(output_dir, "/oracle_roc", b, file.end))
  
  # ##############################################################################
  # # selbal
  # ##############################################################################
  # rownames(X) = paste("Sample", 1:nrow(X), sep = "_")
  # colnames(X) = paste("V", 1:ncol(X), sep = "")
  # 
  # # apply selbal, using CV to select the optimal number of variables
  # start.time = Sys.time()
  # selbal.fit = selbal.cv(x = X, y = as.vector(Y), n.fold = K)
  # end.time = Sys.time()
  # selbal.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # U (transformation) matrix
  # u.selbal = rep(0, p)
  # names(u.selbal) = colnames(X)
  # pba.pos = unlist(subset(
  #   selbal.fit$global.balance, subset = Group == "NUM", select = Taxa))
  # r = length(pba.pos)
  # pba.neg = unlist(subset(
  #   selbal.fit$global.balance, subset = Group == "DEN", select = Taxa))
  # s = length(pba.neg)
  # u.selbal[pba.pos] = 1 / r
  # u.selbal[pba.neg] = -1 / s
  # norm.const = sqrt((r * s) / (r + s))
  # u.selbal = norm.const * u.selbal
  # # check: these are equal
  # # lm(as.vector(Y) ~ log(X) %*% as.matrix(u.selbal))
  # # selbal.fit$glm
  # selbal.thetahat = coefficients(selbal.fit$glm)[2]
  # selbal.betahat = u.selbal %*% as.matrix(selbal.thetahat)
  # 
  # # evaluate model #
  # 
  # # 1. prediction error #
  # # 1a. on training set #
  # # get prediction error on training set
  # selbal.Yhat.train = predict.glm(selbal.fit$glm, 
  #                          newdata = data.frame(X), 
  #                          type = "response")
  # selbal.PE.train = crossprod(Y - selbal.Yhat.train) / n
  # # 1b. on test set #
  # # get prediction error on test set
  # rownames(X.test) = paste("Sample", 1:nrow(X), sep = "_")
  # colnames(X.test) = paste("V", 1:ncol(X), sep = "")
  # selbal.Yhat.test = predict.glm(selbal.fit$glm, newdata = data.frame(X.test), 
  #                         type = "response")
  # selbal.PE.test = crossprod(Y.test - selbal.Yhat.test) / n
  # 
  # # 2. estimation accuracy (i.t.o. beta) #
  # selbal.EA1 = sum(abs(selbal.betahat - beta))
  # selbal.EA2 = sqrt(crossprod(selbal.betahat - beta))
  # selbal.EAInfty = max(abs(selbal.betahat - beta))
  # 
  # # 3. selection accuracy (i.t.o. beta) #
  # selbal.non0.betahat = abs(selbal.betahat) > 10e-8
  # selbal.is0.betahat = selbal.betahat <= 10e-8
  # # FP
  # selbal.FP = sum(is0.beta & selbal.non0.betahat)
  # # FN
  # selbal.FN = sum((non0.beta != cl.non0.betahat) & non0.beta)
  # # TPR
  # selbal.TPR = sum((non0.beta == selbal.non0.betahat) & selbal.non0.betahat) / 
  #   sum(non0.beta)
  # 
  # saveRDS(c(
  #   "PEtr" = selbal.PE.train, 
  #   "PEte" = selbal.PE.test, 
  #   "EA1" = selbal.EA1, 
  #   "EA2" = selbal.EA2, 
  #   "EAInfty" = selbal.EAInfty, 
  #   "FP" = selbal.FP, 
  #   "FN" = selbal.FN, 
  #   "TPR" = selbal.TPR, 
  #   "timing" = selbal.timing,
  #   "betaSparsity" = bspars
  # ), 
  # file = paste0(output_dir, "/selbal_metrics", b, file.end))

  # roc
  # selbal.roc <- apply(slr$bet, 2, function(a) 
  #   roc.for.coef.LR(a, beta, ...))
  # 
  # saveRDS(selbal.roc, file = paste0(output_dir, "/selbal_roc", b, file.end))
}


