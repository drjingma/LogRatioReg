# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Date: 09/23/2021

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
  getMSEyhat = function(y, n, betahat0, betahat, predMat){
    yhat = betahat0 + predMat %*% betahat
    mse = as.vector(crossprod(y - yhat) / n)
    return(mse)
  }
  getEstimationAccuracy = function(true.beta, betahat){
    EA1 = sum(abs(betahat - true.beta))
    EA2 = as.vector(sqrt(crossprod(betahat - true.beta)))
    EAInfty = max(abs(betahat - true.beta))
    return(list(
      EA1 = EA1, 
      EA2 = EA2, 
      EAInfty = EAInfty
    ))
  }
  getSelectionAccuracy = function(is0.true.beta, non0.true.beta, non0.betahat){
    FP = sum(is0.true.beta & non0.betahat)
    FN = sum((non0.true.beta != non0.betahat) & non0.true.beta)
    TPR = sum((non0.true.beta == non0.betahat) & non0.betahat) / 
      sum(non0.true.beta)
    # F-score = precision / recall
    # precision = # true positive results / # of positive results
    #   (including those not identified correctly)
    precision = sum((non0.true.beta == non0.betahat) & non0.betahat) / 
      sum(non0.betahat) 
    # recall = sensitivity = TPR = # true positive results / # of true positives
    Fscore = 2 * precision * TPR / (precision + TPR)
    return(list(
      FP = FP, 
      FN = FN, 
      TPR = TPR, 
      precision = precision, 
      Fscore = Fscore
    ))
  }
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
  sigma.settings = "4blockSigma"
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "2blocks"  
  # "2blocks" => choose j corresp. to two blocks
  #   (one block w/-1, other w/1 in a single contrast)
  # "2blocks2contrasts" => choose two j's corresp to two separate blocks
  # "1block" => choose j corresp. to one block only
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
  sigma_eps = 0.5
  muW = c(rep(log(p), 5), rep(0, p - 5))
  SigmaWblock = matrix(cor_ij, p / 4, p / 4)
  for(i in 1:nrow(SigmaWblock)) SigmaWblock[i, i] = 1
  SigmaW0 = matrix(0, p / 4, p / 4)
  SigmaW = cbind(
    rbind(SigmaWblock, SigmaW0, SigmaW0, SigmaW0), 
    rbind(SigmaW0, SigmaWblock, SigmaW0, SigmaW0), 
    rbind(SigmaW0, SigmaW0, SigmaWblock, SigmaW0),  
    rbind(SigmaW0, SigmaW0, SigmaW0, SigmaWblock)
  )
  SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
  U = getU(btree = SigmaWtree) # transformation matrix
  # plot(SigmaWtree)
  # sbp.fromHclust(SigmaWtree)
  
  # theta settings
  if(theta.settings == "dense"){
    indices.theta = 1
  } else if(theta.settings == "2blocks"){
    SBP = sbp.fromHclust(SigmaWtree)
    # for each column (contrast), find which variables are included (1 or -1)
    contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
    # get the contrasts with length p / 2 -- have 2 blocks of correlated vars
    #   not necessary, but may save on unnecessary computation in the next step
    block.contrasts = which(sapply(contrast.vars, length) == p / 2)
    # pick one such contrast
    if(length(block.contrasts) == 1){
      indices.theta = unname(block.contrasts)
    } else{
      indices.theta = unname(sample(x = block.contrasts, 1))
    }
  } else if(theta.settings == "1block"){
    SBP = sbp.fromHclust(SigmaWtree)
    # for each column (contrast), find which variables are included (1 or -1)
    contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
    # get the contrasts with length p / 2 -- have 2 blocks of correlated vars
    #   not necessary, but may save on unnecessary computation in the next step
    block.contrasts = which(sapply(contrast.vars, length) == p / 4)
    # pick one such contrast
    if(length(block.contrasts) == 1){
      indices.theta = unname(block.contrasts)
    } else{
      indices.theta = unname(sample(x = block.contrasts, 1))
    }
  } else if(theta.settings == "2blocks2contrasts"){
    SBP = sbp.fromHclust(SigmaWtree)
    # for each column (contrast), find which variables are included (1 or -1)
    contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
    # get the contrasts with length p / 2 -- have 2 blocks of correlated vars
    #   not necessary, but may save on unnecessary computation in the next step
    block.contrasts = which(sapply(contrast.vars, length) == p / 4)
    # pick one such contrast
    if(length(block.contrasts) < 2){
      stop("need at least 2 contrasts corresponding to two separate blocks")
    } else{
      indices.theta = unname(sample(x = block.contrasts, 2))
    }
  } else{
    stop("invalid theta.settings")
  }
  # error checking indices.theta found based on theta.settings argument
  if(is.null(indices.theta) | length(indices.theta) != 1){
    stop("invalid indices.theta")
  }
  
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
  colnames(W) <- paste0('s', 1:p)
  XAll <- sweep(W, 1, rowSums(W), FUN='/')
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
  Y <- yAll[1:n, , drop = TRUE]
  Y.test <- yAll[-(1:n), , drop = TRUE]
  
  # about beta
  non0.beta = (beta != 0)
  slr.non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
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
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
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


