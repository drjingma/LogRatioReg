# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit supervised log-ratios method to the data
# Date: 08/29/2021
# Notes: cv.glmnet for predictions, but beta from refit.glmnet, like in aim31

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output_kristyn_metrics_roc"

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
  beta.settings = "new"
  linkage = "average"
  tol = 1e-4
  nlam = 100
  intercept = TRUE
  K = 10
  n = 100
  p = 200
  rho = 0.5 # 0.2, 0.5
  scaling = TRUE
  
  # which beta?
  if(beta.settings == "old" | beta.settings == "linetal2014"){
    beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
  } else{
    beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
  }
  
  # Population parameters
  non0.beta = (beta != 0)
  sigma_eps = 0.5
  seed = 1
  muW = c(rep(log(p), 5), rep(0, p - 5))
  SigmaW <- rgExpDecay(p,rho)$Sigma
  
  # for beta selection accuracy metrics
  slr.non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  
  # beta sparsity
  bspars = sum(non0.beta)
  
  file.end = paste0(
    "_dim", n, "x", p, 
    "_", beta.settings, 
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
  names(beta) <- colnames(W)
  yAll <-  ZAll %*% beta + rnorm(n) * sigma_eps
  
  # subset out training and test sets
  X = XAll[1:n, ]
  X.test = XAll[-(1:n), ]
  Z = ZAll[1:n, ]
  Z.test = ZAll[-(1:n)] 
  Y <- yAll[1:n,]
  Y.test <- yAll[-(1:n),]
  
  ##############################################################################
  
  # tuning parameter
  nlam <- 100
  lambda <- exp(seq(from=2, to=log(1e-4), length.out=100))
  
  # apply supervised log-ratios, using CV to select lambda
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
              rho.type = rho.type, linkage = linkage, standardize = scaling, 
              lambda = lambda)
  btree = slr$btree
  # plot(btree)
  
  # choose lambda
  slr.lam.min.idx = which.min(slr$cvm)
  slr.lam.min = slr$lambda[slr.lam.min.idx]
  slr.a0 = slr$int[slr.lam.min.idx]
  slr.thetahat = slr$bet[, slr.lam.min.idx]
  slr.Uhat = getU(btree = btree)
  slr.betahat = getBeta(slr.thetahat, U = slr.Uhat)
  
  # evaluate model #
  
  # 1. prediction error #
  # 1a. on training set #
  # get prediction error on training set
  slr.Yhat.train = slr.a0 + computeBalances(X, btree) %*% slr.thetahat
  slr.PE.train = as.vector(crossprod(Y - slr.Yhat.train) / n)
  # 1b. on test set #
  # get prediction error on test set
  slr.Yhat.test = slr.a0 + computeBalances(X.test, btree) %*% slr.thetahat
   
  slr.PE.test = as.vector(crossprod(Y.test - slr.Yhat.test) / n)
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  slr.EA1 = sum(abs(slr.betahat - beta))
  slr.EA2 = as.vector(sqrt(crossprod(slr.betahat - beta)))
  slr.EAInfty = max(abs(slr.betahat - beta))
  
  # 3. selection accuracy #
  # 3a. selection of beta #
  ### using SBP matrix
  slr.SBP = sbp.fromHclust(btree)
  slr.non0.thetahat = (slr.thetahat != 0)
  slr.sel.cols.SBP = slr.SBP[, slr.non0.thetahat]
  slr.non0.betahat = apply(slr.sel.cols.SBP, 1, function(row) any(row != 0))
  slr.is0.betahat = !slr.non0.betahat
  # FP
  slr.FP = sum(is0.beta & slr.non0.betahat)
  # FN
  slr.FN = sum((non0.beta != slr.non0.betahat) & non0.beta)
  # TPR
  slr.TPR = sum((non0.beta == slr.non0.betahat) & slr.non0.betahat) / 
    sum(non0.beta)
  
  saveRDS(c(
    "PEtr" = slr.PE.train, 
    "PEte" = slr.PE.test, 
    "EA1" = slr.EA1, 
    "EA2" = slr.EA2, 
    "EAInfty" = slr.EAInfty, 
    "FP" = slr.FP, 
    "FN" = slr.FN, 
    "TPR" = slr.TPR, 
    "betaSparsity" = bspars
  ), 
  file = paste0(output_dir, "/slr_metrics", b, file.end))
  
  # roc
  slr.roc <- apply(slr$bet, 2, function(a) 
    roc.for.coef.LR(a, beta, sbp.fromHclust(btree)))
  
  saveRDS(slr.roc, file = paste0(output_dir, "/slr_roc_samelam", b, file.end))
  
  ##############################################################################
  
  # apply compositional lasso, using CV to select lambda
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = Z, Cmat = matrix(1, p, 1), nlam = nlam, 
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling, 
    lambda = lambda)
  
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
  
  # 3. selection accuracy (i.t.o. beta) #
  cl.non0.betahat = abs(cl.betahat) > 10e-8
  cl.is0.betahat = cl.betahat <= 10e-8
  # FP
  cl.FP = sum(is0.beta & cl.non0.betahat)
  # FN
  cl.FN = sum((non0.beta != cl.non0.betahat) & non0.beta)
  # TPR
  cl.TPR = sum((non0.beta == cl.non0.betahat) & cl.non0.betahat) / 
    sum(non0.beta)
  
  saveRDS(c(
    "PEtr" = cl.PE.train, 
    "PEte" = cl.PE.test, 
    "EA1" = cl.EA1, 
    "EA2" = cl.EA2, 
    "EAInfty" = cl.EAInfty, 
    "FP" = cl.FP, 
    "FN" = cl.FN, 
    "TPR" = cl.TPR, 
    "betaSparsity" = bspars
  ), 
  file = paste0(output_dir, "/classo_metrics", b, file.end))
  
  # roc
  cl.roc <- apply(complasso$bet, 2, function(a) 
    roc.for.coef(a, beta)) # used cv.glmnet beta matrix
  
  saveRDS(cl.roc, file = paste0(output_dir, "/classo_roc_samelam", b, file.end))
}


