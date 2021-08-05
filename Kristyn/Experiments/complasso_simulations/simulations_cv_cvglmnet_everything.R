# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit supervised log-ratios method to the data
# Date: 06/30/2021
# Notes: cvglmnet for beta and theta, like how I (think) I was doing it

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output_cvglmnetall"

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)
library(balance) # for sbp.fromHclust()

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

# Dr. Ma sources
library(Matrix)
library(glmnet)
library(compositions)
library(stats)
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "supervisedlogratiosalpha.R"))

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

# Other simulation settings
numSims = 100
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
muW = c(
  rep(log(p), 5), 
  rep(0, p - 5)
)
SigmaW = matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    SigmaW[i, j] = rho^abs(i - j)
  }
}

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
evals = foreach(
  b = 1:numSims, 
  .combine = cbind
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
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
  muW = c(
    rep(log(p), 5), 
    rep(0, p - 5)
  )
  SigmaW = matrix(0, p, p)
  for(i in 1:p){
    for(j in 1:p){
      SigmaW[i, j] = rho^abs(i - j)
    }
  }
  
  nlam = 100
  
  ##############################################################################
  
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
  
  # # simulate training data #
  # # generate W
  # W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  # V = exp(W)
  # rowsumsV = apply(V, 1, sum)
  # X = V / rowsumsV
  # epsilon = rnorm(n, 0, sigma_eps)
  # Z = log(X)
  # # generate Y
  # Y = Z %*% beta + epsilon
  # 
  # # simulate test data #
  # # simulate independent test set of size n
  # # generate W
  # W.test = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  # V.test = exp(W.test)
  # rowsumsV.test = apply(V.test, 1, sum)
  # X.test = V.test / rowsumsV.test
  # epsilon.test = rnorm(n, 0, sigma_eps)
  # Z.test = log(X.test)
  # # generate Y
  # Y.test = Z.test %*% beta + epsilon.test
  
  # Simulate data like in aim31 #########################################################
  logW <- mvrnorm(n=n*2, mu=muW, Sigma=SigmaW) 
  W <- exp(logW) # basis
  colnames(W) <- paste0('s',1:p)
  XAll <- sweep(W,1,rowSums(W), FUN='/')
  # WCLR <- t(apply(WC,1,clr))
  names(beta) <- colnames(W)
  yAll <-  log(XAll) %*% beta + rnorm(n) * sigma_eps
  Y <- yAll[1:n,]
  Y.test <- yAll[-(1:n),]
  X = XAll[1:n, ]
  X.test = XAll[-(1:n), ]
  ZAll = log(XAll)
  Z = ZAll[1:n, ]
  Z.test = ZAll[-(1:n)]
  
  ##############################################################################
  
  ## Tuning parameters, like in aim31 ###################################################
  nlam <- 100
  maxlam <- 2*max(abs(crossprod(log(XAll[1:n,]),Y)/n))
  lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))
  
  # apply supervised log-ratios, using CV to select lambda
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept,
              lambda = lambda, # like in aim31 ##########################################
              rho.type = rho.type, linkage = linkage, standardize = scaling)
  btree = slr$btree
  # plot(btree)
  
  # choose lambda
  slr.lam.min.idx = which.min(slr$cvm)
  slr.lam.min = slr$lambda[slr.lam.min.idx]
  # slr.a0 = slr$int[slr.lam.min.idx]
  # slr.thetahat = slr$bet[, slr.lam.min.idx]
  slr.Uhat = getU(btree = btree)
  # slr.betahat = getBeta(slr.thetahat, U = slr.Uhat)
  
  # refit, like in run.glmnet() #########################################################
  ba_supervised = computeBalances(XAll, btree = btree)
  cv_exact = slr$cv.glmnet
  refit_exact <- glmnet(
    x=ba_supervised[1:n,],y=Y,family='gaussian',lambda=cv_exact$lambda.min)
  pred_exact <- predict(
    cv_exact, newx = ba_supervised[-(1:n),], type = "response", 
    s = 'lambda.min')
  mse.pred = mean((pred_exact-Y.test)^2) # same as aim31, when settings are same ########
  # get a0, thetahat, betahat
  slr.a0 = cv_exact$glmnet.fit$a0[which.min(cv_exact$cvm)] # slr$int[slr.lam.min.idx] # used for predict in aim31
  slr.thetahat = cv_exact$glmnet.fit$beta[, which.min(cv_exact$cvm)] # slr$bet[, slr.lam.min.idx] # used for predict in aim31
  
  # evaluate model #
  
  # 1. prediction error #
  # 1a. on training set #
  # get prediction error on training set
  slr.Yhat.train = slr.a0 + computeBalances(X, btree) %*% slr.thetahat
  # slr.Yhat.train = refit_exact$a0 +
  #   computeBalances(X, btree = btree) %*% refit_exact$beta # same
  slr.PE.train = as.vector(crossprod(Y - slr.Yhat.train) / n)
  # 1b. on test set #
  # get prediction error on test set
  slr.Yhat.test = slr.a0 + computeBalances(X.test, btree) %*% slr.thetahat
  # slr.Yhat.test = refit_exact$a0 + 
  #   computeBalances(X.test, btree = btree) %*% refit_exact$beta
  # slr.PE.test = as.vector(crossprod(Y.test - slr.Yhat.test) / n)
  slr.PE.test = mean((slr.Yhat.test-Y.test)^2) # same as mse.pred in aim31
  # 2. estimation accuracy #
  # slr.thetahat = as.matrix(refit_exact$beta) # slr$bet[, slr.lam.min.idx] # used refit_exact instead of cv_exact for evaluating in aim31
  slr.betahat = getBeta(slr.thetahat, U = slr.Uhat) # used refit_exact instead of cv_exact for evaluating in aim31
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
  file = paste0(output_dir, "/slr_sim", b, file.end))
  
  ##############################################################################
  
  # apply compositional lasso, using CV to select lambda
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = Z, Cmat = matrix(1, p, 1), nlam = nlam, 
    lambda = lambda, # like aim31 #######################################################
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  
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
  
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  cl.EA1 = sum(abs(cl.betahat - beta))
  cl.EA2 = as.vector(sqrt(crossprod(cl.betahat - beta)))
  cl.EAInfty = max(abs(cl.betahat - beta))
  
  # 3. selection accuracy #
  # 3a. selection of beta #
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
  file = paste0(output_dir, "/complasso_sim", b, file.end))
  
}


