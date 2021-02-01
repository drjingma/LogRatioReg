# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit compositional Lasso to the data
# Date: 1/2021
# Notes: 

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output"

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)
library(stats) # for hclust()
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

# Method Settings
nlam = 200
intercept = TRUE

# Simulation settings
numSims = 100
n = 50
p = 30
rho = 0.2 # 0.2, 0.5
# beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
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

# set.seed(16) # leads to FN = 1
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
  library(balance) # for sbp.fromHclust()
  
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
  nlam = 200
  
  # simulate data
  # generate W
  W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V = exp(W)
  rowsumsV = apply(V, 1, sum)
  X = V / rowsumsV
  epsilon = rnorm(n, 0, sigma_eps)
  Z = log(X)
  # generate Y
  Y = Z %*% matrix(beta) + epsilon
  
  # apply compositional lasso, using GIC to select lambda
  slr0 = fitSLR(Y, X, nlam = nlam, intercept = intercept)
  # transformed
  betahats0.tr = apply(slr0$bet, 2, function(x) LRtoLC(x, slr0$btree))
  non0.betahats = (betahats0.tr != 0) # diff lambda = diff col
  non0ct = apply(non0.betahats, 2, sum) # count of non-zero betas for each lambda
  which0ct.leq3sqrtn = which(non0ct <= 3 * sqrt(n)) # lambda lower bound criteria
  lambda = seq(slr0$lambda[1],
               slr0$lambda[max(which0ct.leq3sqrtn)],
               length.out = nlam)
  nlam = length(lambda)
  slr = fitSLR(Y, X, lambda = lambda, nlam = nlam, intercept = intercept)
  betahats.tr = apply(slr$bet, 2, function(x) LRtoLC(x, slr$btree))
  non0.betahats = (betahats.tr != 0)
  non0ct = apply(non0.betahats, 2, sum)
  
  a0s = slr$int
  betahats = slr$bet
  btree = slr$btree
  GIC = rep(NA, nlam)
  pvn = max(p, n)
  for(m in 1:nlam){
    a0.tmp = a0s[m]
    betahat.tmp = betahats[, m]
    predictSLR.tmp = function(X){
      a0.tmp + computeBalances(X, btree) %*% betahat.tmp
    }
    Y.pred = predictSLR.tmp(X)
    sigmasq.hat = crossprod(Y - Y.pred) / n
    s.lam = sum(non0.betahats[, m])
    GIC[m] = log(sigmasq.hat) + (s.lam - 1) * log(log(n)) * log(pvn) / n
  }
  lam.min.idx = which.min(GIC)
  lam.min = slr$lambda[lam.min.idx]
  a0 = a0s[lam.min.idx]
  betahat = betahats[, lam.min.idx]
  non0.betahat = non0.betahats[, lam.min.idx]
  
  # plot lambda vs. GIC
  # plot(slr$lambda, GIC, type = "l")
  # points(lam.min, min(GIC))
  # text(lam.min, min(GIC),
  #      paste0("(", round(lam.min, 3), ", ", round(min(GIC), 3), ")", sep = ""),
  #      pos = 4)
  
  # evaluate model
  # prediction error
  # simulate independent test set of size n
  # generate W
  W.test = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V.test = exp(W.test)
  rowsumsV.test = apply(V.test, 1, sum)
  X.test = V.test / rowsumsV.test
  epsilon.test = rnorm(n, 0, sigma_eps)
  Z.test = log(X.test)
  # generate Y
  Y.test = Z.test %*% matrix(beta) + epsilon.test
  # get fitted model
  predictSLR = function(X){
    a0 + computeBalances(X, btree) %*% betahat
  }
  Y.pred = predictSLR(X.test)
  PE = crossprod(Y.test - Y.pred) / n
  # estimation accuracy
  # transform betahat to log-contrast space
  betahat.tr = LRtoLC(betahat, btree)
  # check that pred y is equal to that using transformed log-contrasts model
  # a0.tr = as.numeric(mean(Y) - colMeans(Z) %*% betahat.tr)
  # Y.pred.tr = a0.tr + Z.test %*% betahat.tr
  # all.equal(Y.pred, Y.pred.tr)
  # calculate other evaluation metrics
  EA1 = sum(abs(betahat.tr - beta))
  EA2 = crossprod(betahat.tr - beta)
  EAInfty = max(abs(betahat.tr - beta))
  non0.betahat = (betahat.tr != 0)
  # FP
  FP = sum((non0.beta != non0.betahat) & non0.betahat)
  # FN
  FN = sum((non0.beta != non0.betahat) & non0.beta)
  # return
  c(PE, EA1, EA2, EAInfty, FP, FN)
}
rownames(evals) = c("PE", "EA1", "EA2", "EAInfty", "FP", "FN")
eval.means = apply(evals, 1, mean)
eval.sds = apply(evals, 1, sd)
eval.ses = eval.sds / sqrt(numSims)
evals.df = data.frame("mean" = eval.means, "sd" = eval.sds, "se" = eval.ses)
evals.df
saveRDS(
  evals.df, 
  file = paste0(output_dir,
                "/slr_gicLC_simulations", 
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
