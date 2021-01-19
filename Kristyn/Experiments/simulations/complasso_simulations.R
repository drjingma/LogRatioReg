# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit compositional Lasso to the data
# Date: 1/2021
# Notes: 

getwd()

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)

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

# Method Settings
tol = 1e-4
nlam = 200

# Simulation settings
numSims = 100
n = 100
p = 200
rho = 0.5 # 0.2, 0.5
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
  .combine = cbind, 
  .noexport = c("ConstrLassoC0")
) %dopar% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  
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
  complasso0 = ConstrLasso(Y, Z, Cmat = matrix(1, p, 1), nlam = nlam, tol = tol)
  non0.betahats = (complasso0$bet != 0) # diff lambda = diff col
  non0ct = apply(non0.betahats, 2, sum) # count of non-zero betas for each lambda
  which0ct.leq3sqrtn = which(non0ct <= 3 * sqrt(n)) # lambda lower bound criteria
  lambda = seq(complasso0$lambda[1], 
               complasso0$lambda[max(which0ct.leq3sqrtn)], 
               length.out = nlam)
  complasso = ConstrLasso(Y, Z, Cmat = matrix(1, p, 1), lambda = lambda, 
                          nlam = nlam, tol = tol)
  nlam = length(lambda)
  non0.betahats = (complasso$bet != 0)
  non0ct = apply(non0.betahats, 2, sum)
  # which0ct.leq3sqrtn = which(non0ct <= 3 * sqrt(n))
  a0s = complasso$int
  betahats = complasso$bet
  GIC = rep(NA, nlam)
  for(m in 1:nlam){
    a0.tmp = a0s[m]
    betahat.tmp = betahats[, m]
    sigmasq.hat = crossprod(Y - a0.tmp - Z %*% betahat.tmp) / n
    pvn = p
    s.lam = sum(non0.betahats[, m])
    GIC[m] = log(sigmasq.hat) + (s.lam - 1) * log(log(n)) * log(pvn) / n
  }
  lam.min.idx = which.min(GIC)
  lam.min = complasso$lambda[lam.min.idx]
  a0 = complasso$int[lam.min.idx]
  betahat = complasso$bet[, lam.min.idx]
  non0.betahat = non0.betahats[, lam.min.idx]
  
  # plot(complasso$lambda, GIC, type = "l")
  
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
  PE = crossprod(Y.test - a0 - Z.test %*% betahat) / n
  # estimation accuracy
  EA1 = sum(abs(betahat - beta))
  EA2 = crossprod(betahat - beta)
  EAInfty = max(abs(betahat - beta))
  # FP
  FP = sum((non0.beta != non0.betahat) & non0.betahat)
  # FN
  FN = sum((non0.beta != non0.betahat) & non0.beta)
  # return
  c(PE, EA1, EA2, EAInfty, FP, FN)
}
rownames(evals) = c("PE", "EA1", "EA2", "EAInfty", "FP", "FN")
saveRDS(evals,
        file = paste0("Kristyn/Experiments/output",
                      "/complasso_simulations", 
                      "_dim", n, "x", p, 
                      "_rho", rho, 
                      "_seed", rng.seed,
                      ".rds"))
eval.means = apply(evals, 1, mean)
eval.sds = apply(evals, 1, sd)
eval.ses = eval.sds / sqrt(numSims)
evals.df = data.frame("mean" = eval.means, "sd" = eval.sds, "se" = eval.ses)
evals.df
saveRDS(evals.df,
        file = paste0("Kristyn/Experiments/output",
                      "/complasso_simulations_summary", 
                      "_dim", n, "x", p, 
                      "_rho", rho, 
                      "_seed", rng.seed,
                      ".rds"))
