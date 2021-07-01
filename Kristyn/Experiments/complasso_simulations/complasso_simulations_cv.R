# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit compositional Lasso to the data
# Date: 05/10/2021
# Notes: 

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output"

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
  .combine = cbind, 
  .noexport = c("ConstrLassoC0")
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  
  nlam = 100
  
  # simulate training data #
  # generate W
  W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V = exp(W)
  rowsumsV = apply(V, 1, sum)
  X = V / rowsumsV
  epsilon = rnorm(n, 0, sigma_eps)
  Z = log(X)
  # generate Y
  Y = Z %*% beta + epsilon
  
  # simulate test data #
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
  Y.test = Z.test %*% beta + epsilon.test
  
  # apply compositional lasso, using CV to select lambda
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = Z, Cmat = matrix(1, p, 1), nlam = nlam, 
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  
  # choose lambda
  lam.min.idx = which.min(complasso$cvm)
  lam.min = complasso$lambda[lam.min.idx]
  a0 = complasso$int[lam.min.idx]
  betahat = complasso$bet[, lam.min.idx]
  
  # evaluate model #
  
  # 1. prediction error #
  # 1a. on training set #
  # get prediction error on training set
  Yhat.train = a0 + log(X) %*% betahat
  PE.train = as.vector(crossprod(Y - Yhat.train) / n)
  # 1b. on test set #
  # get prediction error on test set
  Yhat.test = a0 + log(X.test) %*% betahat
  PE.test = as.vector(crossprod(Y.test - Yhat.test) / n)
  
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  EA1 = sum(abs(betahat - beta))
  EA2 = as.vector(sqrt(crossprod(betahat - beta)))
  EAInfty = max(abs(betahat - beta))
  
  # 3. selection accuracy #
  # 3a. selection of beta #
  
  # new version #
  non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  non0.betahat = abs(betahat) > 10e-8
  is0.betahat = betahat <= 10e-8
  # FP
  FP = sum(is0.beta & non0.betahat)
  # FN
  FN = sum((non0.beta != non0.betahat) & non0.beta)
  # TPR
  TPR = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
  # beta sparsity
  bspars = sum(non0.beta)
  
  # # old version #
  # non0.beta.old = (beta != 0)
  # is0.beta.old = (beta == 0)
  # non0.betahat.old = (betahat != 0)
  # is0.betahat.old = (betahat == 0)
  # # FP - old version
  # FP.old = sum(is0.beta.old & non0.betahat.old)
  # # FN - old version
  # FN.old = sum((non0.beta.old != non0.betahat.old) & non0.beta.old)
  # # TPR - old version
  # TPR.old = sum((non0.beta.old == non0.betahat.old) & non0.betahat.old) / 
  #   sum(non0.beta.old)
  # # beta sparsity - old version
  # bspars.old = sum(non0.beta.old)
  
  # return
  c(PE.train, PE.test, EA1, EA2, EAInfty, 
    FP, FN, TPR, bspars
    # FP.old, FN.old, TPR.old, bspars.old
    )
}
rownames(evals) = c("PEtr", "PEte", "EA1", "EA2", "EAInfty", 
                    "FP", "FN", "TPR", "betaSparsity"
                    # "FPold", "FNold", "TPRold", "betaSparsityOld"
                    )
eval.means = apply(evals, 1, mean)
eval.sds = apply(evals, 1, sd)
eval.ses = eval.sds / sqrt(numSims)
evals.df = data.frame("mean" = eval.means, "sd" = eval.sds, "se" = eval.ses)
evals.df

file.end = paste0(
  "_dim", n, "x", p, 
  "_", beta.settings, 
  "_rho", rho, 
  "_int", intercept,
  "_scale", scaling,
  "_K", K,
  "_seed", rng.seed,
  ".rds")

saveRDS(evals, file = paste0(output_dir, "/complasso_sims", file.end))
saveRDS(evals.df, file = paste0(output_dir, "/complasso_summaries", file.end))
