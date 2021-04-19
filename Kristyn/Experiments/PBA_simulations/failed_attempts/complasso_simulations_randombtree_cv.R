# attempt at balance tree simulations, for the supervised log-ratios method

getwd()
output_dir = "Kristyn/Experiments/slr_simulations/output"

# libraries
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
tol = 1e-4
nlam = 200
intercept = FALSE
K = 5

# Simulation settings
numSims = 100
n = 50
p = 30
rho = 0.2 # 0.2, 0.5
theta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 1 - 8))
non0.theta = (theta != 0)
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
  U = sbp.fromRandom(log(X))
  Xb = log(X) %*% U
  beta = U %*% as.matrix(theta) # beta = U' theta
  # generate Y
  Y = Xb %*% matrix(theta) + epsilon
  
  # apply compositional lasso, using CV to select lambda
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
    nfolds = K, tol = tol, intercept = intercept)
  
  # choose lambda
  lam.min.idx = which.min(complasso$cvm)
  lam.min = complasso$lambda[lam.min.idx]
  a0 = complasso$int[lam.min.idx]
  betahat = complasso$bet[, lam.min.idx]
  
  # plot(complasso$lambda, complasso$cvm, type = "l")
  
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
  U.test = sbp.fromRandom(log(X.test))
  Xb.test = log(X.test) %*% U.test
  # generate Y
  Y.test = Xb.test %*% matrix(theta) + epsilon.test
  # get prediction error
  Y.pred = a0 + log(X.test) %*% betahat
  PE = crossprod(Y.test - Y.pred) / n
  # estimation accuracy
  non0.beta = (beta != 0)
  non0.betahat = (betahat != 0)
  EA1 = sum(abs(betahat - beta))
  EA2 = crossprod(betahat - beta)
  EAInfty = max(abs(betahat - beta))
  non0.betahat = (betahat != 0)
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
# # when intercept = TRUE
# mean           sd           se
# PE      9774.894517 1.797985e+04 1.797985e+03
# EA1       24.352753 4.235556e+00 4.235556e-01
# EA2       41.804803 1.417365e+01 1.417365e+00
# EAInfty    3.050527 8.050459e-01 8.050459e-02
# FP        18.780000 4.303111e+00 4.303111e-01
# FN         1.340000 1.319780e+00 1.319780e-01
# # when intercept = FALSE
# mean           sd           se
# PE      9726.672253 1.788641e+04 1.788641e+03
# EA1       19.411885 6.414873e+00 6.414873e-01
# EA2       30.614037 1.501114e+01 1.501114e+00
# EAInfty    2.654308 7.686347e-01 7.686347e-02
# FP        14.570000 6.396424e+00 6.396424e-01
# FN         2.270000 1.600852e+00 1.600852e-01
saveRDS(
  evals.df, 
  file = paste0(output_dir,
                "/slr_cv_simulations", 
                # "_", d.method, 
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
