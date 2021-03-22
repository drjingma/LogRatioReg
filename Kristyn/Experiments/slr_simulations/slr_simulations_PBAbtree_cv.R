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
nlam = 200
intercept = TRUE
K = 5
rho.type = "square"

# Simulation settings
numSims = 100
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
generate.theta = 2 # 1 = sparse beta, 2 = not-sparse beta
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

# set.seed(1)
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
  
  # simulate training data #
  # generate W
  W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V = exp(W)
  rowsumsV = apply(V, 1, sum)
  X = V / rowsumsV
  sbp = sbp.fromPBA(X) # contrasts matrix, a.k.a. sbp matrix
  U = getU(sbp = sbp) # U
  epsilon = rnorm(n, 0, sigma_eps)
  Xb = log(X) %*% U # ilr(X)
  # generate theta
  if(generate.theta == 1){ # theta that gives sparse beta
    theta = as.matrix(c(rep(0, p - 2), 1))
  } else if(generate.theta == 2){ # theta that gives not-sparse beta
    theta = as.matrix(c(1, rep(0, p - 2)))
  } else{ # ?
    stop("generate.theta isn't equal to 1 or 2")
  }
  beta = getBeta(theta, sbp = sbp)
  # generate Y
  Y = Xb %*% theta + epsilon
  
  # simulate test data #
  # simulate independent test set of size n
  # generate W
  W.test = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V.test = exp(W.test)
  rowsumsV.test = apply(V.test, 1, sum)
  X.test = V.test / rowsumsV.test
  sbp.test = sbp.fromPBA(X.test) # contrasts matrix, a.k.a. sbp matrix
  U.test = getU(sbp = sbp.test) # U
  epsilon.test = rnorm(n, 0, sigma_eps)
  Xb.test = log(X.test) %*% U.test # ilr(X)
  # generate Y
  Y.test = Xb.test %*% theta + epsilon.test
  
  # apply supervised log-ratios, using CV to select lambda
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
              rho.type = rho.type)
  btree = slr$btree
  
  # choose lambda
  lam.min.idx = which.min(slr$cvm)
  lam.min = slr$lambda[lam.min.idx]
  a0 = slr$int[lam.min.idx]
  thetahat = slr$bet[, lam.min.idx]
  betahat = U %*% thetahat
    
  # evaluate model #
  # 1. prediction error #
  # 1a. on training set #
  # get prediction error on training set
  Yhat.train = a0 + computeBalances(X, btree) %*% thetahat
  PE.train = crossprod(Y - Yhat.train) / n
  # 1b. on test set #
  # get prediction error on test set
  Yhat.test = a0 + computeBalances(X.test, btree) %*% thetahat
  PE.test = crossprod(Y.test - Yhat.test) / n
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  EA1 = sum(abs(betahat - beta))
  EA2 = sqrt(crossprod(betahat - beta))
  EAInfty = max(abs(betahat - beta))
  # 2b. estimation of theta
  # ...
  # 3. selection accuracy #
  # 3a. selection of beta #
  non0.beta = (beta != 0)
  non0s = sum(non0.beta)
  non0.betahat = (betahat != 0)
  # FP
  FP = sum((non0.beta != non0.betahat) & non0.betahat)
  # FN
  FN = sum((non0.beta != non0.betahat) & non0.beta)
  # TPR
  TPR = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
  # 3b. selection of theta
  # ...
  # return
  c(PE.train, PE.test, EA1, EA2, EAInfty, FP, FN, TPR)
}
rownames(evals) = c("PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR")

eval.means = apply(evals, 1, mean)
eval.sds = apply(evals, 1, sd)
eval.ses = eval.sds / sqrt(numSims)
evals.df = data.frame("mean" = eval.means, "sd" = eval.sds, "se" = eval.ses)
evals.df

saveRDS(
  evals, 
  file = paste0(output_dir,
                "/slr_cv_sims", 
                "_PBA", 
                "_theta", generate.theta,
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
saveRDS(
  evals.df, 
  file = paste0(output_dir,
                "/slr_cv_summaries", 
                 "_PBA", 
                "_theta", generate.theta,
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
