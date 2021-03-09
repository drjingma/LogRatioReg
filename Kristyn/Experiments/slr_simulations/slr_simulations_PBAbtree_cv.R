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

# Simulation settings
numSims = 100
n = 50
p = 30
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
  U = sbp.fromPBA(X) # transformation matrix, U
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
  beta = U %*% as.matrix(theta) # beta = U' theta
  # generate Y
  Y = Xb %*% theta + epsilon
  
  # apply supervised log-ratios, using CV to select lambda
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept)
  # transformed
  btree = slr$btree
  
  lam.min.idx = which.min(slr$cvm)
  lam.min = slr$lambda[lam.min.idx]
  a0 = slr$int[lam.min.idx]
  thetahat = slr$bet[, lam.min.idx]
  betahat = U %*% thetahat
    
  # evaluate model
  # prediction error
  # simulate independent test set of size n
  # generate W
  W.test = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V.test = exp(W.test)
  rowsumsV.test = apply(V.test, 1, sum)
  X.test = V.test / rowsumsV.test
  U.test = sbp.fromPBA(X.test) # transformation matrix, U
  epsilon.test = rnorm(n, 0, sigma_eps)
  Xb.test = log(X.test) %*% U # ilr(X)
  # generate Y
  Y.test = Xb.test %*% theta + epsilon.test
  # get fitted model
  predictSLR = function(X){
    a0 + computeBalances(X, btree) %*% thetahat
  }
  Y.pred = predictSLR(X.test)
  PE = crossprod(Y.test - Y.pred) / n
  # # estimation accuracy
  EA1 = sum(abs(betahat - beta))
  EA2 = crossprod(betahat - beta)
  EAInfty = max(abs(betahat - beta))
  non0.beta = (beta != 0)
  non0s = sum(non0.beta)
  print(paste0("number of non-zero elements of beta: ", non0s, " out of ", p))
  non0.betahat = (betahat != 0)
  # # FP
  FP = sum((non0.beta != non0.betahat) & non0.betahat)
  # # FN
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
# # when intercept = TRUE, generate.theta = 1
# mean          sd          se
# PE       0.31496494  0.09164212 0.009164212
# EA1      0.59887596  0.43609169 0.043609169
# EA2      0.07436537  0.07264030 0.007264030
# EAInfty  0.15682232  0.07494411 0.007494411
# FP      10.59000000 10.28787165 1.028787165
# FN       0.00000000  0.00000000 0.000000000
# # when intercept = TRUE, generate.theta = 2
# mean         sd         se
# PE      13.770938 11.7608132 1.17608132
# EA1     34.813792  5.0121014 0.50121014
# EA2     62.003835 17.6985894 1.76985894
# EAInfty  3.256607  0.6599112 0.06599112
# FP       0.000000  0.0000000 0.00000000
# FN       0.160000  0.8005049 0.08005049
saveRDS(
  evals.df, 
  file = paste0(output_dir,
                "/slr_cv_simulations", 
                 "_PBA", 
                "_theta", generate.theta,
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
