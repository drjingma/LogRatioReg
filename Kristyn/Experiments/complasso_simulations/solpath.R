# last updated: 05/11/2021
# simulate a data set using settings in Lin et al 2014, 
# fit Complasso and SLR, and plot the solution path.

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
nlam = 3 # for getting min and max lambda, first
intercept = TRUE
K = 10
n = 100
p = 200
rho = 0.5 # 0.2, 0.5

# Simulation settings
numSims = 50#100
# which beta?
beta.settings = "new"
if(beta.settings == "old" | beta.settings == "linetal2014"){
  beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
} else{
  beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
}
names(beta) = paste0('s',1:p)
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
# Sigma0 = rgExpDecay(p,rho)$Sigma
# all.equal(SigmaW, Sigma0)

################################################################################
# Simulations -- different lambda sequence #
################################################################################

registerDoRNG(rng.seed)
sims = foreach(
  b = 1:numSims
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
  # simulate training data #
  # generate W
  W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  colnames(W) = paste0('s',1:p)
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V = exp(W)
  rowsumsV = apply(V, 1, sum)
  X = V / rowsumsV
  epsilon = rnorm(n, 0, sigma_eps)
  Z = log(X)
  # generate Y
  Y = Z %*% beta + epsilon
  
  # apply compositional lasso
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), 
    nlam = nlam, nfolds = K, tol = tol, intercept = intercept)
  
  # apply supervised log-ratios
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept,
              rho.type = rho.type, linkage = linkage)
  
  # apply slralpha=0.5
  slr0.5 = cvSLRalpha(
    y = Y, X = X, nlam = nlam, nfolds = K, alpha = 0.5, 
    intercept = intercept, rho.type = rho.type, linkage = linkage)
  
  # apply slralpha=1
  slr1 = cvSLRalpha(
    y = Y, X = X, nlam = nlam, nfolds = K, alpha = 1, 
    intercept = intercept, rho.type = rho.type, linkage = linkage)
  
  list(X = X, Y = Y,
       fit.cl = complasso, fit.slr = slr, fit.slr0.5 = slr0.5, fit.slr1 = slr1
  )
}



# organize the simulation results to have more easily-accessible matrices #
# 
lambda.cl.mat = matrix(NA, nlam, numSims)
lambda.slr.mat = matrix(NA, nlam, numSims)
lambda.slr0.5.mat = matrix(NA, nlam, numSims)
lambda.slr1.mat = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  sim.tmp = sims[[i]]
  lambda.cl.mat[, i] = sim.tmp$fit.cl$lambda
  lambda.slr.mat[, i] = sim.tmp$fit.slr$lambda
  lambda.slr0.5.mat[, i] = sim.tmp$fit.slr0.5$lambda
  lambda.slr1.mat[, i] = sim.tmp$fit.slr1$lambda
}



################################################################################
# using the same lambda sequence in each method
nlam = 100#200
lambda.df = data.frame(
  complasso = exp(seq(max(log(lambda.cl.mat)), 
                      min(log(lambda.cl.mat)), length.out = nlam)),
  slr = exp(seq(max(log(lambda.slr.mat)), 
                min(log(lambda.slr.mat)),length.out = nlam)), 
  slr0.5 = exp(seq(max(log(lambda.slr0.5.mat)), 
                min(log(lambda.slr0.5.mat)),length.out = nlam)), 
  slr1 = exp(seq(max(log(lambda.slr1.mat)), 
                min(log(lambda.slr1.mat)),length.out = nlam))
)


registerDoRNG(rng.seed)
sims3 = foreach(
  b = 1:numSims
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
  print(paste0("starting sim ", b))
  
  # get simulated data #
  X = sims[[b]]$X
  Y = sims[[b]]$Y
  
  # apply compositional lasso
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), 
    lambda = lambda.df$complasso, nfolds = K, tol = tol, intercept = intercept)
  # calculate TPR
  cl.res = apply(complasso$bet, 2, function(a) tpr.for.coef(beta, a))
  # cl.res = apply(complasso$bet, 2, function(a) 
  #   getTPR( type = "llc", beta = beta, betahat = a))
  S.hat.cl = cl.res[1, ]
  TPR.cl = cl.res[2, ]
  
  # apply supervised log-ratios
  slr = cvSLR(
    y = Y, X = X, lambda = lambda.df$slr, nfolds = K, intercept = intercept,
    rho.type = rho.type, linkage = linkage)
  # caculate TPR
  btree.slr = slr$btree
  SBP = sbp.fromHclust(btree.slr)
  slr.res = apply(slr$bet, 2, function(a) tpr.for.coef.ilr(beta, a, SBP))
  # slr.res = apply(slr$bet, 2, function(a) 
  #   getTPR( type = "ilr", beta = beta, thetahat = a, sbp = SBP))
  S.hat.slr = slr.res[1, ]
  TPR.slr = slr.res[2, ]
  
  # apply slr0.5
  # apply slralpha=0.5
  slr0.5 = cvSLRalpha(
    y = Y, X = X, lambda = lambda.df$slr0.5, nfolds = K, alpha = 0.5, 
    intercept = intercept, rho.type = rho.type, linkage = linkage)
  # caculate TPR
  btree.slr0.5 = slr0.5$btree
  SBP0.5 = sbp.fromHclust(btree.slr0.5)
  slr0.5.res = apply(slr0.5$bet, 2, function(a) tpr.for.coef.ilr(beta, a, SBP0.5))
  S.hat.slr0.5 = slr0.5.res[1, ]
  TPR.slr0.5 = slr0.5.res[2, ]
  
  # apply slralpha=1
  slr1 = cvSLRalpha(
    y = Y, X = X, lambda = lambda.df$slr1, nfolds = K, alpha = 1, 
    intercept = intercept, rho.type = rho.type, linkage = linkage)
  # caculate TPR
  btree.slr1 = slr1$btree
  SBP1 = sbp.fromHclust(btree.slr1)
  slr1.res = apply(slr1$bet, 2, function(a) tpr.for.coef.ilr(beta, a, SBP1))
  S.hat.slr1 = slr1.res[1, ]
  TPR.slr1 = slr1.res[2, ]
  
  list(X = X, Y = Y, 
       fit.cl = complasso, TPR.cl = TPR.cl, S.hat.cl = S.hat.cl, 
       fit.slr = slr, TPR.slr = TPR.slr, S.hat.slr = S.hat.slr, 
       fit.slr0.5 = slr0.5, TPR.slr0.5 = TPR.slr0.5, S.hat.slr0.5 = S.hat.slr0.5, 
       fit.slr1 = slr1, TPR.slr1 = TPR.slr1, S.hat.slr1 = S.hat.slr1
       )
}

# save results #################################################################

file.end = paste0(
  "_dim", n, "x", p,
  "_", beta.settings,
  "_rho", rho,
  "_type", rho.type,
  "_int", intercept,
  "_K", K,
  "_numSims", numSims,
  "_seed", rng.seed,
  ".rds")

saveRDS(sims3, file = paste0(output_dir, "/solpaths", file.end))