# last updated: 04/11/2021
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

# helper functions
standardizeXY <- function(X, Y){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # Center Y
  Ymean = mean(Y)
  Ytilde = Y - mean(Y)
  
  # Center and scale X
  Xmeans = colMeans(X)
  Xcen = X - matrix(Xmeans, n, p, byrow=T)
  normsX2 = colSums(Xcen^2) / n
  weights = 1 / sqrt(normsX2) # should weights be the vector, or the matix?
  Xtilde = Xcen %*% diag(weights)
  
  # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# Method Settings
tol = 1e-4
nlam = 3 # for testing
intercept = TRUE
K = 5
rho.type = "square"

# Simulation settings
numSims = 100
n = 100
p = 200
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
# if simulations are saved, read them in

# sims3 = readRDS(paste0(output_dir,
#                        "/sims_solpaths",
#                        "_theta_", paste(indices.theta, collapse = "_"),
#                        "_dim", n, "x", p,
#                        "_rho", rho,
#                        "_int", intercept,
#                        "_seed", rng.seed,
#                        "_numSims", numSims,
#                        ".rds"))

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
  
  TPR.cl = rep(NA, nlam)
  S.hat.cl = rep(NA, nlam)
  for(l in 1:nlam){
    a0 = complasso$int[l]
    betahat = complasso$bet[, l]
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR.cl[l] = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat.cl[l] = sum(non0.betahat)
  }
  
  # apply supervised log-ratios
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
              rho.type = rho.type)
  
  btree.slr = slr$btree
  # calculate TPR for supervised log-ratios
  nlam.slr = length(slr$lambda)
  TPR.slr = rep(NA, nlam.slr)
  S.hat.slr = rep(NA, nlam.slr)
  for(l in 1:nlam.slr){
    a0 = slr$int[l]
    thetahat = slr$bet[, l]
    betahat = getBeta(thetahat, btree.slr)
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR.slr[l] = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat.slr[l] = sum(non0.betahat)
  }
  
  list(X = X, Y = Y,
       fit.cl = complasso, TPR.cl = TPR.cl, S.hat.cl = S.hat.cl, 
       fit.slr = slr, TPR.slr = TPR.slr, S.hat.slr = S.hat.slr)
}



# organize the simulation results to have more easily-accessable matrices #
# cl
fit.cl.list = list()
TPR.cl.mat = matrix(NA, nlam, numSims)
S.hat.cl.mat = matrix(NA, nlam, numSims)
lambda.cl.mat = matrix(NA, nlam, numSims)
# slr
fit.slr.list = list()
TPR.slr.mat = matrix(NA, nlam, numSims)
S.hat.slr.mat = matrix(NA, nlam, numSims)
lambda.slr.mat = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  sim.tmp = sims[[i]]
  # cl
  fit.cl.list[[i]] = sim.tmp$fit.cl
  TPR.cl.mat[, i] = sim.tmp$TPR.cl
  S.hat.cl.mat[, i] = sim.tmp$S.hat.cl
  lambda.cl.mat[, i] = sim.tmp$fit.cl$lambda
  # slr
  fit.slr.list[[i]] = sim.tmp$fit.slr
  TPR.slr.mat[, i] = sim.tmp$TPR.slr
  S.hat.slr.mat[, i] = sim.tmp$S.hat.slr
  lambda.slr.mat[, i] = sim.tmp$fit.slr$lambda
}



################################################################################
# using the same lambda sequence in each method
nlam = 200
lambda.df = data.frame(
  complasso = exp(seq(max(log(lambda.cl.mat)), 
                      min(log(lambda.cl.mat)), length.out = nlam)),
  slr = exp(seq(max(log(lambda.slr.mat)), 
                min(log(lambda.slr.mat)),length.out = nlam))
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
  
  TPR.cl = rep(NA, nlam)
  S.hat.cl = rep(NA, nlam)
  for(l in 1:nlam){
    a0 = complasso$int[l]
    betahat = complasso$bet[, l]
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR.cl[l] = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat.cl[l] = sum(non0.betahat)
  }
  
  # apply supervised log-ratios
  slr = cvSLR(
    y = Y, X = X, lambda = lambda.df$slr, nfolds = K, intercept = intercept, 
    rho.type = rho.type)
  
  btree.slr = slr$btree
  # calculate TPR for supervised log-ratios
  nlam.slr = length(slr$lambda)
  TPR.slr = rep(NA, nlam.slr)
  S.hat.slr = rep(NA, nlam.slr)
  for(l in 1:nlam.slr){
    a0 = slr$int[l]
    thetahat = slr$bet[, l]
    betahat = getBeta(thetahat, btree.slr)
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR.slr[l] = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat.slr[l] = sum(non0.betahat)
  }
  
  list(X = X, Y = Y, 
       fit.cl = complasso, TPR.cl = TPR.cl, S.hat.cl = S.hat.cl, 
       fit.slr = slr, TPR.slr = TPR.slr, S.hat.slr = S.hat.slr)
}

# organize the simulation results to have more easily-accessable matrices #
# cl
fit.cl.list = list()
TPR.cl.mat = matrix(NA, nlam, numSims)
S.hat.cl.mat = matrix(NA, nlam, numSims)
lambda.cl.mat = matrix(NA, nlam, numSims)
# slr
fit.slr.list = list()
TPR.slr.mat = matrix(NA, nlam, numSims)
S.hat.slr.mat = matrix(NA, nlam, numSims)
lambda.slr.mat = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  sim.tmp = sims3[[i]]
  # cl
  fit.cl.list[[i]] = sim.tmp$fit.cl
  TPR.cl.mat[, i] = sim.tmp$TPR.cl
  S.hat.cl.mat[, i] = sim.tmp$S.hat.cl
  lambda.cl.mat[, i] = sim.tmp$fit.cl$lambda
  # slr
  fit.slr.list[[i]] = sim.tmp$fit.slr
  TPR.slr.mat[, i] = sim.tmp$TPR.slr
  S.hat.slr.mat[, i] = sim.tmp$S.hat.slr
  lambda.slr.mat[, i] = sim.tmp$fit.slr$lambda
}



# average TPR and S.hat ########################################################

dim(S.hat.cl.mat)
S.hat.cl.avg = apply(S.hat.cl.mat, 1, mean, na.rm = TRUE)
TPR.cl.avg = apply(TPR.cl.mat, 1, mean, na.rm = TRUE)
S.hat.slr.avg = apply(S.hat.slr.mat, 1, mean, na.rm = TRUE)
TPR.slr.avg = apply(TPR.slr.mat, 1, mean, na.rm = TRUE)
S.hat.or.avg = apply(S.hat.or.mat, 1, mean, na.rm = TRUE)
TPR.or.avg = apply(TPR.or.mat, 1, mean, na.rm = TRUE)

# complasso stuff
cl.gg.complete = data.frame(
  "S.hat" = S.hat.cl.avg, 
  "TPR" = TPR.cl.avg)
cl.gg.complete$Type = "CompLasso"
# slr stuff
slr.gg.complete = data.frame(
  "S.hat" = S.hat.slr.avg, 
  "TPR" = TPR.slr.avg)
slr.gg.complete$Type = "SLR"
# or.stuff
or.gg.complete = data.frame(
  "S.hat" = S.hat.or.avg, 
  "TPR" = TPR.or.avg)
or.gg.complete$Type = "Oracle"
# ggplot
gg.complete = rbind(slr.gg.complete, cl.gg.complete, or.gg.complete)
gg.complete$Type = factor(gg.complete$Type, levels = c("CompLasso", "SLR", "Oracle"))
ggplot(gg.complete, aes(x = S.hat, y = TPR)) + 
  geom_line(aes(color = Type), size = 1) +
  theme_bw()


# save results #################################################################

saveRDS(
  sims3,
  file = paste0(output_dir,
                "/sims_solpaths",
                "_dim", n, "x", p,
                "_rho", rho,
                "_int", intercept,
                "_seed", rng.seed,
                "_numSims", numSims,
                ".rds"))
