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

# GIC settings
nlam = 200

# Simulation settings
numSims = 100
n = 100
p = 200
rho = 0.2
beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
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

# set.seed(1)
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
  complasso = ConstrLasso(
    Y, Z, Cmat = matrix(1, p, 1), 
    nlam = nlam, tol = tol)
  non0.betahats = (complasso$bet != 0) # diff lambda = diff col
  non0ct = apply(non0.betahats, 2, sum) # count of non-zero betas for each lambda
  which0ct.leq3sqrtn = which(non0ct <= 3 * sqrt(n)) # lambda lower bound criteria
  lambda = complasso$lambda[which0ct.leq3sqrtn]
  nlam = length(lambda)
  non0.betahats = non0.betahats[, which0ct.leq3sqrtn]
  non0ct = non0ct[which0ct.leq3sqrtn]
  betahats = complasso$bet[, which0ct.leq3sqrtn]
  a0s = complasso$int[which0ct.leq3sqrtn]
  GIC = rep(NA, nlam)
  PE = rep(NA, nlam)
  for(m in 1:nlam){
    a0.tmp = a0s[m]
    betahat.tmp = betahats[, m]
    sigmasq.hat = crossprod(Y - a0.tmp - Z %*% betahat.tmp) / n
    PE[m] = sigmasq.hat
    pvn = p
    s.lam = sum(non0.betahats[, m])
    GIC[m] = log(sigmasq.hat) + (s.lam - 1) * log(log(n)) * log(pvn) / n
  }
  lam.min.idx = which.min(GIC)
  lam.min = complasso$lambda[lam.min.idx]
  a0 = complasso$int[lam.min.idx]
  betahat = complasso$bet[, lam.min.idx]
  non0.betahat = non0.betahats[, lam.min.idx]
  
  # begin plots
  plot(lambda, GIC, type = "l")
  points(x = lam.min, GIC[lam.min.idx])
  text(x = lam.min, y = GIC[lam.min.idx], pos = 4,
       labels = paste0("(", round(lam.min, 3), ",", round(GIC[lam.min.idx], 3), ")"))
  # data.frame(lam = lambda, gic = GIC, non0s = non0ct)
  plot(lambda, non0ct, type = "l")
  # end plots
  
  # evaluate model
  # prediction error
  PE = PE[lam.min.idx]
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
rownames(bs.selected_variables) = colnames(X)

bs.selected_variables_numeric = apply(bs.selected_variables, 2, as.numeric)
bs.selection_percentages = apply(bs.selected_variables_numeric, 1, FUN = 
                                   function(x) sum(x, na.rm = TRUE) / bs.n)
names(bs.selection_percentages) = rownames(bs.selected_variables)
bs.results = list(
  seed = rng.seed,  
  selected_variables = bs.selected_variables, 
  selection_percentages = bs.selection_percentages
)

saveRDS(bs.results,
        file = paste0("Kristyn/Experiments/output",
                      "/complasso_selection", 
                      "_refit",
                      "_B", bs.n, 
                      "_seed", rng.seed,
                      ".rds"))

sort(bs.selection_percentages)


