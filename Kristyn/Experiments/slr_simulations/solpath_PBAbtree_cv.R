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
intercept = TRUE
K = 5
rho.type = "square"

# Simulation settings
numSims = 100
n = 100
p = 1000
rho = 0.2 # 0.2, 0.5
generate.theta = 1 # 1 = sparse beta, 2 = not-sparse beta
sigma_eps = 0.2
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

set.seed(1)

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

# apply compositional lasso
complasso = cv.func(
  method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
  nfolds = K, tol = tol, intercept = intercept)
################################################################################
saveRDS(
  complasso, 
  file = paste0(output_dir,
                "/complasso_1sim", 
                "_PBA", 
                "_theta", generate.theta,
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
################################################################################
# calculate TPR for compositional lasso
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
# plot(x = complasso$lambda, y = TPR.cl, type = "l")
# plot(x = S.hat.cl, y = TPR.cl, type = "l")

# apply supervised log-ratios
slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
            rho.type = rho.type)
################################################################################
saveRDS(
  slr, 
  file = paste0(output_dir,
                "/slr_1sim", 
                "_PBA", 
                "_theta", generate.theta,
                "_dim", n, "x", p, 
                "_rho", rho, 
                "_int", intercept,
                "_seed", rng.seed,
                ".rds"))
################################################################################
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
# plot(x = slr$lambda, y = TPR.slr, type = "l")
# plot(x = S.hat.slr, y = TPR.slr, type = "l")

# library(ggplot2)
# slr.gg = data.frame(
#   lambda = slr$lambda, 
#   S.hat = S.hat.slr, 
#   TPR = TPR.slr, 
#   ID = "SLR")
# slr.gg = slr.gg[order(slr.gg$S.hat), ]
# cl.gg = data.frame(
#   lambda = complasso$lambda, 
#   S.hat = S.hat.cl, 
#   TPR = TPR.cl, 
#   ID = "CompLasso")
# cl.gg = cl.gg[order(cl.gg$S.hat), ]
# data.gg = rbind(slr.gg, cl.gg)
# ggplot(data.gg, aes(x = S.hat, y = TPR)) + 
#   geom_path(aes(color = ID, linetype = ID), size = 1) + 
#   theme_bw()
