# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 3/7/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs"

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores() / 2
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Other simulation settings
numSims = 100

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSims
) %dorng% {
  # rm(list=ls())
  library(mvtnorm)
  
  library(Matrix)
  library(glmnet)
  
  library(balance)
  library(propr)
  
  source("RCode/func_libs.R")
  source("Kristyn/Functions/util.R")
  source("Kristyn/Functions/supervisedlogratioseta.R")
  source("Kristyn/Functions/HSClust.R")
  source("Kristyn/Functions/slrnew.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  source("Kristyn/Functions/helper_functions.R")
  
  # for plots
  library(ggraph) # make dendrogram
  library(igraph) # transform dataframe to graph object: graph_from_data_frame()
  library(tidygraph)
  
  # Tuning parameters###########################################################
  
  # Settings to toggle with
  sigma.settings = "latentVarModel_corX"
  n = 100
  p = 30
  K = 10
  nlam = 100
  neta = p
  intercept = TRUE
  scaling = TRUE
  tol = 1e-4
  sigma_eps1 = 0.1
  sigma_eps2 = 0.1
  # SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  SBP.true = matrix(c(1, -1, -1, -1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 1 # 1, 0.5, 0.25
  theta.value = 1 # weight on a1 -- 1
  a0 = 0 # 0
  rho_alrXj = 0.2
  
  file.end = paste0(
    "_", sigma.settings,
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_dim", n, "x", p, 
    "_noisey", sigma_eps1, 
    "_noisex", sigma_eps2, 
    "_b0", b0, 
    "_b1", b1, 
    "_a0", a0, 
    "_theta", theta.value,
    "_rho", rho_alrXj, 
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # get latent variable
  U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
  # simulate y from latent variable
  y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
  # simulate X: 
  epsj.all = mvrnorm(
    n = 2 * n, mu = rep(0, p - 1), 
    Sigma = sigma_eps2 * rgExpDecay(p - 1, rho_alrXj)$Sigma)
  # epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
  a1 = theta.value * ilrtrans.true$ilr.trans[-p] 
  #   alpha1j = {
  #     c1=theta*ilr.const/k+   if j \in I+
  #     -c2=-theta*ilr.const/k-  if j \in I-
  #     0                       o/w
  #   }
  alrXj.all = a0 + U.all %*% t(a1) + epsj.all #log(Xj/Xp) =alpha0j+alpha1j*U+epsj
  X.all <- alrinv(alrXj.all)
  colnames(X.all) = paste0('s', 1:p)
  
  # subset out training and test sets
  X = X.all[1:n, ]
  X.test = X.all[-(1:n), ]
  Y <- y.all[1:n]
  Y.test <- y.all[-(1:n)]
  
  # about beta
  non0.beta = as.vector(SBP.true != 0)
  is0.beta = !non0.beta
  bspars = sum(non0.beta)
  # solve for beta
  c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
  beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
    as.vector(ilrtrans.true$ilr.trans)
  
  ##############################################################################
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  # ........
}


