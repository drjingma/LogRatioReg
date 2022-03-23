# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 3/7/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_slrs"

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
  source("Kristyn/Functions/supervisedlogratios.R")
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
  sigma.settings = "latentVarModel"
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
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 0.25 # 1, 0.5, 0.25
  theta.value = 1 # weight on a1 -- 1
  a0 = 0 # 0
  
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
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # get latent variable
  U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
  # simulate y from latent variable
  y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
  # simulate X: 
  epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
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
  # plain slr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slr0approx = slr(x = X, y = Y, approx = TRUE)
  end.time = Sys.time()
  slr0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0approx_SBP = slr0approx$sbp
  
  slr0approx.coefs = coefficients(slr0approx$model)
  slr0approx.a0 = slr0approx.coefs[1]
  slr0approx.thetahat = slr0approx.coefs[-1]
  slr0approx.betahat = as.numeric(getBetaFromTheta(
    slr0approx.thetahat, sbp = slr0approx_SBP))

  # compute metrics on the selected model #
  slr0approx.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0approx_SBP),
    ilrX.test = getIlrX(X.test, sbp = slr0approx_SBP),
    n.train = n, n.test = n,
    thetahat0 = slr0approx.a0, thetahat = slr0approx.thetahat, betahat = slr0approx.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    slr0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0approx.thetahat != 0),
    "time" = slr0approx.timing
  ),
  paste0(output_dir, "/slr_approx_metrics", file.end))

  ##############################################################################
  # plain slr method (a balance regression method)
  #   -- spectral clustering (no approximation)
  ##############################################################################
  start.time = Sys.time()
  slr0 = slr(x = X, y = Y, approx = FALSE)
  end.time = Sys.time()
  slr0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0_SBP = slr0$sbp
  
  slr0.coefs = coefficients(slr0$model)
  slr0.a0 = slr0.coefs[1]
  slr0.thetahat = slr0.coefs[-1]
  slr0.betahat = as.numeric(getBetaFromTheta(
    slr0.thetahat, sbp = slr0_SBP))
  
  # compute metrics on the selected model #
  slr0.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0_SBP),
    ilrX.test = getIlrX(X.test, sbp = slr0_SBP),
    n.train = n, n.test = n,
    thetahat0 = slr0.a0, thetahat = slr0.thetahat, betahat = slr0.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0.thetahat != 0),
    "time" = slr0.timing
  ),
  paste0(output_dir, "/slr_metrics", file.end))
  
  ##############################################################################
  # mult slr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  #   -- use Rsq for model selection
  ##############################################################################
  start.time = Sys.time()
  slrmult0approx = slrmult(x = X, y = Y, max.clusters = 5, approx = TRUE)
  end.time = Sys.time()
  slrmult0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  which_fit_slrmult0approx = which.max(slrmult0approx$max.Rsqs)
  slrmult0approx_fit = slrmult0approx$models[[which_fit_slrmult0approx]]
  
  slrmult0approx_SBP = slrmult0approx_fit$sbp
  
  slrmult0approx.coefs = coefficients(slrmult0approx_fit$model)
  slrmult0approx.a0 = slrmult0approx.coefs[1]
  slrmult0approx.thetahat = slrmult0approx.coefs[-1]
  slrmult0approx.betahat = as.numeric(getBetaFromTheta(
    slrmult0approx.thetahat, sbp = slrmult0approx_SBP))
  
  # compute metrics on the selected model #
  slrmult0approx.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrmult0approx_SBP),
    ilrX.test = getIlrX(X.test, sbp = slrmult0approx_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrmult0approx.a0, thetahat = slrmult0approx.thetahat, betahat = slrmult0approx.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrmult0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrmult0approx.thetahat != 0),
    "time" = slrmult0approx.timing
  ),
  paste0(output_dir, "/slrmult_rsq_approx_metrics", file.end))
  
  ##############################################################################
  # mult slr method (a balance regression method)
  #   -- spectral clustering (no approximation)
  #   -- use Rsq for model selection
  ##############################################################################
  start.time = Sys.time()
  slrmult0 = slrmult(x = X, y = Y, max.clusters = 5, approx = FALSE)
  end.time = Sys.time()
  slrmult0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  which_fit_slrmult0 = which.max(slrmult0$max.Rsqs)
  slrmult0_fit = slrmult0$models[[which_fit_slrmult0]]
  
  slrmult0_SBP = slrmult0_fit$sbp
  
  slrmult0.coefs = coefficients(slrmult0_fit$model)
  slrmult0.a0 = slrmult0.coefs[1]
  slrmult0.thetahat = slrmult0.coefs[-1]
  slrmult0.betahat = as.numeric(getBetaFromTheta(
    slrmult0.thetahat, sbp = slrmult0_SBP))
  
  # compute metrics on the selected model #
  slrmult0.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrmult0_SBP),
    ilrX.test = getIlrX(X.test, sbp = slrmult0_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrmult0.a0, thetahat = slrmult0.thetahat, betahat = slrmult0.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrmult0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrmult0.thetahat != 0),
    "time" = slrmult0.timing
  ),
  paste0(output_dir, "/slrmult_rsq_metrics", file.end))
  
  ##############################################################################
  # hierarchical slr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  #   -- use Rsq for model selection
  ##############################################################################
  start.time = Sys.time()
  hslr0approx = hslr(x = X, y = Y, max.levels = 5, approx = TRUE)
  end.time = Sys.time()
  hslr0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  which_fit_hslr0approx = which.max(hslr0approx$max.Rsqs)
  hslr0approx_fit = hslr0approx$models[[which_fit_hslr0approx]]
  
  hslr0approx_SBP = hslr0approx_fit$sbp
  
  hslr0approx.coefs = coefficients(hslr0approx_fit$model)
  hslr0approx.a0 = hslr0approx.coefs[1]
  hslr0approx.thetahat = hslr0approx.coefs[-1]
  hslr0approx.betahat = as.numeric(getBetaFromTheta(
    hslr0approx.thetahat, sbp = hslr0approx_SBP))
  
  # compute metrics on the selected model #
  hslr0approx.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = hslr0approx_SBP),
    ilrX.test = getIlrX(X.test, sbp = hslr0approx_SBP),
    n.train = n, n.test = n,
    thetahat0 = hslr0approx.a0, thetahat = hslr0approx.thetahat, betahat = hslr0approx.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    hslr0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(hslr0approx.thetahat != 0),
    "time" = hslr0approx.timing
  ),
  paste0(output_dir, "/hslr_rsq_approx_metrics", file.end))
  
  ##############################################################################
  # hierarchical slr method (a balance regression method)
  #   -- spectral clustering (no approximation)
  #   -- use Rsq for model selection
  ##############################################################################
  start.time = Sys.time()
  hslr0 = hslr(x = X, y = Y, max.levels = 5, approx = FALSE)
  end.time = Sys.time()
  hslr0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  which_fit_hslr0 = which.max(hslr0$max.Rsqs)
  hslr0_fit = hslr0$models[[which_fit_hslr0]]
  
  hslr0_SBP = hslr0_fit$sbp
  
  hslr0.coefs = coefficients(hslr0_fit$model)
  hslr0.a0 = hslr0.coefs[1]
  hslr0.thetahat = hslr0.coefs[-1]
  hslr0.betahat = as.numeric(getBetaFromTheta(
    hslr0.thetahat, sbp = hslr0_SBP))
  
  # compute metrics on the selected model #
  hslr0.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = hslr0_SBP),
    ilrX.test = getIlrX(X.test, sbp = hslr0_SBP),
    n.train = n, n.test = n,
    thetahat0 = hslr0.a0, thetahat = hslr0.thetahat, betahat = hslr0.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    hslr0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(hslr0.thetahat != 0),
    "time" = hslr0.timing
  ),
  paste0(output_dir, "/hslr_rsq_metrics", file.end))
  
  ##############################################################################
  # mult slr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  #   -- use CV for model selection
  ##############################################################################
  start.time = Sys.time()
  slrmultcv0approx = cv.slr(x = X, y = Y, max.clusters = 5, nfolds = K, approx = TRUE)
  end.time = Sys.time()
  slrmultcv0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrmultcv0approx_fit = slrmultcv0approx$models[[slrmultcv0approx$nclusters_1se_idx]]
  
  slrmultcv0approx_SBP = slrmultcv0approx_fit$sbp
  
  slrmultcv0approx.coefs = coefficients(slrmultcv0approx_fit$model)
  slrmultcv0approx.a0 = slrmultcv0approx.coefs[1]
  slrmultcv0approx.thetahat = slrmultcv0approx.coefs[-1]
  slrmultcv0approx.betahat = as.numeric(getBetaFromTheta(
    slrmultcv0approx.thetahat, sbp = slrmultcv0approx_SBP))
  
  # compute metrics on the selected model #
  slrmultcv0approx.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrmultcv0approx_SBP),
    ilrX.test = getIlrX(X.test, sbp = slrmultcv0approx_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrmultcv0approx.a0, thetahat = slrmultcv0approx.thetahat, betahat = slrmultcv0approx.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrmultcv0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrmultcv0approx.thetahat != 0),
    "time" = slrmultcv0approx.timing
  ),
  paste0(output_dir, "/slrmult_cv_approx_metrics", file.end))
  
  ##############################################################################
  # mult slr method (a balance regression method)
  #   -- spectral clustering (no approximation)
  #   -- use CV for model selection
  ##############################################################################
  start.time = Sys.time()
  slrmultcv0 = cv.slr(x = X, y = Y, max.clusters = 5, nfolds = K, approx = FALSE)
  end.time = Sys.time()
  slrmultcv0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrmultcv0_fit = slrmultcv0$models[[slrmultcv0$nclusters_1se_idx]]
  
  slrmultcv0_SBP = slrmultcv0_fit$sbp
  
  slrmultcv0.coefs = coefficients(slrmultcv0_fit$model)
  slrmultcv0.a0 = slrmultcv0.coefs[1]
  slrmultcv0.thetahat = slrmultcv0.coefs[-1]
  slrmultcv0.betahat = as.numeric(getBetaFromTheta(
    slrmultcv0.thetahat, sbp = slrmultcv0_SBP))
  
  # compute metrics on the selected model #
  slrmultcv0.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrmultcv0_SBP),
    ilrX.test = getIlrX(X.test, sbp = slrmultcv0_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrmultcv0.a0, thetahat = slrmultcv0.thetahat, betahat = slrmultcv0.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrmultcv0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrmultcv0.thetahat != 0),
    "time" = slrmultcv0.timing
  ),
  paste0(output_dir, "/slrmult_cv_metrics", file.end))
  

  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


