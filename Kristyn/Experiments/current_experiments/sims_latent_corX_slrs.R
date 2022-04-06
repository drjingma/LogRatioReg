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
  
  source("RCode/func_libs.R")
  source("Kristyn/Functions/slr.R")
  source("Kristyn/Functions/util.R")
  
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
  slrmax = 10
  sigma_eps1 = 0.1
  sigma_eps2 = 0.1
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 0.5 # 1, 0.5, 0.25
  theta.value = 1 # weight on a1 -- 1
  a0 = 0 # 0
  rho_alrXj = 0.2
  
  file.end = paste0(
    "_", sigma.settings,
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_dim", n, "x", p, 
    "_slrmax", slrmax,
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
  start.time = Sys.time()
  classo = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # cl.lam.idx = which.min(classo$cvm)
  oneSErule = min(classo$cvm) + (classo$cvsd[which.min(classo$cvm)] / K) * 1
  cl.lam.idx = which(classo$cvm <= oneSErule)[1]
  cl.a0 = classo$int[cl.lam.idx]
  cl.betahat = classo$bet[, cl.lam.idx]

  # compute metrics on the selected model #
  cl.metrics = getMetricsLLC(
    y.train = Y, y.test = Y.test,
    logX.train = log(X),
    logX.test = log(X.test),
    n.train = n, n.test = n,
    betahat0 = cl.a0, betahat = cl.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    cl.metrics,
    "betasparsity" = bspars,
    "logratios" = 0,
    "time" = cl.timing
  ),
  paste0(output_dir, "/classo_metrics", file.end))

  ##############################################################################
  # plain slr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slr0approx = slr(x = X, y = Y, approx = TRUE)
  end.time = Sys.time()
  slr0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slr0approx.coefs = getCoefsBM(
    coefs = coefficients(slr0approx$model), sbp = slr0approx$sbp)

  # compute metrics on the selected model #
  slr0approx.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0approx$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0approx$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0approx.coefs$a0, thetahat = slr0approx.coefs$bm.coefs,
    betahat = slr0approx.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    slr0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0approx.coefs$bm.coefs != 0),
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

  slr0.coefs = getCoefsBM(
    coefs = coefficients(slr0$model), sbp = slr0$sbp)

  # compute metrics on the selected model #
  slr0.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0.coefs$a0, thetahat = slr0.coefs$bm.coefs,
    betahat = slr0.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    slr0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0.coefs$bm.coefs != 0),
    "time" = slr0.timing
  ),
  paste0(output_dir, "/slr_metrics", file.end))
  
  ##############################################################################
  # slr1sc method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slr1sc0approx = slr1sc(x = X, y = Y, approx = TRUE)
  end.time = Sys.time()
  slr1sc0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1sc0approx.coefs = getCoefsBM(
    coefs = coefficients(slr1sc0approx$model), sbp = slr1sc0approx$sbp)
  
  # compute metrics on the selected model #
  slr1sc0approx.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1sc0approx$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1sc0approx$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1sc0approx.coefs$a0, thetahat = slr1sc0approx.coefs$bm.coefs,
    betahat = slr1sc0approx.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1sc0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1sc0approx.coefs$bm.coefs != 0),
    "time" = slr1sc0approx.timing
  ),
  paste0(output_dir, "/slr1sc_approx_metrics", file.end))
  
  ##############################################################################
  # slr1sc method (a balance regression method)
  #   -- spectral clustering (no approximation)
  ##############################################################################
  start.time = Sys.time()
  slr1sc0 = slr1sc(x = X, y = Y, approx = FALSE)
  end.time = Sys.time()
  slr1sc0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1sc0.coefs = getCoefsBM(
    coefs = coefficients(slr1sc0$model), sbp = slr1sc0$sbp)
  
  # compute metrics on the selected model #
  slr1sc0.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1sc0$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1sc0$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1sc0.coefs$a0, thetahat = slr1sc0.coefs$bm.coefs,
    betahat = slr1sc0.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1sc0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1sc0.coefs$bm.coefs != 0),
    "time" = slr1sc0.timing
  ),
  paste0(output_dir, "/slr1sc_metrics", file.end))

  ##############################################################################
  # cv.slr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  #   -- use CV to select T = # of clusters in 1st application of spectral
  #       clustering
  ##############################################################################
  start.time = Sys.time()
   slrcv0approx = cv.slr(
    x = X, y = Y, max.clusters = slrmax, nfolds = K, approx = TRUE)
  end.time = Sys.time()
   slrcv0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

   slrcv0approx_fit = slrcv0approx$models[[slrcv0approx$nclusters_1se_idx]]
   slrcv0approx.coefs = getCoefsBM(
     coefs = coefficients(slrcv0approx_fit$model), sbp = slrcv0approx_fit$sbp)

  # compute metrics on the selected model #
   slrcv0approx.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrcv0approx_fit$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrcv0approx_fit$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrcv0approx.coefs$a0, thetahat = slrcv0approx.coefs$bm.coefs,
    betahat = slrcv0approx.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
     slrcv0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrcv0approx.coefs$bm.coefs != 0),
    "time" = slrcv0approx.timing
  ),
  paste0(output_dir, "/slrcv_approx_metrics", file.end))

  ##############################################################################
  # cv.slr method (a balance regression method)
  #   -- spectral clustering (no approximation)
  #   -- use CV to select T = # of clusters in 1st application of spectral
  #       clustering
  ##############################################################################
  start.time = Sys.time()
  slrcv0 = cv.slr(
    x = X, y = Y, max.clusters = slrmax, nfolds = K, approx = FALSE)
  end.time = Sys.time()
   slrcv0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

   slrcv0_fit = slrcv0$models[[slrcv0$nclusters_1se_idx]]
   slrcv0.coefs = getCoefsBM(
    coefs = coefficients( slrcv0_fit$model), sbp = slrcv0_fit$sbp)

  # compute metrics on the selected model #
   slrcv0.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrcv0_fit$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrcv0_fit$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrcv0.coefs$a0, thetahat = slrcv0.coefs$bm.coefs,
    betahat = slrcv0.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
     slrcv0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrcv0.coefs$bm.coefs != 0),
    "time" = slrcv0.timing
  ),
  paste0(output_dir, "/slrcv_metrics", file.end))

  ##############################################################################
  # cv.hslr method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation)
  #   -- use CV to select T = # of levels i.e. hierarchical splits
  ##############################################################################
  start.time = Sys.time()
  hslrcv0approx = cv.hslr(
    x = X, y = Y, max.levels = slrmax, nfolds = K, approx = TRUE)
  end.time = Sys.time()
  hslrcv0approx.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  hslrcv0approx_fit = hslrcv0approx$models[[hslrcv0approx$nclusters_1se_idx]]
  hslrcv0approx.coefs = getCoefsBM(
    coefs = coefficients(hslrcv0approx_fit$model), sbp = hslrcv0approx_fit$sbp)

  # compute metrics on the selected model #
  hslrcv0approx.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = hslrcv0approx_fit$sbp),
    ilrX.test = getIlrX(X.test, sbp = hslrcv0approx_fit$sbp),
    n.train = n, n.test = n,
    thetahat0 = hslrcv0approx.coefs$a0, thetahat = hslrcv0approx.coefs$bm.coefs,
    betahat = hslrcv0approx.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    hslrcv0approx.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(hslrcv0approx.coefs$bm.coefs != 0),
    "time" = hslrcv0approx.timing
  ),
  paste0(output_dir, "/hslrcv_approx_metrics", file.end))

  ##############################################################################
  # cv.hslr method (a balance regression method)
  #   -- spectral clustering (no approximation)
  #   -- use CV to select T = # of clusters in 1st application of spectral
  #       clustering
  ##############################################################################
  start.time = Sys.time()
  hslrcv0 = cv.hslr(
    x = X, y = Y, max.levels = slrmax, nfolds = K, approx = FALSE)
  end.time = Sys.time()
  hslrcv0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  hslrcv0_fit = hslrcv0$models[[hslrcv0$nclusters_1se_idx]]
  hslrcv0.coefs = getCoefsBM(
    coefs = coefficients( hslrcv0_fit$model), sbp = hslrcv0_fit$sbp)

  # compute metrics on the selected model #
  hslrcv0.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = hslrcv0_fit$sbp),
    ilrX.test = getIlrX(X.test, sbp = hslrcv0_fit$sbp),
    n.train = n, n.test = n,
    thetahat0 = hslrcv0.coefs$a0, thetahat = hslrcv0.coefs$bm.coefs,
    betahat = hslrcv0.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    hslrcv0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(hslrcv0.coefs$bm.coefs != 0),
    "time" = hslrcv0.timing
  ),
  paste0(output_dir, "/hslrcv_metrics", file.end))
  
  ##############################################################################
  # cv.hslr1sc method (a balance regression method)
  #   -- spectral clustering (with rank 1 approximation), just once
  #   -- use CV to select T = # of levels i.e. hierarchical splits
  ##############################################################################
  start.time = Sys.time()
  hslr1sccvapprox = cv.hslr1sc(
    x = X, y = Y, max.levels = slrmax, nfolds = K, approx = TRUE)
  end.time = Sys.time()
  hslr1sccvapprox.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  hslr1sccvapprox_fit = hslr1sccvapprox$models[[hslr1sccvapprox$nclusters_1se_idx]]
  hslr1sccvapprox.coefs = getCoefsBM(
    coefs = coefficients(hslr1sccvapprox_fit$model), sbp = hslr1sccvapprox_fit$sbp)

  # compute metrics on the selected model #
  hslr1sccvapprox.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = hslr1sccvapprox_fit$sbp),
    ilrX.test = getIlrX(X.test, sbp = hslr1sccvapprox_fit$sbp),
    n.train = n, n.test = n,
    thetahat0 = hslr1sccvapprox.coefs$a0, thetahat = hslr1sccvapprox.coefs$bm.coefs,
    betahat = hslr1sccvapprox.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    hslr1sccvapprox.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(hslr1sccvapprox.coefs$bm.coefs != 0),
    "time" = hslr1sccvapprox.timing
  ),
  paste0(output_dir, "/hslr1sccv_approx_metrics", file.end))

  ##############################################################################
  # cv.hslr1sc method (a balance regression method)
  #   -- spectral clustering (no approximation), just once
  #   -- use CV to select T = # of clusters in 1st application of spectral
  #       clustering
  ##############################################################################
  start.time = Sys.time()
  hslr1sccv = cv.hslr1sc(
    x = X, y = Y, max.levels = slrmax, nfolds = K, approx = FALSE)
  end.time = Sys.time()
  hslr1sccv.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  hslr1sccv_fit = hslr1sccv$models[[hslr1sccv$nclusters_1se_idx]]
  hslr1sccv.coefs = getCoefsBM(
    coefs = coefficients( hslr1sccv_fit$model), sbp = hslr1sccv_fit$sbp)

  # compute metrics on the selected model #
  hslr1sccv.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = hslr1sccv_fit$sbp),
    ilrX.test = getIlrX(X.test, sbp = hslr1sccv_fit$sbp),
    n.train = n, n.test = n,
    thetahat0 = hslr1sccv.coefs$a0, thetahat = hslr1sccv.coefs$bm.coefs,
    betahat = hslr1sccv.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    hslr1sccv.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(hslr1sccv.coefs$bm.coefs != 0),
    "time" = hslr1sccv.timing
  ),
  paste0(output_dir, "/hslr1sccv_metrics", file.end))
  
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


