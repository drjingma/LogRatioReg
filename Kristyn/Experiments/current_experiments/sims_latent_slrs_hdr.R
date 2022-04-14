# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 4/13/2022

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
  # SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
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
  # slr method (a balance regression method)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0 = slr(
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
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
  # slr method (a balance regression method)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- TRUE, maximal
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0hdrmax = slr(
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, 
    highdegree.regularization.summary = "maximal",
    include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0hdrmax.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0hdrmax.coefs = getCoefsBM(
    coefs = coefficients(slr0hdrmax$model), sbp = slr0hdrmax$sbp)
  
  # compute metrics on the selected model #
  slr0hdrmax.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0hdrmax$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0hdrmax$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0hdrmax.coefs$a0, thetahat = slr0hdrmax.coefs$bm.coefs,
    betahat = slr0hdrmax.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0hdrmax.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0hdrmax.coefs$bm.coefs != 0),
    "time" = slr0hdrmax.timing
  ),
  paste0(output_dir, "/slr_hdr_original_metrics", file.end))
  
  ##############################################################################
  # slr method (a balance regression method)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- TRUE, mean
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0hdrmean = slr(
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, 
    highdegree.regularization.summary = "mean",
    include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0hdrmean.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0hdrmean.coefs = getCoefsBM(
    coefs = coefficients(slr0hdrmean$model), sbp = slr0hdrmean$sbp)
  
  # compute metrics on the selected model #
  slr0hdrmean.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0hdrmean$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0hdrmean$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0hdrmean.coefs$a0, thetahat = slr0hdrmean.coefs$bm.coefs,
    betahat = slr0hdrmean.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0hdrmean.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0hdrmean.coefs$bm.coefs != 0),
    "time" = slr0hdrmean.timing
  ),
  paste0(output_dir, "/slr_hdr_mean_metrics", file.end))
  
  ##############################################################################
  # slr method (a balance regression method)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- TRUE, mean
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0hdrmed = slr(
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, 
    highdegree.regularization.summary = "median",
    include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0hdrmed.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0hdrmed.coefs = getCoefsBM(
    coefs = coefficients(slr0hdrmed$model), sbp = slr0hdrmed$sbp)
  
  # compute metrics on the selected model #
  slr0hdrmed.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0hdrmed$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0hdrmed$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0hdrmed.coefs$a0, thetahat = slr0hdrmed.coefs$bm.coefs,
    betahat = slr0hdrmed.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0hdrmed.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0hdrmed.coefs$bm.coefs != 0),
    "time" = slr0hdrmed.timing
  ),
  paste0(output_dir, "/slr_hdr_median_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  
  ##############################################################################
  # slr method (a balance regression method)
  #   rank 1 approximation -- TRUE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0amap = slr(
    x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0amap.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0amap.coefs = getCoefsBM(
    coefs = coefficients(slr0amap$model), sbp = slr0amap$sbp)
  
  # compute metrics on the selected model #
  slr0amap.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0amap$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0amap$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0amap.coefs$a0, thetahat = slr0amap.coefs$bm.coefs,
    betahat = slr0amap.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0amap.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0amap.coefs$bm.coefs != 0),
    "time" = slr0amap.timing
  ),
  paste0(output_dir, "/slr_amini_approx_metrics", file.end))
  
  ##############################################################################
  # slr method (a balance regression method)
  #   rank 1 approximation -- TRUE
  #   amini regularization -- TRUE
  #   high degree regularization -- TRUE, maximal
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0hdrmaxamap = slr(
    x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
    highdegree.regularization = TRUE, 
    highdegree.regularization.summary = "maximal",
    include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0hdrmaxamap.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0hdrmaxamap.coefs = getCoefsBM(
    coefs = coefficients(slr0hdrmaxamap$model), sbp = slr0hdrmaxamap$sbp)
  
  # compute metrics on the selected model #
  slr0hdrmaxamap.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0hdrmaxamap$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0hdrmaxamap$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0hdrmaxamap.coefs$a0, thetahat = slr0hdrmaxamap.coefs$bm.coefs,
    betahat = slr0hdrmaxamap.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0hdrmaxamap.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0hdrmaxamap.coefs$bm.coefs != 0),
    "time" = slr0hdrmaxamap.timing
  ),
  paste0(output_dir, "/slr_hdr_original_amini_approx_metrics", file.end))
  
  ##############################################################################
  # slr method (a balance regression method)
  #   rank 1 approximation -- TRUE
  #   amini regularization -- TRUE
  #   high degree regularization -- TRUE, mean
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0hdrmeanamap = slr(
    x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
    highdegree.regularization = TRUE, 
    highdegree.regularization.summary = "mean",
    include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0hdrmeanamap.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0hdrmeanamap.coefs = getCoefsBM(
    coefs = coefficients(slr0hdrmeanamap$model), sbp = slr0hdrmeanamap$sbp)
  
  # compute metrics on the selected model #
  slr0hdrmeanamap.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0hdrmeanamap$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0hdrmeanamap$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0hdrmeanamap.coefs$a0, thetahat = slr0hdrmeanamap.coefs$bm.coefs,
    betahat = slr0hdrmeanamap.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0hdrmeanamap.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0hdrmeanamap.coefs$bm.coefs != 0),
    "time" = slr0hdrmeanamap.timing
  ),
  paste0(output_dir, "/slr_hdr_mean_amini_approx_metrics", file.end))
  
  ##############################################################################
  # slr method (a balance regression method)
  #   rank 1 approximation -- TRUE
  #   amini regularization -- TRUE
  #   high degree regularization -- TRUE, mean
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0hdrmedamap = slr(
    x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
    highdegree.regularization = TRUE, 
    highdegree.regularization.summary = "median",
    include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0hdrmedamap.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0hdrmedamap.coefs = getCoefsBM(
    coefs = coefficients(slr0hdrmedamap$model), sbp = slr0hdrmedamap$sbp)
  
  # compute metrics on the selected model #
  slr0hdrmedamap.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0hdrmedamap$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0hdrmedamap$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0hdrmedamap.coefs$a0, thetahat = slr0hdrmedamap.coefs$bm.coefs,
    betahat = slr0hdrmedamap.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0hdrmedamap.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0hdrmedamap.coefs$bm.coefs != 0),
    "time" = slr0hdrmedamap.timing
  ),
  paste0(output_dir, "/slr_hdr_median_amini_approx_metrics", file.end))
  
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}

