# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 4/27/2022

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
  source("Kristyn/Functions/slr1sc.R")
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
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  # SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 1 # 1, 0.5, 0.25
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
  
  cl.lam.min.idx = which.min(classo$cvm)
  cl.a0 = classo$int[cl.lam.min.idx]
  cl.betahat = classo$bet[, cl.lam.min.idx]
  
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
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   spectral.clustering.algorithm == "kmeans"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrk = slr(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "kmeans")
  end.time = Sys.time()
  slrk.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrk.coefs = getCoefsBM(
    coefs = coefficients(slrk$model), sbp = slrk$sbp)
  
  # compute metrics on the selected model #
  slrk.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrk$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrk$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrk.coefs$a0, thetahat = slrk.coefs$bm.coefs,
    betahat = slrk.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrk.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrk.coefs$bm.coefs != 0),
    "time" = slrk.timing
  ),
  paste0(output_dir, "/slr_kmeans_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   spectral.clustering.algorithm == "kmeans"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrk = slr(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "kmeans")
  end.time = Sys.time()
  slrk.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrk.coefs = getCoefsBM(
    coefs = coefficients(slrk$model), sbp = slrk$sbp)
  
  # compute metrics on the selected model #
  slrk.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrk$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrk$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrk.coefs$a0, thetahat = slrk.coefs$bm.coefs,
    betahat = slrk.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrk.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrk.coefs$bm.coefs != 0),
    "time" = slrk.timing
  ),
  paste0(output_dir, "/slr_kmeans_amini_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   spectral.clustering.algorithm == "kmeans"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- TRUE
  ##############################################################################
  start.time = Sys.time()
  slrk = slr(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, spectral.clustering.method = "kmeans")
  end.time = Sys.time()
  slrk.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrk.coefs = getCoefsBM(
    coefs = coefficients(slrk$model), sbp = slrk$sbp)
  
  # compute metrics on the selected model #
  slrk.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrk$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrk$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrk.coefs$a0, thetahat = slrk.coefs$bm.coefs,
    betahat = slrk.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrk.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrk.coefs$bm.coefs != 0),
    "time" = slrk.timing
  ),
  paste0(output_dir, "/slr_kmeans_hdr_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   spectral.clustering.algorithm == "cut"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrc = slr(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "cut")
  end.time = Sys.time()
  slrc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrc.coefs = getCoefsBM(
    coefs = coefficients(slrc$model), sbp = slrc$sbp)
  
  # compute metrics on the selected model #
  slrc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrc.coefs$a0, thetahat = slrc.coefs$bm.coefs,
    betahat = slrc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrc.coefs$bm.coefs != 0),
    "time" = slrc.timing
  ),
  paste0(output_dir, "/slr_cut_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   spectral.clustering.algorithm == "cut"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrc = slr(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "cut")
  end.time = Sys.time()
  slrc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrc.coefs = getCoefsBM(
    coefs = coefficients(slrc$model), sbp = slrc$sbp)
  
  # compute metrics on the selected model #
  slrc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrc.coefs$a0, thetahat = slrc.coefs$bm.coefs,
    betahat = slrc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrc.coefs$bm.coefs != 0),
    "time" = slrc.timing
  ),
  paste0(output_dir, "/slr_cut_amini_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   spectral.clustering.algorithm == "cut"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrc = slr(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, spectral.clustering.method = "cut")
  end.time = Sys.time()
  slrc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrc.coefs = getCoefsBM(
    coefs = coefficients(slrc$model), sbp = slrc$sbp)
  
  # compute metrics on the selected model #
  slrc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrc.coefs$a0, thetahat = slrc.coefs$bm.coefs,
    betahat = slrc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrc.coefs$bm.coefs != 0),
    "time" = slrc.timing
  ),
  paste0(output_dir, "/slr_cut_hdr_metrics", file.end))
  
  
  
  
  
  
  
  
  
  
  
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   spectral.clustering.algorithm == "kmeans"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrk = slr(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "kmeans")
  end.time = Sys.time()
  slrk.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrk.coefs = getCoefsBM(
    coefs = coefficients(slrk$model), sbp = slrk$sbp)
  
  # compute metrics on the selected model #
  slrk.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrk$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrk$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrk.coefs$a0, thetahat = slrk.coefs$bm.coefs,
    betahat = slrk.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrk.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrk.coefs$bm.coefs != 0),
    "time" = slrk.timing
  ),
  paste0(output_dir, "/slr_maxGamma_kmeans_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   spectral.clustering.algorithm == "kmeans"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrk = slr(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "kmeans")
  end.time = Sys.time()
  slrk.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrk.coefs = getCoefsBM(
    coefs = coefficients(slrk$model), sbp = slrk$sbp)
  
  # compute metrics on the selected model #
  slrk.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrk$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrk$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrk.coefs$a0, thetahat = slrk.coefs$bm.coefs,
    betahat = slrk.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrk.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrk.coefs$bm.coefs != 0),
    "time" = slrk.timing
  ),
  paste0(output_dir, "/slr_maxGamma_kmeans_amini_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   spectral.clustering.algorithm == "kmeans"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- TRUE
  ##############################################################################
  start.time = Sys.time()
  slrk = slr(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, spectral.clustering.method = "kmeans")
  end.time = Sys.time()
  slrk.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrk.coefs = getCoefsBM(
    coefs = coefficients(slrk$model), sbp = slrk$sbp)
  
  # compute metrics on the selected model #
  slrk.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrk$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrk$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrk.coefs$a0, thetahat = slrk.coefs$bm.coefs,
    betahat = slrk.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrk.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrk.coefs$bm.coefs != 0),
    "time" = slrk.timing
  ),
  paste0(output_dir, "/slr_maxGamma_kmeans_hdr_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   spectral.clustering.algorithm == "cut"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrc = slr(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "cut")
  end.time = Sys.time()
  slrc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrc.coefs = getCoefsBM(
    coefs = coefficients(slrc$model), sbp = slrc$sbp)
  
  # compute metrics on the selected model #
  slrc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrc.coefs$a0, thetahat = slrc.coefs$bm.coefs,
    betahat = slrc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrc.coefs$bm.coefs != 0),
    "time" = slrc.timing
  ),
  paste0(output_dir, "/slr_maxGamma_cut_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   spectral.clustering.algorithm == "cut"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrc = slr(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, spectral.clustering.method = "cut")
  end.time = Sys.time()
  slrc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrc.coefs = getCoefsBM(
    coefs = coefficients(slrc$model), sbp = slrc$sbp)
  
  # compute metrics on the selected model #
  slrc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrc.coefs$a0, thetahat = slrc.coefs$bm.coefs,
    betahat = slrc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrc.coefs$bm.coefs != 0),
    "time" = slrc.timing
  ),
  paste0(output_dir, "/slr_maxGamma_cut_amini_metrics", file.end))
  
  ##############################################################################
  # slr method
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   spectral.clustering.algorithm == "cut"
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slrc = slr(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = TRUE, spectral.clustering.method = "cut")
  end.time = Sys.time()
  slrc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrc.coefs = getCoefsBM(
    coefs = coefficients(slrc$model), sbp = slrc$sbp)
  
  # compute metrics on the selected model #
  slrc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slrc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slrc.coefs$a0, thetahat = slrc.coefs$bm.coefs,
    betahat = slrc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrc.coefs$bm.coefs != 0),
    "time" = slrc.timing
  ),
  paste0(output_dir, "/slr_maxGamma_cut_hdr_metrics", file.end))
  
  
  
  
  
  
  
  
  
  
  
  
  ##############################################################################
  # slr method using k-means spectral clustering with K = 3
  #   similarity.matrix -- TRUE
  #   maxGamma -- FALSE
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr1sc = slr1sc(
    similarity.matrix = TRUE, maxGamma = FALSE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE)
  end.time = Sys.time()
  slr1sc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1sc.coefs = getCoefsBM(
    coefs = coefficients(slr1sc$model), sbp = slr1sc$sbp)
  
  # compute metrics on the selected model #
  slr1sc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1sc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1sc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1sc.coefs$a0, thetahat = slr1sc.coefs$bm.coefs,
    betahat = slr1sc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1sc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1sc.coefs$bm.coefs != 0),
    "time" = slr1sc.timing
  ),
  paste0(output_dir, "/slr_1sc_metrics", file.end))
  
  ##############################################################################
  # slr method using k-means spectral clustering with K = 3
  #   similarity.matrix -- TRUE
  #   maxGamma -- TRUE
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr1sc = slr1sc(
    similarity.matrix = TRUE, maxGamma = TRUE,
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE)
  end.time = Sys.time()
  slr1sc.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1sc.coefs = getCoefsBM(
    coefs = coefficients(slr1sc$model), sbp = slr1sc$sbp)
  
  # compute metrics on the selected model #
  slr1sc.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1sc$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1sc$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1sc.coefs$a0, thetahat = slr1sc.coefs$bm.coefs,
    betahat = slr1sc.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1sc.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1sc.coefs$bm.coefs != 0),
    "time" = slr1sc.timing
  ),
  paste0(output_dir, "/slr_1sc_maxGamma_metrics", file.end))
  
}

