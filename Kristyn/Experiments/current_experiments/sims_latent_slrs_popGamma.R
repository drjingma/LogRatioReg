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
  source("Kristyn/Functions/popGamma.R")
  
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
  # population Gamma
  popGamma1 = popGamma(
    alpha1 = c(a1, 0), beta1 = b1, var_epsilon = sigma_eps1^2, 
    var_epsilon2 = sigma_eps2^2, U = U.all[1:n, ])
  rownames(popGamma1) <- colnames(popGamma1) <- colnames(X)
  
  ##############################################################################
  # slr method using population Gamma (subtractFrom1 = FALSE)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0 = slr_popGamma(
    popGamma = popGamma1, subtractFrom1 = FALSE, 
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0.coefs = getCoefsBM(
    coefs = coefficients(slr0$model), sbp = slr0$sbp)
  
  # compute metrics on the selected model #
  slr0.metrics = c(
    
  )
  
  saveRDS(c(
    slr0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0.coefs$bm.coefs != 0),
    "time" = slr0.timing
  ),
  paste0(output_dir, "/slr_popGamma_metrics", file.end))
  
  ##############################################################################
  # slr method using population Gamma (subtractFrom1 = TRUE)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- FALSE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr1 = slr_popGamma(
    popGamma = popGamma1, subtractFrom1 = TRUE, 
    x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1.coefs = getCoefsBM(
    coefs = coefficients(slr1$model), sbp = slr1$sbp)
  
  # compute metrics on the selected model #
  slr1.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1.coefs$a0, thetahat = slr1.coefs$bm.coefs,
    betahat = slr1.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1.coefs$bm.coefs != 0),
    "time" = slr1.timing
  ),
  paste0(output_dir, "/slr_popGamma_subtractfrom1_metrics", file.end))
  
  ##############################################################################
  # slr method using population Gamma (subtractFrom1 = FALSE)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0am = slr_popGamma(
    popGamma = popGamma1, subtractFrom1 = FALSE, 
    x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr0am.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0am.coefs = getCoefsBM(
    coefs = coefficients(slr0am$model), sbp = slr0am$sbp)
  
  # compute metrics on the selected model #
  slr0am.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0am$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0am$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0am.coefs$a0, thetahat = slr0am.coefs$bm.coefs,
    betahat = slr0am.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0am.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0am.coefs$bm.coefs != 0),
    "time" = slr0am.timing
  ),
  paste0(output_dir, "/slr_popGamma_amini_metrics", file.end))
  
  ##############################################################################
  # slr method using population Gamma (subtractFrom1 = TRUE)
  #   rank 1 approximation -- FALSE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr1am = slr_popGamma(
    popGamma = popGamma1, subtractFrom1 = TRUE, 
    x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr1am.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1am.coefs = getCoefsBM(
    coefs = coefficients(slr1am$model), sbp = slr1am$sbp)
  
  # compute metrics on the selected model #
  slr1am.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1am$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1am$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1am.coefs$a0, thetahat = slr1am.coefs$bm.coefs,
    betahat = slr1am.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1am.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1am.coefs$bm.coefs != 0),
    "time" = slr1am.timing
  ),
  paste0(output_dir, "/slr_popGamma_subtractfrom1_amini_metrics", file.end))
  
  ##############################################################################
  # slr method using population Gamma (subtractFrom1 = FALSE)
  #   rank 1 approximation -- TRUE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr0amap = slr_popGamma(
    popGamma = popGamma1, subtractFrom1 = FALSE, 
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
  paste0(output_dir, "/slr_popGamma_amini_approx_metrics", file.end))
  
  ##############################################################################
  # slr method using population Gamma (subtractFrom1 = TRUE)
  #   rank 1 approximation -- TRUE
  #   amini regularization -- TRUE
  #   high degree regularization -- FALSE
  #   include leading eigenvector -- FALSE
  ##############################################################################
  start.time = Sys.time()
  slr1amap = slr_popGamma(
    popGamma = popGamma1, subtractFrom1 = TRUE, 
    x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
    highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)
  end.time = Sys.time()
  slr1amap.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr1amap.coefs = getCoefsBM(
    coefs = coefficients(slr1amap$model), sbp = slr1amap$sbp)
  
  # compute metrics on the selected model #
  slr1amap.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr1amap$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr1amap$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr1amap.coefs$a0, thetahat = slr1amap.coefs$bm.coefs,
    betahat = slr1amap.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr1amap.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr1amap.coefs$bm.coefs != 0),
    "time" = slr1amap.timing
  ),
  paste0(output_dir, "/slr_popGamma_subtractfrom1_amini_approx_metrics", file.end))
}

