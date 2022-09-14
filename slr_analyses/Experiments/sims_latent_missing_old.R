# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 8/16/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Experiments/outputs/metrics_unlabeled"

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
numSims = 25 #100

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
  source("slr_analyses/Functions/slrs.R")
  source("slr_analyses/Functions/codalasso.R")
  source("slr_analyses/Functions/util.R")
  
  # Tuning parameters ##########################################################
  
  # Settings to toggle with
  settings.name = "ContinuousResponseUnlabeled"
  n = 30
  p = 30
  K = 10
  nlam = 100
  intercept = TRUE
  scaling = TRUE
  tol = 1e-4
  sigma_y = 0.001 # sigma (for y)
  sigma_x = 0.15 # sigma_j (for x)
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  # SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  # (b1 = 0.5, theta.value = 0.5, a0 = 0, prop.missing = 0.75, ulimit = 0.5)
  # (b1 = 1, theta.value = 0.3, a0 = 0, prop.missing = 0.70, 0.75, ulimit = 0.5)
  b1 = 1 # 0.5
  theta.value = 0.5 # weight on a1 -- 0.5
  a0 = 0 # 0
  n.unlabeled = n * 100
  ulimit = 0.5
  
  file.end = paste0(
    "_", settings.name,
    "_n", n.unlabeled, "unlabeled",
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_dim", n, "x", p, 
    "_ulimit", ulimit,
    "_noisey", sigma_y, 
    "_noisex", sigma_x, 
    "_b0", b0, 
    "_b1", b1, 
    "_a0", a0, 
    "_theta", theta.value,
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # get latent variable
  U.all = matrix(runif(min = -ulimit, max = ulimit, 2 * n + n.unlabeled), ncol = 1)
  # simulate y from latent variable
  y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n + n.unlabeled) * sigma_y)
  # simulate X: 
  epsj.all = matrix(rnorm((2 * n + n.unlabeled) * (p - 1)), nrow = (2 * n + n.unlabeled)) * sigma_x
  a1 = theta.value * ilrtrans.true$ilr.trans[-p] 
  #   alpha1j = {
  #     c1=theta*ilr.const/k+   if j \in I+
  #     -c2=-theta*ilr.const/k-  if j \in I-
  #     0                       o/w
  #   }
  alrXj.all = a0 + U.all %*% t(a1) + epsj.all #log(Xj/Xp) =alpha0j+alpha1j*U+epsj
  X.all <- alrinv(alrXj.all)
  colnames(X.all) = paste0('s', 1:p)
  
  # subset out labeled and unlabeled sets
  X.labeled = X.all[1:(2 * n), ]
  Y.labeled = y.all[1:(2 * n)]
  X2 = X.all[-(1:(2 * n)), ]
  # subset out training and test sets from the labeled data
  X = X.labeled[1:n, ]
  X.test = X.labeled[-(1:n), ]
  Y <- Y.labeled[1:n]
  Y.test <- Y.labeled[-(1:n)]
  
  # about linear log-contrast models' coefficients
  llc.coefs.non0 = as.vector(SBP.true != 0)
  # solve for beta
  c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
  llc.coefs.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
    as.vector(ilrtrans.true$ilr.trans)
  
  saveRDS(list(
    X = X, Y = Y, X.test = X.test, Y.test = Y.test, 
    SBP.true = SBP.true, beta.true = llc.coefs.true, 
    llc.coefs.non0 = llc.coefs.non0
  ),
  paste0(output_dir, "/data", file.end))
  
  ##############################################################################
  # about the chosen settings
  
  # # Aitchison variation when j != k 
  # # when j != k, aitchison var is Sjk = (c1 + c2)^2 Var[U] + 2 * sigma_eps2
  # varU = (2 * ulimit)^2 / 12
  # c1plusc2^2 * varU # term 1
  # 2 * sigma_x^2 # term 2 (want this term to dominate)
  # 
  # # Correlation bt clr(Xj) & y 
  # covclrXy = a1 * b1 * varU # covariance, in numerator
  # varclrX = a1^2 * varU + (1 - (1 / (p))) * sigma_x^2 # variance of clrX
  # vary = b1^2 * varU + sigma_y^2 # variance of y
  # # population correlations?
  # covclrXy / (sqrt(varclrX) * sqrt(vary))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  # fit models
  
  ##############################################################################
  # compositional lasso
  # -- fits a linear log contrast model
  ##############################################################################
  start.time = Sys.time()
  classo = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1),
    nlam = nlam, nfolds = K, tol = tol, intercept = intercept,
    scaling = scaling)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # cl.lam.idx = which.min(classo$cvm)
  oneSErule = min(classo$cvm) + classo$cvsd[which.min(classo$cvm)] * 1
  cl.lam.idx = which(classo$cvm <= oneSErule)[1]
  cl.a0 = classo$int[cl.lam.idx]
  cl.betahat = classo$bet[, cl.lam.idx]
  
  # compute metrics on the selected model #
  # prediction error
  cl.Yhat.test = cl.a0 + log(X.test) %*% cl.betahat
  cl.MSE.test = as.vector(crossprod(Y.test - cl.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  cl.metrics = getMetricsLLC(
    est.llc.coefs = cl.betahat,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = cl.MSE.test,
    cl.metrics,
    "logratios" = 0,
    "time" = cl.timing,
    "randindex" = NA,
    "adjrandindex" = NA
  ),
  paste0(output_dir, "/classo_metrics", file.end))
  
  ##############################################################################
  # slr
  #   screen.method = "wald"
  #   cluster.method = "spectral"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slrspeccv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  slrspec = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = slrspeccv$threshold[slrspeccv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  slrspec.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrspec.fullSBP) = colnames(X)
  slrspec.fullSBP[match(
    names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp
  slrspec.coefs = getCoefsBM(
    coefs = coefficients(slrspec$fit), sbp = slrspec.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  slrspec.Yhat.test = slrspec.coefs$a0 + 
    slr.fromContrast(X.test, slrspec.fullSBP) * slrspec.coefs$bm.coefs
  slrspec.MSE.test = as.vector(crossprod(Y.test - slrspec.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  slrspec.metrics = getMetricsBM(
    est.llc.coefs = slrspec.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = slrspec.MSE.test,
    slrspec.metrics,
    "logratios" = sum(slrspec.coefs$bm.coefs != 0),
    "time" = slrspec.timing, 
    "randindex" = randidx(
      SBP.true, slrspec.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, slrspec.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/slr_spectral_metrics", file.end))
  
  ##############################################################################
  # slr
  #   screen.method = "wald"
  #   cluster.method = "hierarchical"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slrhiercv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  slrhier = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = slrhiercv$threshold[slrhiercv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  slrhier.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrhier.fullSBP) = colnames(X)
  slrhier.fullSBP[match(
    names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
  slrhier.coefs = getCoefsBM(
    coefs = coefficients(slrhier$fit), sbp = slrhier.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  slrhier.Yhat.test = slrhier.coefs$a0 + 
    slr.fromContrast(X.test, slrhier.fullSBP) * slrhier.coefs$bm.coefs
  slrhier.MSE.test = as.vector(crossprod(Y.test - slrhier.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  slrhier.metrics = getMetricsBM(
    est.llc.coefs = slrhier.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = slrhier.MSE.test,
    slrhier.metrics,
    "logratios" = sum(slrhier.coefs$bm.coefs != 0),
    "time" = slrhier.timing, 
    "randindex" = randidx(
      SBP.true, slrhier.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, slrhier.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/slr_hierarchical_metrics", file.end))
  
  ##############################################################################
  # semislr - don't use unlabeled data in CV
  #   screen.method = "wald"
  #   cluster.method = "spectral"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  semislrspec0cv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  semislrspec0 = slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, 
    screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = semislrspec0cv$threshold[semislrspec0cv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  semislrspec0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  semislrspec0.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(semislrspec0.fullSBP) = colnames(X)
  semislrspec0.fullSBP[match(
    names(semislrspec0$sbp), rownames(semislrspec0.fullSBP))] = semislrspec0$sbp
  semislrspec0.coefs = getCoefsBM(
    coefs = coefficients(semislrspec0$fit), sbp = semislrspec0.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  semislrspec0.Yhat.test = semislrspec0.coefs$a0 + 
    slr.fromContrast(X.test, semislrspec0.fullSBP) * semislrspec0.coefs$bm.coefs
  semislrspec0.MSE.test = as.vector(crossprod(Y.test - semislrspec0.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  semislrspec0.metrics = getMetricsBM(
    est.llc.coefs = semislrspec0.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = semislrspec0.MSE.test,
    semislrspec0.metrics,
    "logratios" = sum(semislrspec0.coefs$bm.coefs != 0),
    "time" = semislrspec0.timing, 
    "randindex" = randidx(
      SBP.true, semislrspec0.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, semislrspec0.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/semislr_spectral_noCVUnlabeled_metrics", file.end))
  
  ##############################################################################
  # semislr - don't use unlabeled data in CV
  #   screen.method = "wald"
  #   cluster.method = "hierarchical"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  semislrhier0cv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  semislrhier0 = slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, 
    screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = semislrhier0cv$threshold[semislrhier0cv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  semislrhier0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  semislrhier0.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(semislrhier0.fullSBP) = colnames(X)
  semislrhier0.fullSBP[match(
    names(semislrhier0$sbp), rownames(semislrhier0.fullSBP))] = semislrhier0$sbp
  semislrhier0.coefs = getCoefsBM(
    coefs = coefficients(semislrhier0$fit), sbp = semislrhier0.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  semislrhier0.Yhat.test = semislrhier0.coefs$a0 + 
    slr.fromContrast(X.test, semislrhier0.fullSBP) * semislrhier0.coefs$bm.coefs
  semislrhier0.MSE.test = as.vector(crossprod(Y.test - semislrhier0.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  semislrhier0.metrics = getMetricsBM(
    est.llc.coefs = semislrhier0.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = semislrhier0.MSE.test,
    semislrhier0.metrics,
    "logratios" = sum(semislrhier0.coefs$bm.coefs != 0),
    "time" = semislrhier0.timing, 
    "randindex" = randidx(
      SBP.true, semislrhier0.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, semislrhier0.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/semislr_hierarchical_noCVUnlabeled_metrics", file.end))
  
  
  
  ##############################################################################
  # semislr - use unlabeled data in CV (but don't fold it)
  #   screen.method = "wald"
  #   cluster.method = "spectral"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  semislrspec1cv = cv.slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, fold.unlabeled = FALSE,
    screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  semislrspec1 = slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, 
    screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = semislrspec1cv$threshold[semislrspec1cv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  semislrspec1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  semislrspec1.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(semislrspec1.fullSBP) = colnames(X)
  semislrspec1.fullSBP[match(
    names(semislrspec1$sbp), rownames(semislrspec1.fullSBP))] = semislrspec1$sbp
  semislrspec1.coefs = getCoefsBM(
    coefs = coefficients(semislrspec1$fit), sbp = semislrspec1.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  semislrspec1.Yhat.test = semislrspec1.coefs$a0 + 
    slr.fromContrast(X.test, semislrspec1.fullSBP) * semislrspec1.coefs$bm.coefs
  semislrspec1.MSE.test = as.vector(crossprod(Y.test - semislrspec1.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  semislrspec1.metrics = getMetricsBM(
    est.llc.coefs = semislrspec1.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = semislrspec1.MSE.test,
    semislrspec1.metrics,
    "logratios" = sum(semislrspec1.coefs$bm.coefs != 0),
    "time" = semislrspec1.timing, 
    "randindex" = randidx(
      SBP.true, semislrspec1.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, semislrspec1.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/semislr_spectral_CVUnlabeledNotFolded_metrics", file.end))
  
  ##############################################################################
  # semislr - use unlabeled data in CV (but don't fold it)
  #   screen.method = "wald"
  #   cluster.method = "hierarchical"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  semislrhier1cv = cv.slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, fold.unlabeled = FALSE,
    screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  semislrhier1 = slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, 
    screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = semislrhier1cv$threshold[semislrhier1cv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  semislrhier1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  semislrhier1.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(semislrhier1.fullSBP) = colnames(X)
  semislrhier1.fullSBP[match(
    names(semislrhier1$sbp), rownames(semislrhier1.fullSBP))] = semislrhier1$sbp
  semislrhier1.coefs = getCoefsBM(
    coefs = coefficients(semislrhier1$fit), sbp = semislrhier1.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  semislrhier1.Yhat.test = semislrhier1.coefs$a0 + 
    slr.fromContrast(X.test, semislrhier1.fullSBP) * semislrhier1.coefs$bm.coefs
  semislrhier1.MSE.test = as.vector(crossprod(Y.test - semislrhier1.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  semislrhier1.metrics = getMetricsBM(
    est.llc.coefs = semislrhier1.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = semislrhier1.MSE.test,
    semislrhier1.metrics,
    "logratios" = sum(semislrhier1.coefs$bm.coefs != 0),
    "time" = semislrhier1.timing, 
    "randindex" = randidx(
      SBP.true, semislrhier1.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, semislrhier1.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/semislr_hierarchical_CVUnlabeledNotFolded_metrics", file.end))
  
  
  
  ##############################################################################
  # semislr - use unlabeled data in CV (and fold it)
  #   screen.method = "wald"
  #   cluster.method = "spectral"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  semislrspec2cv = cv.slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, fold.unlabeled = TRUE,
    screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  semislrspec2 = slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, 
    screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = semislrspec2cv$threshold[semislrspec2cv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  semislrspec2.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  semislrspec2.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(semislrspec2.fullSBP) = colnames(X)
  semislrspec2.fullSBP[match(
    names(semislrspec2$sbp), rownames(semislrspec2.fullSBP))] = semislrspec2$sbp
  semislrspec2.coefs = getCoefsBM(
    coefs = coefficients(semislrspec2$fit), sbp = semislrspec2.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  semislrspec2.Yhat.test = semislrspec2.coefs$a0 + 
    slr.fromContrast(X.test, semislrspec2.fullSBP) * semislrspec2.coefs$bm.coefs
  semislrspec2.MSE.test = as.vector(crossprod(Y.test - semislrspec2.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  semislrspec2.metrics = getMetricsBM(
    est.llc.coefs = semislrspec2.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = semislrspec2.MSE.test,
    semislrspec2.metrics,
    "logratios" = sum(semislrspec2.coefs$bm.coefs != 0),
    "time" = semislrspec2.timing, 
    "randindex" = randidx(
      SBP.true, semislrspec2.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, semislrspec2.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/semislr_spectral_CVUnlabeledFolded_metrics", file.end))
  
  ##############################################################################
  # semislr - use unlabeled data in CV (and fold it)
  #   screen.method = "wald"
  #   cluster.method = "hierarchical"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  semislrhier2cv = cv.slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, fold.unlabeled = TRUE,
    screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    scale = scaling, trace.it = FALSE)
  semislrhier2 = slr(
    x = X, y = Y, 
    x.unlabeled = X2, use.unlabeled = TRUE, 
    screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = semislrhier2cv$threshold[semislrhier2cv$index["1se",]], 
    positive.slope = TRUE)
  end.time = Sys.time()
  semislrhier2.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  semislrhier2.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(semislrhier2.fullSBP) = colnames(X)
  semislrhier2.fullSBP[match(
    names(semislrhier2$sbp), rownames(semislrhier2.fullSBP))] = semislrhier2$sbp
  semislrhier2.coefs = getCoefsBM(
    coefs = coefficients(semislrhier2$fit), sbp = semislrhier2.fullSBP)
  
  # compute metrics on the selected model #
  # prediction error
  semislrhier2.Yhat.test = semislrhier2.coefs$a0 + 
    slr.fromContrast(X.test, semislrhier2.fullSBP) * semislrhier2.coefs$bm.coefs
  semislrhier2.MSE.test = as.vector(crossprod(Y.test - semislrhier2.Yhat.test) / n)
  # estimation accuracy, selection accuracy #
  semislrhier2.metrics = getMetricsBM(
    est.llc.coefs = semislrhier2.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = semislrhier2.MSE.test,
    semislrhier2.metrics,
    "logratios" = sum(semislrhier2.coefs$bm.coefs != 0),
    "time" = semislrhier2.timing,
    "randindex" = randidx(
      SBP.true, semislrhier2.fullSBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, semislrhier2.fullSBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/semislr_hierarchical_CVUnlabeledFolded_metrics", file.end))
  
  ##############################################################################
  # CoDaCoRe
  # -- fits a balance regression model with possibly multiple balances
  ##############################################################################
  library(codacore)
  
  if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
    reticulate::use_condaenv("anaconda3")
  }
  
  start.time = Sys.time()
  codacore0 = codacore(
    x = X, y = Y, logRatioType = "ILR",
    objective = "regression", cvParams = list(numFolds = K)) 
  end.time = Sys.time()
  codacore0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  if(length(codacore0$ensemble) > 0){ # at least 1 log-ratio found
    codacore0_SBP = matrix(0, nrow = p, ncol = length(codacore0$ensemble))
    codacore0_coeffs = rep(NA, length(codacore0$ensemble))
    for(col.idx in 1:ncol(codacore0_SBP)){
      codacore0_SBP[
        codacore0$ensemble[[col.idx]]$hard$numerator, col.idx] = 1
      codacore0_SBP[
        codacore0$ensemble[[col.idx]]$hard$denominator, col.idx] = -1
      codacore0_coeffs[col.idx] = codacore0$ensemble[[col.idx]]$slope
    }
    
    names(codacore0_coeffs) = paste(
      "balance", 1:length(codacore0_coeffs), sep = "")
    rownames(codacore0_SBP) = colnames(X)
    codacore0.coefs2 = getCoefsBM(
      coefs = codacore0_coeffs * codacore0$yScale, 
      sbp = codacore0_SBP)
    codacore0.betahat = codacore0.coefs2$llc.coefs
    
    # compute metrics on the selected model #
    # prediction error
    codacore0.Yhat.test = predict(codacore0, X.test)
    
  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0_coeffs = c()
    codacore0model = stats::glm(Y ~ 1, family = "gaussian")
    codacore0.betahat = rep(0, p)
    
    # compute metrics on the selected model #
    # prediction error
    codacore0.Yhat.test = predict(codacore0model, X.test)
  }
  codacore0.MSE.test = crossprod(Y.test - codacore0.Yhat.test) / n
  
  # beta estimation accuracy, selection accuracy #
  codacore0.metrics = getMetricsBM(
    est.llc.coefs = codacore0.betahat,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true, metrics = c("estimation", "selection"))
  
  saveRDS(c(
    "mse" = codacore0.MSE.test,
    codacore0.metrics,
    "logratios" = length(codacore0_coeffs), 
    "time" = codacore0.timing,
    "randindex" = randidx(
      SBP.true, codacore0_SBP[, 1, drop = FALSE], adjusted = FALSE),
    "adjrandindex" = randidx(
      SBP.true, codacore0_SBP[, 1, drop = FALSE], adjusted = TRUE)
  ),
  paste0(output_dir, "/codacore_metrics", file.end))
  
  
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}
