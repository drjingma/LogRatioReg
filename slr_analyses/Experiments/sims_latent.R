# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 8/16/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Experiments/outputs/metrics"

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
  source("slr_analyses/Functions/slrs.R")
  source("slr_analyses/Functions/codalasso.R")
  source("slr_analyses/Functions/util.R")
  
  # Tuning parameters###########################################################
  
  # Settings to toggle with
  sigma.settings = "latentVarModel"
  n = 100
  p = 30
  K = 10
  nlam = 100
  intercept = TRUE
  scaling = TRUE
  tol = 1e-4
  sigma_y = 0.1 # sigma (for y)
  sigma_x = 0.1 # sigma_j (for x)
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  # SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 0.5 # 0.5
  theta.value = 1 # weight on a1 -- 1
  a0 = 0 # 0
  ulimit = 0.5
  
  file.end = paste0(
    "_", sigma.settings,
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
  U.all = matrix(runif(min = -ulimit, max = ulimit, 2 * n), ncol = 1)
  # simulate y from latent variable
  y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_y)
  # simulate X: 
  epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_x
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
  bspars = sum(non0.beta)
  # solve for beta
  c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
  beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
    as.vector(ilrtrans.true$ilr.trans)
  
  saveRDS(list(
    X = X, Y = Y, X.test = X.test, Y.test = Y.test, 
    SBP.true = SBP.true, beta.true = beta.true, 
    non0.beta = non0.beta
  ),
  paste0(output_dir, "/data", file.end))
  
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
  cl.metrics = getMetricsLLC(
    y.train = Y, y.test = Y.test,
    logX.train = log(X),
    logX.test = log(X.test),
    n.train = n, n.test = n,
    betahat0 = cl.a0, betahat = cl.betahat,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    cl.metrics,
    "betasparsity" = bspars,
    "logratios" = 0,
    "time" = cl.timing
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
    threshold = slrspeccv$threshold[slrspeccv$index["1se",]])
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
  slrspec.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrspec.fullSBP),
    ilrX.test = getIlrX(X.test, sbp = slrspec.fullSBP),
    n.train = n, n.test = n,
    thetahat0 = slrspec.coefs$a0, thetahat = slrspec.coefs$bm.coefs,
    betahat = slrspec.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrspec.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrspec.coefs$bm.coefs != 0),
    "time" = slrspec.timing
  ),
  paste0(output_dir, "/slr_spectral_metrics", file.end))
  
  ##############################################################################
  # slr
  #   screen.method = "wald"
  #   cluster.method = "hieararchical"
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
    threshold = slrhiercv$threshold[slrhiercv$index["1se",]])
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
  slrhier.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrhier.fullSBP),
    ilrX.test = getIlrX(X.test, sbp = slrhier.fullSBP),
    n.train = n, n.test = n,
    thetahat0 = slrhier.coefs$a0, thetahat = slrhier.coefs$bm.coefs,
    betahat = slrhier.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrhier.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrhier.coefs$bm.coefs != 0),
    "time" = slrhier.timing
  ),
  paste0(output_dir, "/slr_hierarchical_metrics", file.end))
  
  ##############################################################################
  # selbal method (a balance regression method)
  # -- fits a balance regression model with one balance
  ##############################################################################
  library(selbal) # masks stats::cor()
  slbl.data = getSelbalData(X = X, y = Y, classification = FALSE)
  
  start.time = Sys.time()
  slbl = selbal.cv(x = slbl.data$X, y = slbl.data$y, n.fold = K)
  end.time = Sys.time()
  slbl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slbl.coefs = getCoefsSelbal(
    X = slbl.data$X, y = slbl.data$y, selbal.fit = slbl, classification = FALSE, 
    check = TRUE)
  
  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slbl.Yhat.train = predict.glm(
    slbl$glm, 
    newdata = data.frame(V1 = balance::balance.fromSBP(
      x = slbl.data$X, y = slbl.coefs$sbp)), 
    type = "response")
  slbl.PE.train = crossprod(Y - slbl.Yhat.train) / n
  # get prediction error on test set
  slbl.test.data = getSelbalData(X = X.test, y = Y.test, classification = FALSE)
  slbl.Yhat.test = predict.glm(
    slbl$glm, 
    newdata = data.frame(V1 = balance::balance.fromSBP(
      x = slbl.test.data$X, y = slbl.coefs$sbp)), 
    type = "response")
  slbl.PE.test = crossprod(Y.test - slbl.Yhat.test) / n
  # beta estimation accuracy, selection accuracy #
  slbl.metrics = getMetricsBM(
    thetahat = slbl.coefs$bm.coefs, betahat = slbl.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slbl.metrics = c(PEtr = slbl.PE.train, PEte = slbl.PE.test, slbl.metrics)
  
  saveRDS(c(
    slbl.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slbl.coefs$bm.coefs != 0), 
    "time" = slbl.timing
  ),
  paste0(output_dir, "/selbal_metrics", file.end))
  
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
    # prediction errors
    # get prediction error on training set
    codacore0.Yhat.train = predict(codacore0, X)
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0, X.test)
    
  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0_coeffs = c()
    codacore0model = stats::glm(Y ~ 1, family = "gaussian")
    codacore0.betahat = rep(0, p)
    
    # compute metrics on the selected model #
    # prediction errors
    # get prediction error on training set
    codacore0.Yhat.train = predict(codacore0model, X)
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0model, X.test)
  }
  codacore0.PE.train = crossprod(Y - codacore0.Yhat.train) / n
  codacore0.PE.test = crossprod(Y.test - codacore0.Yhat.test) / n
  
  # beta estimation accuracy, selection accuracy #
  codacore0.metrics = getMetricsBM(
    betahat = codacore0.betahat,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  codacore0.metrics = c(
    PEtr = codacore0.PE.train, PEte = codacore0.PE.test, codacore0.metrics)
  
  saveRDS(c(
    codacore0.metrics,
    "betasparsity" = bspars,
    "logratios" = length(codacore0_coeffs), 
    "time" = codacore0.timing
  ),
  paste0(output_dir, "/codacore_metrics", file.end))
  
  ##############################################################################
  # Log-Ratio Lasso
  # -- regresses on pairwise log-ratios
  ##############################################################################
  library(logratiolasso)
  source("slr_analyses/Functions/logratiolasso.R")
  Wc = scale(log(X), center = TRUE, scale = FALSE)
  Yc = Y - mean(Y)

  start.time = Sys.time()
  lrl <- cv_two_stage(z = Wc, y = Yc, n_folds = K)
  end.time = Sys.time()
  lrl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  lrl.Yhat.train = Wc %*% lrl$beta_min
  lrl.PE.train = crossprod(Y - lrl.Yhat.train) / n
  # get prediction error on test set
  Wc.test = scale(log(X.test), center = TRUE, scale = FALSE)
  Yc.test = Y.test - mean(Y.test)
  lrl.Yhat.test = Wc.test %*% lrl$beta_min
  lrl.PE.test = crossprod(Yc.test - lrl.Yhat.test) / n

  # beta estimation accuracy, selection accuracy #
  lrl.metrics = getMetricsBM(
    betahat = lrl$beta_min, # don't back-scale bc only centered X (didn't scale)
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  lrl.metrics = c(
    PEtr = lrl.PE.train, PEte = lrl.PE.test, lrl.metrics)

  saveRDS(c(
    lrl.metrics,
    "betasparsity" = bspars,
    "logratios" = NA,
    "time" = lrl.timing
  ),
  paste0(output_dir, "/lrlasso_metrics", file.end))
  
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}
