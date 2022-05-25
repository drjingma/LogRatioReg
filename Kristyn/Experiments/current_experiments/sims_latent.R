# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 5/24/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics"

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
  
  source("RCode/func_libs_1.R") # for classo to work
  
  source("Kristyn/Functions/slr.R")
  source("Kristyn/Functions/util.R")
  source("Kristyn/Functions/slrscreen.R")
  
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
  sigma_eps1 = 0.1 # sigma (for y)
  sigma_eps2 = 0.1 # sigma_j (for x)
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  # SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 0.5 # 1, 0.5, 0.25
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
  
  saveRDS(list(
    X = X, Y = Y, X.test = X.test, Y.test = Y.test, 
    SBP.true = SBP.true, beta.true = beta.true, 
    is0.beta = is0.beta, non0.beta = non0.beta
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
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta, 
    true.beta = beta.true)
  
  saveRDS(c(
    cl.metrics,
    "betasparsity" = bspars,
    "logratios" = 0,
    "time" = cl.timing, 
    "adhoc" = NA
  ),
  paste0(output_dir, "/classo_metrics", file.end))
  
  ##############################################################################
  # slr method using k-means spectral clustering with K = 3
  #   alpha = 0.05
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slr0.05 = slr(
    x = X, y = Y, alpha = 0.05)
  end.time = Sys.time()
  slr0.05.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0.05.coefs = getCoefsBM(
    coefs = coefficients(slr0.05$model), sbp = slr0.05$sbp)
  
  # compute metrics on the selected model #
  slr0.05.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0.05$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0.05$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0.05.coefs$a0, thetahat = slr0.05.coefs$bm.coefs,
    betahat = slr0.05.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0.05.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0.05.coefs$bm.coefs != 0),
    "time" = slr0.05.timing, 
    "adhoc" = slr0.05$adhoc.invoked
  ),
  paste0(output_dir, "/slr_alpha0.05_metrics", file.end))
  
  ##############################################################################
  # slr method using k-means spectral clustering with K = 3
  #   alpha = 0.01
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slr0.01 = slr(x = X, y = Y, alpha = 0.01)
  end.time = Sys.time()
  slr0.01.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0.01.coefs = getCoefsBM(
    coefs = coefficients(slr0.01$model), sbp = slr0.01$sbp)
  
  # compute metrics on the selected model #
  slr0.01.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slr0.01$sbp),
    ilrX.test = getIlrX(X.test, sbp = slr0.01$sbp),
    n.train = n, n.test = n,
    thetahat0 = slr0.01.coefs$a0, thetahat = slr0.01.coefs$bm.coefs,
    betahat = slr0.01.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slr0.01.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0.01.coefs$bm.coefs != 0),
    "time" = slr0.01.timing, 
    "adhoc" = slr0.01$adhoc.invoked
  ),
  paste0(output_dir, "/slr_alpha0.01_metrics", file.end))
  
  ##############################################################################
  # slr method with screening step
  #   method = "wald"
  #   response.type = "continuous"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slrscreen0cv = cv.slr.screen(
    x = X, y = Y, method = "wald", 
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrscreen0 = slr.screen(
    x = X, y = Y, method = "wald", 
    response.type = "continuous", s0.perc = 0, zeta = 0, 
    threshold = slrscreen0cv$threshold[slrscreen0cv$index["1se",]])
  end.time = Sys.time()
  slrscreen0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrscreen0.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrscreen0.fullSBP) = colnames(X)
  slrscreen0.fullSBP[match(
    names(slrscreen0$sbp), rownames(slrscreen0.fullSBP))] = slrscreen0$sbp
  
  slrscreen0.coefs = getCoefsBM(
    coefs = coefficients(slrscreen0$fit), sbp = slrscreen0.fullSBP)
  
  # compute metrics on the selected model #
  slrscreen0.metrics = getMetricsBM(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = slrscreen0.fullSBP),
    ilrX.test = getIlrX(X.test, sbp = slrscreen0.fullSBP),
    n.train = n, n.test = n,
    thetahat0 = slrscreen0.coefs$a0, thetahat = slrscreen0.coefs$bm.coefs,
    betahat = slrscreen0.coefs$llc.coefs,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)
  
  saveRDS(c(
    slrscreen0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slrscreen0.coefs$bm.coefs != 0),
    "time" = slrscreen0.timing, 
    "adhoc" = slrscreen0$adhoc.invoked
  ),
  paste0(output_dir, "/slrscreen_metrics", file.end))
  
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
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slbl.metrics = c(PEtr = slbl.PE.train, PEte = slbl.PE.test, slbl.metrics)
  
  saveRDS(c(
    slbl.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slbl.coefs$bm.coefs != 0), 
    "time" = slbl.timing, 
    "adhoc" = NA
  ),
  paste0(output_dir, "/selbal_metrics", file.end))
  
  ##############################################################################
  # CoDaCoRe
  # -- fits a balance regression model with possibly multiple balances
  ##############################################################################
  library(codacore)
  
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
    
    codacore0.betahat = getBetaFromCodacore(
      SBP_codacore = codacore0_SBP, coeffs_codacore = codacore0_coeffs, p = p)
    
    # compute metrics on the selected model #
    # prediction errors
    # get prediction error on training set
    codacore0.Yhat.train = predict(codacore0, X)
    codacore0.PE.train = crossprod(Y - codacore0.Yhat.train) / n
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0, X.test)
    codacore0.PE.test = crossprod(Y - codacore0.Yhat.test) / n
    
  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0_coeffs = c()
    codacore0model = stats::glm(Y ~ 1, family = "gaussian")
    codacore0.betahat = rep(0, p)
    
    # compute metrics on the selected model #
    # prediction errors
    # get prediction error on training set
    codacore0.Yhat.train = predict(codacore0model, data.frame(X))
    codacore0.PE.train = crossprod(Y - codacore0.Yhat.train) / n
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0model, data.frame(X.test))
    codacore0.PE.test = crossprod(Y - codacore0.Yhat.test) / n
    
  }
  # beta estimation accuracy, selection accuracy #
  codacore0.metrics = getMetricsBM(
    betahat = codacore0.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  codacore0.metrics = c(
    PEtr = codacore0.PE.train, PEte = codacore0.PE.test, codacore0.metrics)
  
  saveRDS(c(
    codacore0.metrics,
    "betasparsity" = bspars,
    "logratios" = length(codacore0_coeffs), 
    "time" = codacore0.timing, 
    "adhoc" = NA
  ),
  paste0(output_dir, "/codacore_metrics", file.end))
  
  # ##############################################################################
  # # Log-Ratio Lasso
  # # -- regresses on pairwise log-ratios
  # ##############################################################################
  # library(logratiolasso)
  # source("Kristyn/Functions/logratiolasso.R")
  # Wc = scale(log(X), center = TRUE, scale = FALSE)
  # Yc = Y - mean(Y)
  # 
  # start.time = Sys.time()
  # lrl <- cv_two_stage(z = Wc, y = Yc, n_folds = K)
  # end.time = Sys.time()
  # lrl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # compute metrics on the selected model #
  # # prediction errors
  # # get prediction error on training set
  # lrl.Yhat.train = Wc %*% lrl$beta_min
  # lrl.PE.train = crossprod(Y - lrl.Yhat.train) / n
  # # get prediction error on test set
  # Wc.test = scale(log(X.test), center = TRUE, scale = FALSE)
  # Yc.test = Y.test - mean(Y.test)
  # lrl.Yhat.test = Wc.test %*% lrl$beta_min
  # lrl.PE.test = crossprod(Y - lrl.Yhat.test) / n
  # 
  # # beta estimation accuracy, selection accuracy #
  # lrl.metrics = getMetricsBM(
  #   betahat = lrl$beta_min, # don't back-scale bc only centered X (didn't scale)
  #   true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
  #   true.beta = beta.true, metrics = c("betaestimation", "selection"))
  # lrl.metrics = c(
  #   PEtr = lrl.PE.train, PEte = lrl.PE.test, lrl.metrics)
  # 
  # saveRDS(c(
  #   lrl.metrics,
  #   "betasparsity" = bspars,
  #   "logratios" = NA, 
  #   "time" = lrl.timing, 
  #   "adhoc" = NA
  # ),
  # paste0(output_dir, "/lrlasso_metrics", file.end))
  
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}
