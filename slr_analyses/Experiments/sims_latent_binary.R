# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 6/16/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Experiments/outputs/metrics_binary"

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
  source("RCode/SLR.R")
  source("slr_analyses/Functions/codalasso.R")
  source("slr_analyses/Functions/util.R")
  
  # Tuning parameters###########################################################
  
  # Settings to toggle with
  sigma.settings = "latentVarModel_binary"
  n = 100
  p = 30
  K = 10
  nlam = 100
  neta = p
  intercept = TRUE
  scaling = TRUE
  tol = 1e-4
  # sigma_eps1 = 0.1
  sigma_eps2 = 0.1 # 0.1, 0.01
  # SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 6 # 6
  theta.value = 1 # weight on a1 -- 1
  a0 = 0 # 0
  
  file.end = paste0(
    "_", sigma.settings,
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_dim", n, "x", p, 
    # "_noisey", sigma_eps1, 
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
  y.all = rbinom(n = 2 * n, size = 1, p = as.vector(sigmoid(b0 + b1 * U.all)))
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
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  start.time = Sys.time()
  classo = codalasso(X, Y, numFolds = K)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  cl.betahat = classo$cll$betas[-1]
  
  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  classo.Yhat.train = predict(classo, X)
  classo.AUC.train = pROC::roc(
    Y, classo.Yhat.train, levels = c(0, 1), direction = "<")$auc
  # get prediction error on test set
  classo.Yhat.test = predict(classo, X.test)
  classo.AUC.test = pROC::roc(
    Y.test, classo.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  cl.metrics = getMetricsLLC(
    betahat = cl.betahat,
    true.sbp = SBP.true, non0.true.beta = non0.beta, 
    true.beta = beta.true, 
    metrics = c("betaestimation", "selection"), classification = TRUE)
  cl.metrics = c(
    AUCtr = classo.AUC.train, AUCte = classo.AUC.test, cl.metrics)

  saveRDS(c(
    cl.metrics,
    "betasparsity" = bspars,
     "logratios" = 0,
    "time" = cl.timing
  ),
  paste0(output_dir, "/classo_metrics", file.end))

  ##############################################################################
  # slr method with screening step
  #   method = "wald"
  #   response.type = "binary"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "mse"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slr0cv = cv.slr(
    x = X, y = Y, method = "wald", 
    response.type = "binary", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slr0 = slr(
    x = X, y = Y, method = "wald", 
    response.type = "binary", s0.perc = 0, zeta = 0, 
    threshold = slr0cv$threshold[slr0cv$index["1se",]])
  end.time = Sys.time()
  slr0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slr0.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slr0.fullSBP) = colnames(X)
  slr0.fullSBP[match(
    names(slr0$sbp), rownames(slr0.fullSBP))] = slr0$sbp
  
  slr0.coefs = getCoefsBM(
    coefs = coefficients(slr0$fit), sbp = slr0.fullSBP)
  
  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slr0.Yhat.train = predict(
    slr0$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = X, y = slr0.fullSBP)), 
    type = "response")
  slr0.AUC.train = pROC::roc(
    Y, slr0.Yhat.train, levels = c(0, 1), direction = "<")$auc
  # get prediction error on test set
  slr0.Yhat.test = predict(
    slr0$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = X.test, y = slr0.fullSBP)), 
    type = "response")
  slr0.AUC.test = pROC::roc(
    Y.test, slr0.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  slr0.metrics = getMetricsBM(
    thetahat0 = slr0.coefs$a0, thetahat = slr0.coefs$bm.coefs,
    betahat = slr0.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true,
    metrics = c("betaestimation", "selection"), classification = TRUE)
  slr0.metrics = c(
    AUCtr = slr0.AUC.train, AUCte = slr0.AUC.test, 
    slr0.metrics)
  
  saveRDS(c(
    slr0.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slr0.coefs$bm.coefs != 0),
    "time" = slr0.timing, 
    "adhoc" = slr0$adhoc.invoked
  ),
  paste0(output_dir, "/slr_metrics", file.end))

  ##############################################################################
  # selbal method (a balance regression method)
  # -- fits a balance regression model with one balance
  ##############################################################################
  slbl.data = getSelbalData(
    X = X, y = Y, classification = TRUE, levels = c(0, 1), labels = c(0, 1))
  
  start.time = Sys.time()
  slbl = selbal::selbal.cv(x = slbl.data$X, y = slbl.data$y, n.fold = K)
  end.time = Sys.time()
  slbl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slbl.coefs = getCoefsSelbal(
    X = slbl.data$X, y = slbl.data$y, selbal.fit = slbl, classification = TRUE, 
    check = TRUE)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slbl.Yhat.train = predict.glm(
    slbl$glm, 
    newdata = data.frame(V1 = balance::balance.fromSBP(
      x = slbl.data$X, y = slbl.coefs$sbp)), 
    type = "response")
  slbl.AUC.train = pROC::roc(
    slbl.data$y, slbl.Yhat.train, levels = c(0, 1), direction = "<")$auc
  # get prediction error on test set
  slbl.test.data = getSelbalData(
    X = X.test, y = Y.test, classification = TRUE, 
    levels = c(0, 1), labels = c(0, 1))
  slbl.Yhat.test = predict.glm(
    slbl$glm, 
    newdata = data.frame(V1 = balance::balance.fromSBP(
      x = slbl.test.data$X, y = slbl.coefs$sbp)), 
    type = "response")
  slbl.AUC.test = pROC::roc(
    slbl.test.data$y, slbl.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  slbl.metrics = getMetricsBM(
    thetahat = slbl.coefs$bm.coefs, betahat = slbl.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slbl.metrics = c(AUCtr = slbl.AUC.train, AUCte = slbl.AUC.test, slbl.metrics)
  
  saveRDS(c(
    slbl.metrics,
    "betasparsity" = bspars,
    "logratios" = sum(slbl.coefs$bm.coefs != 0), 
    "time" = slbl.timing
  ),
  paste0(output_dir, "/selbal_metrics", file.end))
  
  ##############################################################################
  # CoDaCoRe (a balance regression method)
  # -- fits a balance regression model with possibly multiple balances
  ##############################################################################
  start.time = Sys.time()
  codacore0 = codacore::codacore(
    x = X, y = Y, logRatioType = "ILR", 
    objective = "binary classification", cvParams = list(numFolds = K))
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
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0, X.test)
    
  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0_coeffs = c()
    codacore0model = stats::glm(Y ~ 1, family = "binomial")
    codacore0.betahat = rep(0, p)
    
    # compute metrics on the selected model #
    # prediction errors
    # get prediction error on training set
    codacore0.Yhat.train = predict(codacore0model, X)
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0model, X.test)
  }
  codacore0.AUC.train = pROC::roc(
    Y, codacore0.Yhat.train, levels = c(0, 1), direction = "<")$auc
  codacore0.AUC.test = pROC::roc(
    Y.test, codacore0.Yhat.test, levels = c(0, 1), direction = "<")$auc
  
  # beta estimation accuracy, selection accuracy #
  codacore0.metrics = getMetricsBM(
    betahat = codacore0.betahat,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  codacore0.metrics = c(
    AUCtr = codacore0.AUC.train, AUCte = codacore0.AUC.test, codacore0.metrics)
  
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
  
  start.time = Sys.time()
  lrl <- cv_two_stage(
    z = Wc, y = Y, n_folds = K, family="binomial")
  end.time = Sys.time()
  lrl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  lrl.Yhat.train = as.numeric(Wc %*% lrl$beta_min)
  lrl.AUC.train = pROC::roc(
    Y, lrl.Yhat.train, levels = c(0, 1), direction = "<")$auc
  # get prediction error on test set
  Wc.test = scale(log(X.test), center = TRUE, scale = FALSE)
  lrl.Yhat.test = as.numeric(Wc.test %*% lrl$beta_min)
  lrl.AUC.test = pROC::roc(
    Y.test, lrl.Yhat.test, levels = c(0, 1), direction = "<")$auc
  
  # beta estimation accuracy, selection accuracy #
  lrl.metrics = getMetricsBM(
    betahat = lrl$beta_min, # don't back-scale bc only centered X (didn't scale)
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  lrl.metrics = c(
    PEtr = lrl.AUC.train, PEte = lrl.AUC.test, lrl.metrics)
  
  saveRDS(c(
    lrl.metrics,
    "betasparsity" = bspars,
    "logratios" = NA,
    "time" = lrl.timing,
    "adhoc" = NA
  ),
  paste0(output_dir, "/lrlasso_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


