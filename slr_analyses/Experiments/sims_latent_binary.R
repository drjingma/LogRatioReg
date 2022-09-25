# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 7/19/2022
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
  print(b)
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
  settings.name = "BinaryResponse"
  hparam = "1se"
  n = 100
  p = 30
  K = 10
  nlam = 100
  scaling = TRUE
  sigma_x = 0.1
  # SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 6 # 6
  c.value = 1 # a1 = c.value / k+ or c.value / k- or 0
  a0 = 0 # 0
  ulimit = 0.5
  
  file.end = paste0(
    "_", settings.name,
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_hparam", hparam,
    "_dim", n, "x", p, 
    "_ulimit", ulimit,
    # "_noisey", sigma_y, 
    "_noisex", sigma_x, 
    "_b0", b0, 
    "_b1", b1, 
    "_a0", a0, 
    "_c", c.value,
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # if(file.exists(paste0(output_dir, "/data", file.end))){
  #   data.tmp = readRDS(paste0(output_dir, "/data", file.end))
  #   X = data.tmp$X
  #   Y = data.tmp$Y
  #   X.test = data.tmp$X.test
  #   Y.test = data.tmp$Y.test
  #   SBP.true = data.tmp$SBP.true
  #   llc.coefs.true = data.tmp$llc.coefs.true
  #   llc.coefs.non0 = data.tmp$llc.coefs.non0
  # } else{
    # get latent variable
    U.all = matrix(runif(min = -ulimit, max = ulimit, 2 * n), ncol = 1)
    # simulate y from latent variable
    y.all = rbinom(n = 2 * n, size = 1, p = as.vector(sigmoid(b0 + b1 * U.all)))
    # simulate X: 
    epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_x
    a1 = c.value * ilrtrans.true$ilr.trans.unscaled[-p] 
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
    
    # about linear log-contrast models' coefficients
    llc.coefs.non0 = as.vector(SBP.true != 0)
    # solve for beta
    theta.value = c.value / ilrtrans.true$const
    c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
    llc.coefs.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
      as.vector(ilrtrans.true$ilr.trans)
    
    saveRDS(list(
      X = X, Y = Y, X.test = X.test, Y.test = Y.test, 
      SBP.true = SBP.true, llc.coefs.true = llc.coefs.true, 
      llc.coefs.non0 = llc.coefs.non0
    ),
    paste0(output_dir, "/data", file.end))
  # }

  ##############################################################################
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  start.time = Sys.time()
  if(hparam == "min"){
    classo = codalasso(X, Y, numFolds = K, gamma = 0)
  } else if(hparam == "1se"){
    classo = codalasso(X, Y, numFolds = K, gamma = 1)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  cl.betahat = classo$cll$betas[-1]

  # compute metrics on the selected model #
  # prediction error
  cl.Yhat.test = predict(classo, X.test)
  cl.AUC.test = pROC::roc(
    Y.test, cl.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # estimation accuracy, selection accuracy #
  cl.metrics = getMetricsLLC(
    est.llc.coefs = cl.betahat,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))

  saveRDS(c(
    "auc" = cl.AUC.test,
    cl.metrics,
    "logratios" = 0,
    "time" = cl.timing
  ),
  paste0(output_dir, "/classo_metrics", file.end))
  
  ##############################################################################
  # slr
  #   screening.method = "wald"
  #   cluster.method = "spectral"
  #   response.type = "binary"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "auc"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slrspec1cv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "auc",
     scale = scaling, trace.it = FALSE)
  if(hparam == "min"){
    slrspec1 = slr(
      x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrspec1cv$threshold[slrspec1cv$index["min",]], 
      positive.slope = TRUE)
  } else if(hparam == "1se"){
    slrspec1 = slr(
      x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrspec1cv$threshold[slrspec1cv$index["1se",]], 
      positive.slope = TRUE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slrspec1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrspec1.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrspec1.fullSBP) = colnames(X)
  slrspec1.fullSBP[match(
    names(slrspec1$sbp), rownames(slrspec1.fullSBP))] = slrspec1$sbp

  slrspec1.coefs = getCoefsBM(
    coefs = coefficients(slrspec1$fit), sbp = slrspec1.fullSBP)

  # compute metrics on the selected model #
  # prediction error
  slrspec1.Yhat.test = predict(
    slrspec1$fit,
    data.frame(balance = slr.fromContrast(X.test, slrspec1.fullSBP)),
    type = "response")
  slrspec1.AUC.test = pROC::roc(
    Y.test, slrspec1.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  slrspec1.metrics = getMetricsBM(
    est.llc.coefs = slrspec1.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))

  saveRDS(c(
    "auc" = slrspec1.AUC.test,
    slrspec1.metrics,
    "logratios" = sum(slrspec1.coefs$bm.coefs != 0),
    "time" = slrspec1.timing
  ),
  paste0(output_dir, "/slr_spectral_metrics", file.end))

  ##############################################################################
  # slr
  #   screening.method = "wald"
  #   cluster.method = "hierarchical"
  #   response.type = "binary"
  #   s0.perc = 0
  #   zeta = 0
  #   type.measure = "auc"
  # -- fits a balance regression model with one balance
  ##############################################################################
  start.time = Sys.time()
  slrhier1cv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "auc",
     scale = scaling, trace.it = FALSE)
  if(hparam == "min"){
    slrhier1 = slr(
      x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrhier1cv$threshold[slrhier1cv$index["min",]], 
      positive.slope = TRUE)
  } else if(hparam == "1se"){
    slrhier1 = slr(
      x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrhier1cv$threshold[slrhier1cv$index["1se",]], 
      positive.slope = TRUE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slrhier1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrhier1.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrhier1.fullSBP) = colnames(X)
  slrhier1.fullSBP[match(
    names(slrhier1$sbp), rownames(slrhier1.fullSBP))] = slrhier1$sbp

  slrhier1.coefs = getCoefsBM(
    coefs = coefficients(slrhier1$fit), sbp = slrhier1.fullSBP)

  # compute metrics on the selected model #
  # prediction error
  slrhier1.Yhat.test = predict(
    slrhier1$fit,
    data.frame(balance = slr.fromContrast(X.test, slrhier1.fullSBP)),
    type = "response")
  slrhier1.AUC.test = pROC::roc(
    Y.test, slrhier1.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  slrhier1.metrics = getMetricsBM(
    est.llc.coefs = slrhier1.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true,
    metrics = c("estimation", "selection"))

  saveRDS(c(
    "auc" = slrhier1.AUC.test,
    slrhier1.metrics,
    "logratios" = sum(slrhier1.coefs$bm.coefs != 0),
    "time" = slrhier1.timing
  ),
  paste0(output_dir, "/slr_hierarchical_metrics", file.end))

  ##############################################################################
  # selbal method (a balance regression method)
  # -- fits a balance regression model with one balance
  ##############################################################################
  library(selbal) # masks stats::cor()
  slbl.data = getSelbalData(
    X = X, y = Y, classification = TRUE, levels = c(0, 1), labels = c(0, 1))

  start.time = Sys.time()
  if(hparam == "min"){
    slbl = selbal::selbal.cv(
      x = slbl.data$X, y = slbl.data$y, n.fold = K, opt.cri = "min")
  } else if(hparam == "1se"){
    slbl = selbal::selbal.cv(
      x = slbl.data$X, y = slbl.data$y, n.fold = K, opt.cri = "1se")
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slbl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slbl.coefs = getCoefsSelbal(
    X = slbl.data$X, y = slbl.data$y, selbal.fit = slbl, classification = TRUE,
    check = TRUE)

  # compute metrics on the selected model #
  # prediction error
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
    est.llc.coefs = slbl.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true, metrics = c("estimation", "selection"))

  saveRDS(c(
    "auc" = slbl.AUC.test,
    slbl.metrics,
    "logratios" = sum(slbl.coefs$bm.coefs != 0),
    "time" = slbl.timing
  ),
  paste0(output_dir, "/selbal_metrics", file.end))

  ##############################################################################
  # CoDaCoRe (a balance regression method)
  # -- fits a balance regression model with possibly multiple balances
  ##############################################################################
  library(codacore)

  if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
    reticulate::use_condaenv("anaconda3")
  }

  start.time = Sys.time()
  if(hparam == "min"){
    codacore0 = codacore::codacore(
      x = X, y = Y, logRatioType = "ILR",
      objective = "binary classification", cvParams = list(numFolds = K), 
      lambda = 0) 
  } else if(hparam == "1se"){
    codacore0 = codacore::codacore(
      x = X, y = Y, logRatioType = "ILR",
      objective = "binary classification", cvParams = list(numFolds = K), 
      lambda = 1) 
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
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
    # prediction error
    codacore0.Yhat.test = predict(codacore0, X.test)

  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0_coeffs = c()
    codacore0model = stats::glm(Y ~ 1, family = "binomial")
    codacore0.betahat = rep(0, p)

    # compute metrics on the selected model #
    # prediction error
    codacore0.Yhat.test = predict(codacore0model, X.test)
  }
  codacore0.AUC.test = pROC::roc(
    Y.test, codacore0.Yhat.test, levels = c(0, 1), direction = "<")$auc

  # beta estimation accuracy, selection accuracy #
  codacore0.metrics = getMetricsBM(
    est.llc.coefs = codacore0.betahat,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true, metrics = c("estimation", "selection"))

  saveRDS(c(
    "auc" = codacore0.AUC.test,
    codacore0.metrics,
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
  if(hparam == "min"){
    lrl <- cv_two_stage(
      z = Wc, y = Y, n_folds = K, family="binomial", gamma = 0)
    lrl.betahat = lrl$beta_min
  } else if(hparam == "1se"){
    lrl <- cv_two_stage(
      z = Wc, y = Y, n_folds = K, family="binomial", gamma = 1)
    lrl.betahat = lrl$beta_gammase
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  lrl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # compute metrics on the selected model #
  # prediction error
  Wc.test = scale(log(X.test), center = TRUE, scale = FALSE)
  lrl.Yhat.test = as.numeric(Wc.test %*% lrl.betahat)
  lrl.AUC.test = pROC::roc(
    Y.test, lrl.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  lrl.metrics = getMetricsBM(
    est.llc.coefs = lrl.betahat,
    true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
    true.llc.coefs = llc.coefs.true, metrics = c("estimation", "selection"))

  saveRDS(c(
    "auc" = lrl.AUC.test,
    lrl.metrics,
    "logratios" = NA,
    "time" = lrl.timing
  ),
  paste0(output_dir, "/lrlasso_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
  # print(paste0("sim ", b, " completed successfully."))
}


