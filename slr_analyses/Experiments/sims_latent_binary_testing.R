# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 7/19/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Experiments/outputs/metrics_binary"

################################################################################
# Simulations #
################################################################################
seed = 1234

for(b in 1:100){
  # print(b)
  # rm(list=ls())
  set.seed(seed + b)
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
  sigma.settings = "latentVarModel_binary"
  n = 100
  p = 30
  K = 10
  nlam = 100
  scaling = TRUE
  sigma_eps2 = 0.1
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  # SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
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
  if(file.exists(paste0(output_dir, "/data", file.end))){
    data.tmp = readRDS(paste0(output_dir, "/data", file.end))
    X = data.tmp$X
    Y = data.tmp$Y
    X.test = data.tmp$X.test
    Y.test = data.tmp$Y.test
    SBP.true = data.tmp$SBP.true
    beta.true = data.tmp$beta.true
    non0.beta = data.tmp$non0.beta
    bspars = sum(non0.beta)
  } else{
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
  }
  
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
  slrspeccv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "auc", 
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrspec = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0, 
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
  # prediction errors
  # get prediction error on training set
  slrspec.Yhat.train = predict(
    slrspec$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = X, y = slrspec.fullSBP)), 
    type = "response")
  slrspec.AUC.train = pROC::roc(
    Y, slrspec.Yhat.train, levels = c(0, 1), direction = "<")$auc
  # get prediction error on test set
  slrspec.Yhat.test = predict(
    slrspec$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = X.test, y = slrspec.fullSBP)), 
    type = "response")
  slrspec.AUC.test = pROC::roc(
    Y.test, slrspec.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  slrspec.metrics = getMetricsBM(
    thetahat0 = slrspec.coefs$a0, thetahat = slrspec.coefs$bm.coefs,
    betahat = slrspec.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true,
    metrics = c("betaestimation", "selection"), classification = TRUE)
  slrspec.metrics = c(
    AUCtr = slrspec.AUC.train, AUCte = slrspec.AUC.test, 
    slrspec.metrics)
  
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
  slrhiercv = cv.slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "auc", 
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrhier = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0, 
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
  # prediction errors
  # get prediction error on training set
  slrhier.Yhat.train = predict(
    slrhier$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = X, y = slrhier.fullSBP)), 
    type = "response")
  slrhier.AUC.train = pROC::roc(
    Y, slrhier.Yhat.train, levels = c(0, 1), direction = "<")$auc
  # get prediction error on test set
  slrhier.Yhat.test = predict(
    slrhier$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = X.test, y = slrhier.fullSBP)), 
    type = "response")
  slrhier.AUC.test = pROC::roc(
    Y.test, slrhier.Yhat.test, levels = c(0, 1), direction = "<")$auc
  # beta estimation accuracy, selection accuracy #
  slrhier.metrics = getMetricsBM(
    thetahat0 = slrhier.coefs$a0, thetahat = slrhier.coefs$bm.coefs,
    betahat = slrhier.coefs$llc.coefs,
    true.sbp = SBP.true, non0.true.beta = non0.beta,
    true.beta = beta.true,
    metrics = c("betaestimation", "selection"), classification = TRUE)
  slrhier.metrics = c(
    AUCtr = slrhier.AUC.train, AUCte = slrhier.AUC.test, 
    slrhier.metrics)
  
  ### fin ###
  print(paste0("sim ", b, " completed successfully."))
}


