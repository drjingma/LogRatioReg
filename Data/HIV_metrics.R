# Purpose: compare slr to other methods on data sets
# Date: 6/28/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "Data/outputs_metrics"

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
numSplits = 20

################################################################################
# 70/30 train/test splits with sparsity filtering at 80%
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSplits
) %dorng% {
  # rm(list=ls())
  library(mvtnorm)
  
  library(Matrix)
  library(glmnet)
  
  source("Functions/slrs.R")
  source("Functions/util.R")
  
  getF1 = function(y, yhat){
    TPR = sum((yhat > 0) & (y == 1)) / sum(y == 1) # TPR (recall)
    prec = sum((yhat > 0) & (y == 1)) / sum(yhat > 0) # precision
    F1 = 2 / (1/TPR + 1/prec) # f1 is harmonic mean of precision & recall
    return(F1)
  }
  
  # tuning parameter settings
  hparam = "1se"
  K = 10
  scaling = TRUE
  
  filter.perc = 0.8 # 0.8, 1
  split.perc = 0.7 # 0.7, 0.8
  
  file.end = paste0(
    "_HIV", 
    "_split", split.perc, 
    "_filter", filter.perc, 
    "_hparam", hparam,
    "_gbm",
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # HIV: one of the HIV data sets in selbal package -- has binary response
  #   n = 155 samples, 
  #   p = 60 taxa (counts for microbial taxa at genus level), 
  #   1 covariate (MSM), 
  #   1 response (HIV_Status - binary)
  W = selbal::HIV[, 1:60]
  W.origin <- W
  W <- W[,apply(W==0,2,mean)<filter.perc]
  covar = data.frame(MSM = selbal::HIV[, 61])
  X = sweep(W, 1, rowSums(W), FUN='/')
  Y = selbal::HIV[, 62]
  # levels(Y) # (control, case)
  Y2 = ifelse(Y == "Pos", 1, 0)
  p = ncol(W)
  
  ##############################################################################
  # 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
  X_gbm = selbal::cmultRepl2(W, zero.rep = "bayes")
  
  ##############################################################################
  # Train/Test Split
  #   Following Gordon-Rodriguez et al. 2022, fit each method on 20 random
  #     train/test splits, sampled with stratification by case-control.
  numObs = nrow(X_gbm)
  inputDim = ncol(X_gbm)
  if(file.exists(paste0(output_dir, "/data", file.end))){
    data.tmp = readRDS(paste0(output_dir, "/data", file.end))
    XTr = data.tmp$XTr
    XTe = data.tmp$XTe
    YTr = data.tmp$YTr
    YTe = data.tmp$YTe
    Y2Tr = data.tmp$Y2Tr
    Y2Te = data.tmp$Y2Te
    covarTr = data.tmp$covarTr
    covarTe = data.tmp$covarTe
  } else{
    # stratified sampling
    trainIdx = 1:numObs
    # caseIdx = sample(cut(1:sum(Y2), breaks=5, labels=F))
    caseIdx = 1 - 
      (seq(sum(Y2)) %in% sample(1:sum(Y2), size=ceiling(split.perc * sum(Y2))))
    # controlIdx = sample(cut(1:sum(1 - Y2), breaks=5, labels=F))
    controlIdx = 1 - 
      (seq(sum(1 - Y2)) %in% 
         sample(1:sum(1 - Y2), size = ceiling(split.perc * sum(1 - Y2))))
    trainIdx[Y2 == 1] = caseIdx
    trainIdx[Y2 == 0] = controlIdx
    XTr = X_gbm[trainIdx != 1,]
    YTr = Y[trainIdx != 1]
    Y2Tr = Y2[trainIdx != 1]
    covarTr = data.frame(MSM = covar[trainIdx != 1, ])
    XTe = X_gbm[trainIdx == 1,]
    YTe = Y[trainIdx == 1]
    Y2Te = Y2[trainIdx == 1]
    covarTe = data.frame(MSM = covar[trainIdx == 1,])
    
    saveRDS(list(
      XTr = XTr, YTr = YTr, Y2Tr = Y2Tr, covarTr = covarTr,
      XTe = XTe, YTe = YTe, Y2Te = Y2Te, covarTe = covarTe
    ),
    paste0(output_dir, "/data", file.end))
  }
  
  ##############################################################################
  # fit methods
  ##############################################################################
  
  # classo #####################################################################
  source("Functions/codalasso.R")
  start.time = Sys.time()
  if(hparam == "min"){
    classo1 = codalasso(
      XTr, Y2Tr, numFolds = K, gamma = 0, type.measure = "AUC", 
      stratify = FALSE)
  } else if(hparam == "1se"){
    classo1 = codalasso(
      XTr, Y2Tr, numFolds = K, gamma = 1, type.measure = "AUC", 
      stratify = FALSE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  cl1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get prediction error on test set
  classo1.Yhat.test = predict(classo1, XTe) # before sigmoid

  cl1.metrics = c(
    acc = mean((classo1.Yhat.test > 0) == Y2Te),
    auc = pROC::roc(
      Y2Te, classo1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(abs(classo1$cll$betas[-1]) > 10e-8) / p,
    f1 = getF1(Y2Te, classo1.Yhat.test > 0),
    time = cl1.timing
  )

  saveRDS(
    cl1.metrics,
    paste0(output_dir, "/classo_metrics", file.end))

  # slr - spectral clustering with auc #########################################
  start.time = Sys.time()
  slrspec1cv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "auc",
    scale = scaling, trace.it = FALSE)
  if(hparam == "min"){
    slrspec1 = slr(
      x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrspec1cv$threshold[slrspec1cv$index["min",]],
      positive.slope = TRUE)
  } else if(hparam == "1se"){
    slrspec1 = slr(
      x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrspec1cv$threshold[slrspec1cv$index["1se",]],
      positive.slope = TRUE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slrspec1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get SBP
  slrspec1.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrspec1.fullSBP) = colnames(XTr)
  slrspec1.fullSBP[match(
    names(slrspec1$sbp), rownames(slrspec1.fullSBP))] = slrspec1$sbp

  # get prediction error on test set
  slrspec1.Yhat.test = predict(
    slrspec1$fit,
    data.frame(balance = slr.fromContrast(XTe, slrspec1.fullSBP)),
    type = "response")

  slrspec1.metrics = c(
    acc = mean((slrspec1.Yhat.test > 0.5) == Y2Te),
    auc = pROC::roc(
      Y2Te, slrspec1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(slrspec1.fullSBP > 0) / p,
    f1 = getF1(Y2Te, slrspec1.Yhat.test > 0.5),
    time = slrspec1.timing
  )

  saveRDS(
    slrspec1.metrics,
    paste0(output_dir, "/slr_spectral_metrics", file.end))

  if(!all(slrspec1.fullSBP == 0) & slrspec1$theta[2] < 0){
    slrspec1.fullSBP = -slrspec1.fullSBP
  }
  saveRDS(
    slrspec1.fullSBP,
    paste0(output_dir, "/slr_spectral_sbp", file.end)
  )

  # slr - hierarchical clustering with auc #####################################
  start.time = Sys.time()
  slrhier1cv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "auc",
    scale = scaling, trace.it = FALSE)
  if(hparam == "min"){
    slrhier1 = slr(
      x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrhier1cv$threshold[slrhier1cv$index["min",]],
      positive.slope = TRUE)
  } else if(hparam == "1se"){
    slrhier1 = slr(
      x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
      response.type = "binary", s0.perc = 0, zeta = 0,
      threshold = slrhier1cv$threshold[slrhier1cv$index["1se",]],
      positive.slope = TRUE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slrhier1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get SBP
  slrhier1.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrhier1.fullSBP) = colnames(XTr)
  slrhier1.fullSBP[match(
    names(slrhier1$sbp), rownames(slrhier1.fullSBP))] = slrhier1$sbp

  # get prediction error on test set
  slrhier1.Yhat.test = predict(
    slrhier1$fit,
    data.frame(balance = slr.fromContrast(XTe, slrhier1.fullSBP)),
    type = "response")

  slrhier1.metrics = c(
    acc = mean((slrhier1.Yhat.test > 0.5) == Y2Te),
    auc = pROC::roc(
      Y2Te, slrhier1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(slrhier1.fullSBP > 0) / p,
    f1 = getF1(Y2Te, slrhier1.Yhat.test > 0.5),
    time = slrhier1.timing
  )

  saveRDS(
    slrhier1.metrics,
    paste0(output_dir, "/slr_hierarchical_metrics", file.end))

  if(!all(slrhier1.fullSBP == 0) & slrhier1$theta[2] < 0){
    slrhier1.fullSBP = -slrhier1.fullSBP
  }
  saveRDS(
    slrhier1.fullSBP,
    paste0(output_dir, "/slr_hierarchical_sbp", file.end)
  )

  # selbal #####################################################################
  start.time = Sys.time()
  if(hparam == "min"){
    slbl0 = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K, opt.cri = "min")
  } else if(hparam == "1se"){
    slbl0 = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K, opt.cri = "1se")
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slbl0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get theta-hat and gamma-hat
  slbl0.coefs = getCoefsSelbal(
    X = XTr, y = YTr, selbal.fit = slbl0, classification = TRUE,
    check = TRUE)

  # get prediction error on test set
  slbl0.Yhat.test = predict.glm(
    slbl0$glm,
    newdata = data.frame(V1 = balance::balance.fromSBP(XTe, slbl0.coefs$sbp)),
    type = "response")

  slbl0.metrics = c(
    acc = mean((slbl0.Yhat.test > 0.5) == Y2Te),
    # < 0.5 bc order of levels = c(case, control) instead of c(control, case)
    auc = pROC::roc(
      YTe, slbl0.Yhat.test, levels = c("Neg", "Pos"), direction = "<")$auc,
    percselected = sum(slbl0.coefs$sbp > 0) / p,
    f1 = getF1(Y2Te, slbl0.Yhat.test > 0.5),
    time = slbl0.timing
  )

  saveRDS(
    slbl0.metrics,
    paste0(output_dir, "/selbal_metrics", file.end))

  slbl0_sbp = slbl0.coefs$sbp
  if(slbl0$glm$coefficients[2] < 0){
    slbl0_sbp = -slbl0_sbp
  }
  saveRDS(
    slbl0_sbp,
    paste0(output_dir, "/selbal_sbp", file.end)
  )

  # selbal with extra covariate ################################################
  start.time = Sys.time()
  if(hparam == "min"){
    slbl1 = selbal::selbal.cv(
      x = XTr, y = YTr, n.fold = K, covar = covarTr, opt.cri = "min")
  } else if(hparam == "1se"){
    slbl1 = selbal::selbal.cv(
      x = XTr, y = YTr, n.fold = K, covar = covarTr, opt.cri = "1se")
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slbl1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get theta-hat and gamma-hat
  slbl1.coefs = getCoefsSelbal(
    X = XTr, y = YTr, selbal.fit = slbl1, classification = TRUE,
    check = TRUE, covar = covarTr)

  # get prediction error on test set
  slbl1.Yhat.test = predict.glm(
    slbl1$glm,
    newdata = data.frame(
      V1 = balance::balance.fromSBP(x = XTe, y = slbl1.coefs$sbp),
      MSM = covarTe
    ),
    type = "response")

  slbl1.metrics = c(
    acc = mean((slbl1.Yhat.test > 0.5) == Y2Te),
    # < 0.5 bc order of levels = c(case, control) instead of c(control, case)
    auc = pROC::roc(
      YTe, slbl1.Yhat.test, levels = c("Neg", "Pos"), direction = "<")$auc,
    percselected = sum(slbl1.coefs$sbp > 0) / p,
    f1 = getF1(Y2Te, slbl1.Yhat.test > 0.5),
    time = slbl1.timing
  )

  saveRDS(
    slbl1.metrics,
    paste0(output_dir, "/selbal_covar_metrics", file.end))

  slbl1_sbp = slbl1.coefs$sbp
  if(slbl1$glm$coefficients[2] < 0){
    slbl1_sbp = -slbl1_sbp
  }
  saveRDS(
    slbl1_sbp,
    paste0(output_dir, "/selbal_covar_sbp", file.end)
  )
  
  # codacore - 1 balance #######################################################
  library(codacore)
  
  start.time = Sys.time()
  if(hparam == "min"){
    codacore1 = codacore::codacore(
      x = XTr, y = Y2Tr, logRatioType = "ILR",
      objective = "binary classification", cvParams = list(numFolds = K),
      maxBaseLearners = 1,
      lambda = 0)
  } else if(hparam == "1se"){
    codacore1 = codacore::codacore(
      x = XTr, y = Y2Tr, logRatioType = "ILR",
      objective = "binary classification", cvParams = list(numFolds = K),
      maxBaseLearners = 1,
      lambda = 1)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  codacore1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get prediction error on test set and gamma-hat
  if(length(codacore1$ensemble) > 0){ # at least 1 log-ratio found
    codacore1_SBP = matrix(0, nrow = p, ncol = 1)
    codacore1_SBP[codacore1$ensemble[[1]]$hard$numerator, 1] = 1
    codacore1_SBP[codacore1$ensemble[[1]]$hard$denominator, 1] = -1
    codacore1_coeffs = codacore1$ensemble[[1]]$slope
    # codacore1.betahat = getBetaFromCodacore(
    #   SBP_codacore = codacore1_SBP, coeffs_codacore = codacore1_coeffs, p = p)
    names(codacore1_coeffs) = "balance1"
    codacore1.coefs2 = getCoefsBM(
      coefs = codacore1_coeffs,
      sbp = codacore1_SBP)
    codacore1.betahat = codacore1.coefs2$llc.coefs
    codacore1.Yhat.test = predict(codacore1, XTe)
    # adjust codacore_SBP to correspond to positive theta-hats
    if(codacore1_coeffs[1] < 0){
      codacore1_SBP[, 1] = -codacore1_SBP[, 1]
    }
  } else{
    print(paste0("sim ", b, " -- codacore has no log-ratios"))
    codacore1_coeffs = c()
    codacore1_SBP = matrix(0, nrow = p, ncol = 1) 
    codacore1model = stats::glm(YTr ~ 1, family = "gaussian")
    codacore1.betahat = rep(0, p)
    codacore1.Yhat.test = predict(codacore1model, XTe)
  }
  rownames(codacore1_SBP) = colnames(XTr) 
  
  codacore1.metrics = c(
    acc = mean((codacore1.Yhat.test > 0) == Y2Te),
    auc = pROC::roc(
      Y2Te, codacore1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(abs(codacore1.betahat) > 10e-8) / p,
    f1 = getF1(Y2Te, codacore1.Yhat.test > 0),
    time = codacore1.timing
  )
  
  saveRDS(
    codacore1.metrics,
    paste0(output_dir, "/codacore1_metrics", file.end))
  
  saveRDS(
    codacore1_SBP,
    paste0(output_dir, "/codacore1_sbp", file.end)
  )

  # log-ratio lasso ############################################################
  library(logratiolasso)
  source("Functions/logratiolasso.R")
  WTr.c = scale(log(XTr), center = TRUE, scale = FALSE)

  start.time = Sys.time()
  if(hparam == "min"){
    lrl <- cv_two_stage(
      z = WTr.c, y = Y2Tr, n_folds = K, family="binomial", gamma = 0)
    lrl.betahat = lrl$beta_min
  } else if(hparam == "1se"){
    lrl <- cv_two_stage(
      z = WTr.c, y = Y2Tr, n_folds = K, family="binomial", gamma = 1)
    lrl.betahat = lrl$beta_gammase
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  lrl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get prediction error on test set
  WTe.c = scale(log(XTe), center = TRUE, scale = FALSE)
  lrl.Yhat.test = as.numeric(WTe.c %*% lrl.betahat) # before sigmoid

  lrl.metrics = c(
    acc = mean((lrl.Yhat.test > 0) == Y2Te),
    auc = pROC::roc(
      Y2Te, lrl.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(abs(lrl.betahat) > 10e-8) / p,
    f1 = getF1(Y2Te, lrl.Yhat.test > 0),
    time = lrl.timing
  )

  saveRDS(
    lrl.metrics,
    paste0(output_dir, "/lrlasso_metrics", file.end))

}

