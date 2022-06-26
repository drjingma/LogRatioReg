# Purpose: compare slr to other methods on data sets
# Date: 6/19/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs_mse"

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
# 80/20 train/test splits #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSplits
) %dorng% {
  # rm(list=ls())
  library(mvtnorm)
  
  library(Matrix)
  library(glmnet)
  
  source("RCode/func_libs.R")
  source("slr_analyses/Functions/slrs.R")
  source("slr_analyses/Functions/codalasso.R")
  source("slr_analyses/Functions/util.R")
  
  getF1 = function(y, yhat){
    TPR = sum((yhat > 0) & (y == 1)) / sum(y == 1) # TPR (recall)
    prec = sum((yhat > 0) & (y == 1)) / sum(yhat > 0) # precision
    F1 = 2 / (1/TPR + 1/prec) # f1 is harmonic mean of precision & recall
    return(F1)
  }
  
  # tuning parameter settings
  K = 10
  scaling = TRUE
  
  file.end = paste0(
    "_sim", b,
    "_HIV", 
    "_gbm",
    ".rds")
  
  ##############################################################################
  # HIV: one of the HIV data sets in selbal package -- has binary response
  #   n = 155 samples, 
  #   p = 60 taxa (counts for microbial taxa at genus level), 
  #   1 covariate (MSM), 
  #   1 response (HIV_Status - binary)
  W = selbal::HIV[, 1:60]
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
  #   Following Gordon-Rodriguez et al. 2022, fit each method on 20 random 80/20
  #     train/test splits, sampled with stratification by case-control.
  numObs = nrow(X_gbm)
  inputDim = ncol(X_gbm)
  # stratified sampling
  trainIdx = 1:numObs
  caseIdx = sample(cut(1:sum(Y2), breaks=5, labels=F))
  controlIdx = sample(cut(1:sum(1 - Y2), breaks=5, labels=F))
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
    XTr = XTr, YTr = YTr, Y2Tr = Y2Tr,
    XTe = XTe, YTe = YTe, Y2Te = Y2Te
  ),
  paste0(output_dir, "/data", file.end))
  
  ##############################################################################
  # fit methods
  ##############################################################################
  
  # classo #####################################################################
  start.time = Sys.time()
  classo = codalasso(XTr, Y2Tr, numFolds = K)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get prediction error on test set
  classo.Yhat.test = predict(classo, XTe) # before sigmoid

  cl.metrics = c(
    acc = mean((classo.Yhat.test > 0) == Y2Te),
    auc = pROC::roc(
      Y2Te, classo.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(abs(classo$cll$betas[-1]) > 10e-8) / p,
    f1 = getF1(Y2Te, classo.Yhat.test > 0),
    time = cl.timing
  )

  saveRDS(
    cl.metrics,
    paste0(output_dir, "/classo_metrics", file.end))
  
  # slr - spectral clustering ##################################################
  start.time = Sys.time()
  slrspeccv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "mse",
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrspec = slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    threshold = slrspeccv$threshold[slrspeccv$index["1se",]])
  end.time = Sys.time()
  slrspec.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get SBP
  slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrspec.fullSBP) = colnames(XTr)
  slrspec.fullSBP[match(
    names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp
  
  # get prediction error on test set
  slrspec.Yhat.test = predict(
    slrspec$fit,
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrspec.fullSBP)),
    type = "response")
  
  slrspec.metrics = c(
    acc = mean((slrspec.Yhat.test > 0.5) == Y2Te),
    auc = pROC::roc(
      Y2Te, slrspec.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(slrspec.fullSBP > 0) / p,
    f1 = getF1(Y2Te, slrspec.Yhat.test > 0.5),
    time = slrspec.timing
  )
  
  saveRDS(
    slrspec.metrics,
    paste0(output_dir, "/slr_spectral_metrics", file.end))
  
  # slr - hierarchical clustering ##############################################
  start.time = Sys.time()
  slrhiercv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "mse",
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrhier = slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    threshold = slrhiercv$threshold[slrhiercv$index["1se",]])
  end.time = Sys.time()
  slrhier.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get SBP
  slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrhier.fullSBP) = colnames(XTr)
  slrhier.fullSBP[match(
    names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
  
  # get prediction error on test set
  slrhier.Yhat.test = predict(
    slrhier$fit,
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrhier.fullSBP)),
    type = "response")
  
  slrhier.metrics = c(
    acc = mean((slrhier.Yhat.test > 0.5) == Y2Te),
    auc = pROC::roc(
      Y2Te, slrhier.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(slrhier.fullSBP > 0) / p,
    f1 = getF1(Y2Te, slrhier.Yhat.test > 0.5),
    time = slrhier.timing
  )
  
  saveRDS(
    slrhier.metrics,
    paste0(output_dir, "/slr_hierarchical_metrics", file.end))
  
  # selbal #####################################################################
  start.time = Sys.time()
  slbl0 = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K)
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
    newdata = data.frame(V1 = balance::balance.fromSBP(
      x = XTe, y = slbl0.coefs$sbp)),
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

  # selbal with extra covariate ################################################
  start.time = Sys.time()
  slbl1 = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K, covar = covarTr)
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

  # codacore ###################################################################
  start.time = Sys.time()
  codacore0 = codacore::codacore(
    x = XTr, y = Y2Tr, logRatioType = "ILR",
    objective = "binary classification", cvParams = list(numFolds = K))
  end.time = Sys.time()
  codacore0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get prediction error on test set and gamma-hat
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
    codacore0.Yhat.test = predict(codacore0, XTe) # before sigmoid
  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0_coeffs = c()
    codacore0model = stats::glm(Y2Tr ~ 1, family = "binomial")
    codacore0.betahat = rep(0, p)
    codacore0.Yhat.test = predict(codacore0model, XTe) # before sigmoid
  }

  codacore0.metrics = c(
    acc = mean((codacore0.Yhat.test > 0) == Y2Te),
    auc = pROC::roc(
      Y2Te, codacore0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(abs(codacore0.betahat) > 10e-8) / p,
    f1 = getF1(Y2Te, codacore0.Yhat.test > 0),
    time = codacore0.timing
  )

  saveRDS(
    codacore0.metrics,
    paste0(output_dir, "/codacore_metrics", file.end))

  # log-ratio lasso ############################################################
  library(logratiolasso)
  source("slr_analyses/Functions/logratiolasso.R")
  WTr.c = scale(log(XTr), center = TRUE, scale = FALSE)

  start.time = Sys.time()
  lrl <- cv_two_stage(
    z = WTr.c, y = Y2Tr, n_folds = K, family="binomial")
  end.time = Sys.time()
  lrl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # get prediction error on test set
  WTe.c = scale(log(XTe), center = TRUE, scale = FALSE)
  lrl.Yhat.test = as.numeric(WTe.c %*% lrl$beta_min) # before sigmoid

  lrl.metrics = c(
    acc = mean((lrl.Yhat.test > 0) == Y2Te),
    auc = pROC::roc(
      Y2Te, lrl.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(abs(lrl$beta_min) > 10e-8) / p,
    f1 = getF1(Y2Te, lrl.Yhat.test > 0),
    time = lrl.timing
  )

  saveRDS(
    lrl.metrics,
    paste0(output_dir, "/lrlasso_metrics", file.end))
  
}
