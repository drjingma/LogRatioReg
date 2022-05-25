# Purpose: compare slr to selbal on data sets
# Date: 5/25/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_binary"

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
  
  library(balance)
  library(selbal)
  
  library(pROC)
  
  source("RCode/func_libs.R")
  source("Kristyn/Functions/slr.R")
  source("Kristyn/Functions/slrscreen.R")
  source("Kristyn/Functions/codalasso.R")
  source("Kristyn/Functions/util.R")
  
  # tuning parameter settings
  K = 10
  scaling = TRUE
  
  file.end = paste0(
    "_sim", b,
    "_Crohns", 
    "_gbm",
    ".rds")
  
  ##############################################################################
  # Crohn: a data set in selbal package
  #   n = 975 samples, 
  #   p = 48 taxa (counts for microbial taxa at genus level), 
  #   1 response (y - binary)
  W = selbal::Crohn[, 1:48]
  X = sweep(W, 1, rowSums(W), FUN='/')
  Y = selbal::Crohn[, 49]
  Y2 = ifelse(Y == "CD", 1, 0)
  
  ##############################################################################
  # 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
  X_gbm = cmultRepl2(W, zero.rep = "bayes")
  
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
  XTe = X_gbm[trainIdx == 1,]
  YTe = Y[trainIdx == 1]
  Y2Te = Y2[trainIdx == 1]
  
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
  classo.Yhat.test = predict(classo, XTe)
  
  cl.metrics = c(
    acc = mean((classo.Yhat.test > 0) == Y2Te),
    auc = roc(
      Y2Te, classo.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    time = cl.timing
  )
  
  saveRDS(
    cl.metrics,
    paste0(output_dir, "/classo_metrics", file.end))
  
  # slr -- alpha = 0.01 ########################################################
  start.time = Sys.time()
  slr0.01 = slr(x = XTr, y = Y2Tr, alpha = 0.01, classification = TRUE)
  end.time = Sys.time()
  slr0.01.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get prediction error on test set
  slr0.01.Yhat.test = predict(
    slr0.01$model, 
    data.frame(V1 = balance::balance.fromSBP(x = XTe, y = slr0.01$sbp)), 
    type = "response")
  
  slr0.01.metrics = c(
    acc = mean((slr0.01.Yhat.test > 0.5) == Y2Te),
    auc = roc(
      Y2Te, slr0.01.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    time = slr0.01.timing
  )
  
  saveRDS(
    slr0.01.metrics,
    paste0(output_dir, "/slr_alpha0.01_metrics", file.end))
  
  # slr -- alpha = 0.05 ########################################################
  start.time = Sys.time()
  slr0.05 = slr(x = XTr, y = Y2Tr, alpha = 0.05, classification = TRUE)
  end.time = Sys.time()
  slr0.05.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get prediction error on test set
  slr0.05.Yhat.test = predict(
    slr0.05$model, 
    data.frame(V1 = balance::balance.fromSBP(x = XTe, y = slr0.05$sbp)), 
    type = "response")
  
  slr0.05.metrics = c(
    acc = mean((slr0.05.Yhat.test > 0.5) == Y2Te),
    auc = roc(
      Y2Te, slr0.05.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    time = slr0.05.timing
  )
  
  saveRDS(
    slr0.05.metrics,
    paste0(output_dir, "/slr_alpha0.05_metrics", file.end))
  
  # slr -- screen ##############################################################
  start.time = Sys.time()
  slrscreen0cv = cv.slr.screen(
    x = XTr, y = Y2Tr, method = "wald", 
    response.type = "binary", s0.perc = 0, zeta = 0, 
    nfolds = K, type.measure = "mse", 
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrscreen0 = slr.screen(
    x = XTr, y = Y2Tr, method = "wald", 
    response.type = "binary", s0.perc = 0, zeta = 0, 
    threshold = slrscreen0cv$threshold[slrscreen0cv$index["1se",]])
  end.time = Sys.time()
  slrscreen0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slrscreen0.fullSBP = matrix(0, nrow = ncol(XTr), ncol = 1)
  rownames(slrscreen0.fullSBP) = colnames(XTr)
  slrscreen0.fullSBP[match(
    names(slrscreen0$sbp), rownames(slrscreen0.fullSBP))] = slrscreen0$sbp
  
  # get prediction error on test set
  slrscreen0.Yhat.test = predict(
    slrscreen0$fit, 
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrscreen0.fullSBP)), 
    type = "response")
  
  slrscreen0.metrics = c(
    acc = mean((slrscreen0.Yhat.test > 0.5) == Y2Te),
    auc = roc(
      Y2Te, slrscreen0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    time = slrscreen0.timing
  )
  
  saveRDS(
    slrscreen0.metrics,
    paste0(output_dir, "/slrscreen_metrics", file.end))
  
  # selbal #####################################################################
  start.time = Sys.time()
  slbl = selbal.cv(x = XTr, y = YTr, n.fold = K)
  end.time = Sys.time()
  slbl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  slbl.coefs = getCoefsSelbal(
    X = XTr, y = YTr, selbal.fit = slbl, classification = TRUE, 
    check = TRUE)
  
  # get prediction error on test set
  slbl.Yhat.test = predict.glm(
    slbl$glm, 
    newdata = data.frame(V1 = balance::balance.fromSBP(
      x = XTe, y = slbl.coefs$sbp)), 
    type = "response")
  
  slbl.metrics = c(
    acc = mean((slbl.Yhat.test < 0.5) == Y2Te),
    auc = roc(
      YTe, slbl.Yhat.test, levels = levels(YTe), direction = "<")$auc,
    time = slbl.timing
  )
  
  saveRDS(
    slbl.metrics,
    paste0(output_dir, "/selbal_metrics", file.end))
  
  # codacore ###################################################################
  library(codacore)
  start.time = Sys.time()
  codacore0 = codacore(
    x = XTr, y = Y2Tr, logRatioType = "ILR", 
    objective = "binary classification", cvParams = list(numFolds = K))
  end.time = Sys.time()
  codacore0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  if(length(codacore0$ensemble) > 0){ # at least 1 log-ratio found
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0, XTe)
    
  } else{
    print(paste0("sim ", i, " -- codacore has no log-ratios"))
    codacore0model = stats::glm(Y ~ 1, family = "binomial")
    
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0model, data.frame(XTe))
  }
  
  codacore0.metrics = c(
    acc = mean((codacore0.Yhat.test > 0.5) == Y2Te),
    auc = roc(
        Y2Te, codacore0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    time = codacore0.timing
  )
  
  saveRDS(
    codacore0.metrics,
    paste0(output_dir, "/codacore_metrics", file.end))
  
}