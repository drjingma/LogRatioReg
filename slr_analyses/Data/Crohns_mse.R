# Purpose: compare slr to other methods on data sets
# Date: 6/18/2022
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
  
  getF1 = function(y, yhat){ # y and yhat are binary
      TPR = sum((yhat == 1) & (y == 1)) / sum(y == 1) # TPR (recall)
      prec = sum((yhat == 1) & (y == 1)) / sum(yhat  == 1) # precision
    F1 = 2 / (1/TPR + 1/prec) # f1 is harmonic mean of precision & recall
    return(F1)
  }
  
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
  levels(Y) = c("no", "CD") # (control, case)
  Y2 = ifelse(Y == "CD", 1, 0)
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
  if(file.exists(paste0(output_dir, "/data", file.end))){
    data.tmp = readRDS(paste0(output_dir, "/data", file.end))
    XTr = data.tmp$XTr
    XTe = data.tmp$XTe
    YTr = data.tmp$YTr
    YTe = data.tmp$YTe
    Y2Tr = data.tmp$Y2Tr
    Y2Te = data.tmp$Y2Te
  } else{
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
  }
  
  ##############################################################################
  # fit methods
  ##############################################################################
  
  # # classo #####################################################################
  # start.time = Sys.time()
  # classo = codalasso(XTr, Y2Tr, numFolds = K)
  # end.time = Sys.time()
  # cl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get prediction error on test set
  # classo.Yhat.test = predict(classo, XTe) # before sigmoid
  # 
  # cl.metrics = c(
  #   acc = mean((classo.Yhat.test > 0) == Y2Te),
  #   auc = pROC::roc(
  #     Y2Te, classo.Yhat.test, levels = c(0, 1), direction = "<")$auc,
  #   percselected = sum(abs(classo$cll$betas[-1]) > 10e-8) / p,
  #   f1 = getF1(Y2Te, classo.Yhat.test > 0),
  #   time = cl.timing
  # )
  # 
  # saveRDS(
  #   cl.metrics,
  #   paste0(output_dir, "/classo_metrics", file.end))
  
  # slr - spectral clustering - accuracy #######################################
  start.time = Sys.time()
  slrspec0cv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "accuracy",
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrspec0 = slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    threshold = slrspec0cv$threshold[slrspec0cv$index["1se",]])
  end.time = Sys.time()
  slrspec0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get SBP
  slrspec0.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrspec0.fullSBP) = colnames(XTr)
  slrspec0.fullSBP[match(
    names(slrspec0$sbp), rownames(slrspec0.fullSBP))] = slrspec0$sbp
  
  # get prediction error on test set
  slrspec0.Yhat.test = predict(
    slrspec0$fit,
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrspec0.fullSBP)),
    type = "response")
  
  slrspec0.metrics = c(
    acc = mean((slrspec0.Yhat.test > 0.5) == Y2Te),
    auc = pROC::roc(
      Y2Te, slrspec0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(slrspec0.fullSBP > 0) / p,
    f1 = getF1(Y2Te, slrspec0.Yhat.test > 0.5),
    time = slrspec0.timing
  )
  
  saveRDS(
    slrspec0.metrics,
    paste0(output_dir, "/slr_spectral_accuracy_metrics", file.end))
  
  if(slrspec0$theta[2] < 0){
    slrspec0.fullSBP = -slrspec0.fullSBP
  }
  saveRDS(
    slrspec0.fullSBP, 
    paste0(output_dir, "/slr_spectral_accuracy_sbp", file.end)
  )
  
  # slr - hierarchical clustering - accuracy ###################################
  start.time = Sys.time()
  slrhier0cv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "accuracy",
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrhier0 = slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    threshold = slrhier0cv$threshold[slrhier0cv$index["1se",]])
  end.time = Sys.time()
  slrhier0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get SBP
  slrhier0.fullSBP = matrix(0, nrow = p, ncol = 1)
  rownames(slrhier0.fullSBP) = colnames(XTr)
  slrhier0.fullSBP[match(
    names(slrhier0$sbp), rownames(slrhier0.fullSBP))] = slrhier0$sbp
  
  # get prediction error on test set
  slrhier0.Yhat.test = predict(
    slrhier0$fit,
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrhier0.fullSBP)),
    type = "response")
  
  slrhier0.metrics = c(
    acc = mean((slrhier0.Yhat.test > 0.5) == Y2Te),
    auc = pROC::roc(
      Y2Te, slrhier0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    percselected = sum(slrhier0.fullSBP > 0) / p,
    f1 = getF1(Y2Te, slrhier0.Yhat.test > 0.5),
    time = slrhier0.timing
  )
  
  saveRDS(
    slrhier0.metrics,
    paste0(output_dir, "/slr_hierarchical_accuracy_metrics", file.end))
  
  if(slrhier0$theta[2] < 0){
    slrhier0.fullSBP = -slrhier0.fullSBP
  }
  saveRDS(
    slrhier0.fullSBP, 
    paste0(output_dir, "/slr_hierarchical_accuracy_sbp", file.end)
  )
  
  # slr - spectral clustering - auc ############################################
  start.time = Sys.time()
  slrspec1cv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "auc",
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrspec1 = slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "spectral",
    response.type = "binary", s0.perc = 0, zeta = 0,
    threshold = slrspec1cv$threshold[slrspec1cv$index["1se",]])
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
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrspec1.fullSBP)),
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
    paste0(output_dir, "/slr_spectral_auc_metrics", file.end))
  
  if(slrspec1$theta[2] < 0){
    slrspec1.fullSBP = -slrspec1.fullSBP
  }
  saveRDS(
    slrspec1.fullSBP, 
    paste0(output_dir, "/slr_spectral_auc_sbp", file.end)
  )
  
  # slr - hierarchical clustering - auc ########################################
  start.time = Sys.time()
  slrhier1cv = cv.slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    nfolds = K, type.measure = "auc",
    parallel = FALSE, scale = scaling, trace.it = FALSE)
  slrhier1 = slr(
    x = XTr, y = Y2Tr, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "binary", s0.perc = 0, zeta = 0,
    threshold = slrhier1cv$threshold[slrhier1cv$index["1se",]])
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
    data.frame(balance = balance::balance.fromSBP(
      x = XTe, y = slrhier1.fullSBP)),
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
    paste0(output_dir, "/slr_hierarchical_auc_metrics", file.end))
  
  if(slrhier1$theta[2] < 0){
    slrhier1.fullSBP = -slrhier1.fullSBP
  }
  saveRDS(
    slrhier1.fullSBP, 
    paste0(output_dir, "/slr_hierarchical_auc_sbp", file.end)
  )
  
  # # selbal #####################################################################
  # start.time = Sys.time()
  # slbl = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K)
  # end.time = Sys.time()
  # slbl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get theta-hat and gamma-hat
  # slbl.coefs = getCoefsSelbal(
  #   X = XTr, y = YTr, selbal.fit = slbl, classification = TRUE,
  #   check = TRUE)
  # 
  # # get prediction error on test set
  # slbl.Yhat.test = predict.glm(
  #   slbl$glm,
  #   newdata = data.frame(V1 = balance::balance.fromSBP(
  #     x = XTe, y = slbl.coefs$sbp)),
  #   type = "response")
  # 
  # slbl.metrics = c(
  #   acc = mean((slbl.Yhat.test > 0.5) == Y2Te),
  #   # < 0.5 bc order of levels = c(case, control) instead of c(control, case)
  #   auc = pROC::roc(
  #     YTe, slbl.Yhat.test, levels = c("no", "CD"), direction = "<")$auc,
  #   percselected = sum(slbl.coefs$sbp > 0) / p,
  #   f1 = getF1(Y2Te, slbl.Yhat.test > 0.5),
  #   time = slbl.timing
  # )
  # 
  # saveRDS(
  #   slbl.metrics,
  #   paste0(output_dir, "/selbal_metrics", file.end))
  # 
  # slbl_sbp = slbl.coefs$sbp
  # if(slbl$glm$coefficients[2] < 0){
  #   slbl_sbp = -slbl_sbp
  # }
  # saveRDS(
  #   slbl_sbp,
  #   paste0(output_dir, "/selbal_sbp", file.end)
  # )
  # 
  # # codacore ###################################################################
  # start.time = Sys.time()
  # codacore0 = codacore::codacore(
  #   x = XTr, y = Y2Tr, logRatioType = "ILR",
  #   objective = "binary classification", cvParams = list(numFolds = K))
  # end.time = Sys.time()
  # codacore0.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get prediction error on test set and gamma-hat
  # if(length(codacore0$ensemble) > 0){ # at least 1 log-ratio found
  #   codacore0_SBP = matrix(0, nrow = p, ncol = length(codacore0$ensemble))
  #   codacore0_coeffs = rep(NA, length(codacore0$ensemble))
  #   for(col.idx in 1:ncol(codacore0_SBP)){
  #     codacore0_SBP[
  #       codacore0$ensemble[[col.idx]]$hard$numerator, col.idx] = 1
  #     codacore0_SBP[
  #       codacore0$ensemble[[col.idx]]$hard$denominator, col.idx] = -1
  #     codacore0_coeffs[col.idx] = codacore0$ensemble[[col.idx]]$slope
  #   }
  #   codacore0.betahat = getBetaFromCodacore(
  #     SBP_codacore = codacore0_SBP, coeffs_codacore = codacore0_coeffs, p = p)
  #   codacore0.Yhat.test = predict(codacore0, XTe) # before sigmoid
  #   # adjust codacore_SBP to correspond to positive theta-hats #################
  #   for(col in 1:ncol(codacore0_SBP)){
  #     if(codacore0_coeffs[col] < 0){
  #       codacore0_SBP[, col] = -codacore0_SBP[, col]
  #     }
  #   }
  # } else{
  #   print(paste0("sim ", i, " -- codacore has no log-ratios"))
  #   codacore0_coeffs = c()
  #   SBP_codacore = matrix(0, nrow = p, ncol = 1) ###############################
  #   codacore0model = stats::glm(Y2Tr ~ 1, family = "binomial")
  #   codacore0.betahat = rep(0, p)
  #   codacore0.Yhat.test = predict(codacore0model, XTe) # before sigmoid
  # }
  # rownames(codacore0_SBP) = colnames(XTr) ######################################
  # 
  # codacore0.metrics = c(
  #   acc = mean((codacore0.Yhat.test > 0) == Y2Te),
  #   auc = pROC::roc(
  #     Y2Te, codacore0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
  #   percselected = sum(abs(codacore0.betahat) > 10e-8) / p,
  #   f1 = getF1(Y2Te, codacore0.Yhat.test > 0),
  #   time = codacore0.timing
  # )
  # 
  # saveRDS(
  #   codacore0.metrics,
  #   paste0(output_dir, "/codacore_metrics", file.end))
  # 
  # saveRDS(
  #   codacore0_SBP, 
  #   paste0(output_dir, "/codacore_sbp", file.end)
  # )
  # 
  # # log-ratio lasso ############################################################
  # library(logratiolasso)
  # source("slr_analyses/Functions/logratiolasso.R")
  # WTr.c = scale(log(XTr), center = TRUE, scale = FALSE)
  # 
  # start.time = Sys.time()
  # lrl <- cv_two_stage(
  #   z = WTr.c, y = Y2Tr, n_folds = K, family="binomial")
  # end.time = Sys.time()
  # lrl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get prediction error on test set
  # WTe.c = scale(log(XTe), center = TRUE, scale = FALSE)
  # lrl.Yhat.test = as.numeric(WTe.c %*% lrl$beta_min) # before sigmoid
  # 
  # lrl.metrics = c(
  #   acc = mean((lrl.Yhat.test > 0) == Y2Te),
  #   auc = pROC::roc(
  #     Y2Te, lrl.Yhat.test, levels = c(0, 1), direction = "<")$auc,
  #   percselected = sum(abs(lrl$beta_min) > 10e-8) / p,
  #   f1 = getF1(Y2Te, lrl.Yhat.test > 0),
  #   time = lrl.timing
  # )
  # 
  # saveRDS(
  #   lrl.metrics,
  #   paste0(output_dir, "/lrlasso_metrics", file.end))
  
}
