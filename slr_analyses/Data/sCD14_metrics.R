# Purpose: compare slr to other methods on data sets
# Date: 6/29/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs_metrics"

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
  
  library(pROC)
  
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
  hparam = "1se"
  K = 10
  nlam = 100
  intercept = TRUE
  scaling = TRUE
  tol = 1e-4
  
  filter.perc = 0.8
  split.perc = 0.7
  
  file.end = paste0(
    "_sCD14", 
    "_hparam", hparam,
    "_gbm",
    "_sim", b,
    ".rds")
  
  ################################################################################
  # sCD14: another HIV data set in selbal package
  #   n = 151 samples (a subset from sCD14 data set), 
  #   p = 60 taxa (counts for microbial taxa at genus level), 
  #   1 response (sCD14 - continuous)
  W = selbal::sCD14[, 1:60]
  W.origin <- W
  W <- W[,apply(W==0,2,mean)<filter.perc]
  X = sweep(W, 1, rowSums(W), FUN='/')
  Y = selbal::sCD14[, 61]
  p = ncol(W)
  
  ##############################################################################
  # 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
  X_gbm = selbal::cmultRepl2(W, zero.rep = "bayes")
  
  ##############################################################################
  # Train/Test Split
  #   Following Gordon-Rodriguez et al. 2022, fit each method on 20 random 80/20
  #     train/test splits, 
  #     -- since response is continuous, no need to stratify by case-control.
  numObs = nrow(X_gbm)
  inputDim = ncol(X_gbm)
  if(file.exists(paste0(
    output_dir, "/data", "_split", split.perc, "filter", filter.perc, 
    file.end))){
    data.tmp = readRDS(paste0(
      output_dir, "/data", "_split", split.perc, "filter", filter.perc, 
      file.end))
    XTr = data.tmp$XTr
    XTe = data.tmp$XTe
    YTr = data.tmp$YTr
    YTe = data.tmp$YTe
  } else{
    # trainIdx = sample(cut(1:numObs, breaks=5, labels=F))
    trainIdx = 1-seq(numObs) %in% sample(1:numObs,size=ceiling(split.perc*numObs))
    XTr = X_gbm[trainIdx != 1,]
    YTr = Y[trainIdx != 1]
    XTe = X_gbm[trainIdx == 1,]
    YTe = Y[trainIdx == 1]
    
    saveRDS(list(
      XTr = XTr, YTr = YTr,
      XTe = XTe, YTe = YTe
    ),
    paste0(output_dir, "/data", file.end))
  }
  
  ##############################################################################
  # fit methods
  ##############################################################################
  
  # classo #####################################################################
  start.time = Sys.time()
  classo = cv.func(
    method="ConstrLasso", y = YTr, x = log(XTr), Cmat = matrix(1, p, 1),
    nlam = nlam, nfolds = K, tol = tol, intercept = intercept,
    scaling = scaling)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # get gamma-hat
  if(hparam == "min"){
    cl.lam.idx = which.min(classo$cvm)
  } else if(hparam == "1se"){
    oneSErule = min(classo$cvm) + classo$cvsd[which.min(classo$cvm)] * 1
    cl.lam.idx = which(classo$cvm <= oneSErule)[1]
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  cl.a0 = classo$int[cl.lam.idx]
  cl.betahat = classo$bet[, cl.lam.idx]
  
  # get prediction error on test set
  classo.Yhat.test = cl.a0 + log(as.matrix(XTe)) %*% cl.betahat
  
  cl.metrics = c(
    mse = as.vector(crossprod(YTe - classo.Yhat.test) / nrow(XTe)),
    percselected = sum(abs(cl.betahat) > 10e-8) / p,
    time = cl.timing
  )
  
  saveRDS(
    cl.metrics,
    paste0(output_dir, "/classo_metrics", file.end))
  
  # # slr - spectral #############################################################
  # start.time = Sys.time()
  # slrspeccv = cv.slr(
  #   x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
  #   response.type = "continuous", s0.perc = 0, zeta = 0,
  #   nfolds = K, type.measure = "mse",
  #   scale = scaling, trace.it = FALSE)
  # if(hparam == "min"){
  #   slrspec = slr(
  #     x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
  #     response.type = "continuous", s0.perc = 0, zeta = 0,
  #     threshold = slrspeccv$threshold[slrspeccv$index["min",]],
  #     positive.slope = TRUE)
  # } else if(hparam == "1se"){
  #   slrspec = slr(
  #     x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
  #     response.type = "continuous", s0.perc = 0, zeta = 0,
  #     threshold = slrspeccv$threshold[slrspeccv$index["1se",]],
  #     positive.slope = TRUE)
  # } else{
  #   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  # }
  # end.time = Sys.time()
  # slrspec.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get SBP
  # slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
  # rownames(slrspec.fullSBP) = colnames(X)
  # slrspec.fullSBP[match(
  #   names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp
  # 
  # # get prediction error on test set
  # slrspec.Yhat.test = predict(
  #   slrspec$fit,
  #   data.frame(balance = slr.fromContrast(XTe, slrspec.fullSBP)),
  #   type = "response")
  # 
  # slrspec.metrics = c(
  #   mse = as.vector(crossprod(YTe - slrspec.Yhat.test) / nrow(XTe)),
  #   percselected = sum(slrspec.fullSBP > 0) / p,
  #   time = slrspec.timing
  # )
  # 
  # saveRDS(
  #   slrspec.metrics,
  #   paste0(output_dir, "/slr_spectral_metrics", file.end))
  # 
  # if(!all(slrspec.fullSBP == 0) & slrspec$theta[2] < 0){
  #   slrspec.fullSBP = -slrspec.fullSBP
  # }
  # saveRDS(
  #   slrspec.fullSBP,
  #   paste0(output_dir, "/slr_spectral_sbp", file.end)
  # )
  # 
  # # slr - hierarchical #########################################################
  # start.time = Sys.time()
  # slrhiercv = cv.slr(
  #   x = XTr, y = YTr, screen.method = "wald", cluster.method = "hierarchical",
  #   response.type = "continuous", s0.perc = 0, zeta = 0,
  #   nfolds = K, type.measure = "mse",
  #   scale = scaling, trace.it = FALSE)
  # if(hparam == "min"){
  #   slrhier = slr(
  #     x = XTr, y = YTr, screen.method = "wald", cluster.method = "hierarchical",
  #     response.type = "continuous", s0.perc = 0, zeta = 0,
  #     threshold = slrhiercv$threshold[slrhiercv$index["min",]],
  #     positive.slope = TRUE)
  # } else if(hparam == "1se"){
  #   slrhier = slr(
  #     x = XTr, y = YTr, screen.method = "wald", cluster.method = "hierarchical",
  #     response.type = "continuous", s0.perc = 0, zeta = 0,
  #     threshold = slrhiercv$threshold[slrhiercv$index["1se",]],
  #     positive.slope = TRUE)
  # } else{
  #   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  # }
  # end.time = Sys.time()
  # slrhier.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get SBP
  # slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
  # rownames(slrhier.fullSBP) = colnames(X)
  # slrhier.fullSBP[match(
  #   names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
  # 
  # # get prediction error on test set
  # slrhier.Yhat.test = predict(
  #   slrhier$fit,
  #   data.frame(balance = slr.fromContrast(XTe, slrhier.fullSBP)),
  #   type = "response")
  # 
  # slrhier.metrics = c(
  #   mse = as.vector(crossprod(YTe - slrhier.Yhat.test) / nrow(XTe)),
  #   percselected = sum(slrhier.fullSBP > 0) / p,
  #   time = slrhier.timing
  # )
  # 
  # saveRDS(
  #   slrhier.metrics,
  #   paste0(output_dir, "/slr_hierarchical_metrics", file.end))
  # 
  # if(!all(slrhier.fullSBP == 0) & slrhier$theta[2] < 0){
  #   slrhier.fullSBP = -slrhier.fullSBP
  # }
  # saveRDS(
  #   slrhier.fullSBP,
  #   paste0(output_dir, "/slr_hierarchical_sbp", file.end)
  # )
  # 
  # # selbal #####################################################################
  # start.time = Sys.time()
  # if(hparam == "min"){
  #   slbl = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K, opt.cri = "max")
  # } else if(hparam == "1se"){
  #   slbl = selbal::selbal.cv(x = XTr, y = YTr, n.fold = K, opt.cri = "1se")
  # } else{
  #   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  # }
  # end.time = Sys.time()
  # slbl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get theta-hat and gamma-hat
  # slbl.coefs = getCoefsSelbal(
  #   X = XTr, y = YTr, selbal.fit = slbl, classification = FALSE,
  #   check = TRUE)
  # 
  # # if log-ratio pair is selected, test its significance (F test of null model)
  # not.signif = FALSE
  # if(sum(slbl.coefs$sbp == 1) == 1 & sum(slbl.coefs$sbp == -1) == 1){
  #   ilrU = getIlrTrans(sbp = slbl.coefs$sbp)
  #   data.slbl.tr = data.frame(cbind(YTr, log(as.matrix(XTr)) %*% matrix(ilrU)))
  #   colnames(data.slbl.tr)[ncol(data.slbl.tr)] <- "V1"
  #   supposed.fit = lm(YTr ~ ., data = data.slbl.tr)
  #   supposed.fit.summary = summary(supposed.fit)
  #   pval.tmp = 1 - pf(supposed.fit.summary$fstatistic["value"], supposed.fit.summary$fstatistic["numdf"], supposed.fit.summary$fstatistic["dendf"])
  #   null.fit = lm(YTr ~ 1)
  #   Ftest = anova(null.fit, supposed.fit, test = "F")
  #   if(Ftest$`Pr(>F)`[2] > 0.05){ # if it's not significant, select null model
  #     not.signif = TRUE
  #   }
  # }
  # 
  # if(not.signif){
  #   
  #   # get prediction error on test set
  #   slbl.Yhat.test = predict.glm(
  #     null.fit,
  #     newdata = data.frame(V1 = balance::balance.fromSBP(
  #       x = XTe, y = slbl.coefs$sbp)),
  #     type = "response")
  #   
  #   slbl.metrics = c(
  #     mse = as.vector(crossprod(YTe - slbl.Yhat.test) / nrow(XTe)),
  #     percselected = sum(slbl.coefs$sbp > 0) / p,
  #     time = slbl.timing,
  #     percselected.testlrpair = 0
  #   )
  # } else{
  #   
  #   # get prediction error on test set
  #   slbl.Yhat.test = predict.glm(
  #     slbl$glm,
  #     newdata = data.frame(V1 = balance::balance.fromSBP(
  #       x = XTe, y = slbl.coefs$sbp)),
  #     type = "response")
  #   
  #   slbl.metrics = c(
  #     mse = as.vector(crossprod(YTe - slbl.Yhat.test) / nrow(XTe)),
  #     percselected = sum(slbl.coefs$sbp > 0) / p,
  #     time = slbl.timing,
  #     percselected.testlrpair = sum(slbl.coefs$sbp > 0) / p
  #   )
  # }
  # 
  # saveRDS(
  #   slbl.metrics,
  #   paste0(output_dir, "/selbal_metrics", file.end))
  # 
  # slbl_sbp = slbl.coefs$sbp
  # if((length(slbl$glm$coefficients) > 1) & slbl$glm$coefficients[2] < 0){
  #   slbl_sbp = -slbl_sbp
  # }
  # saveRDS(
  #   slbl_sbp,
  #   paste0(output_dir, "/selbal_sbp", file.end)
  # )
  # 
  # # codacore ###################################################################
  # library(codacore)
  # if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
  #   reticulate::use_condaenv("anaconda3")
  # }
  # 
  # start.time = Sys.time()
  # if(hparam == "min"){
  #   codacore0 = codacore(
  #     x = XTr, y = YTr, logRatioType = "ILR",
  #     objective = "regression", cvParams = list(numFolds = K), 
  #     lambda = 0) 
  # } else if(hparam == "1se"){
  #   codacore0 = codacore(
  #     x = XTr, y = YTr, logRatioType = "ILR",
  #     objective = "regression", cvParams = list(numFolds = K), 
  #     lambda = 1) 
  # } else{
  #   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  # }
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
  #   
  #   names(codacore0_coeffs) = paste(
  #     "balance", 1:length(codacore0_coeffs), sep = "")
  #   rownames(codacore0_SBP) = colnames(X)
  #   codacore0.coefs2 = getCoefsBM(
  #     coefs = codacore0_coeffs * codacore0$yScale,
  #     sbp = codacore0_SBP)
  #   codacore0.betahat = codacore0.coefs2$llc.coefs
  #   
  #   codacore0.Yhat.test = predict(codacore0, XTe)
  #   # adjust codacore_SBP to correspond to positive theta-hats #################
  #   for(col in 1:ncol(codacore0_SBP)){
  #     if(codacore0_coeffs[col] < 0){
  #       codacore0_SBP[, col] = -codacore0_SBP[, col]
  #     }
  #   }
  # } else{
  #   print(paste0("sim ", b, " -- codacore has no log-ratios"))
  #   codacore0_coeffs = c()
  #   codacore0_SBP = matrix(0, nrow = p, ncol = 1) ###############################
  #   codacore0model = stats::glm(YTr ~ 1, family = "gaussian")
  #   codacore0.betahat = rep(0, p)
  #   codacore0.Yhat.test = predict(codacore0model, XTe)
  # }
  # rownames(codacore0_SBP) = colnames(XTr) ######################################
  # 
  # codacore0.metrics = c(
  #   mse = as.vector(crossprod(YTe - codacore0.Yhat.test) / nrow(XTe)),
  #   percselected = sum(abs(codacore0.betahat) > 10e-8) / p,
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
  # YTr.c = YTr - mean(YTr)
  # 
  # start.time = Sys.time()
  # if(hparam == "min"){
  #   lrl <- cv_two_stage(z = WTr.c, y = YTr.c, n_folds = K, gamma = 0)
  #   lrl.betahat = lrl$beta_min
  # } else if(hparam == "1se"){
  #   lrl <- cv_two_stage(z = WTr.c, y = YTr.c, n_folds = K, gamma = 1) 
  #   lrl.betahat = lrl$beta_gammase
  # } else{
  #   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  # }
  # end.time = Sys.time()
  # lrl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # get prediction error on test set
  # WTe.c = scale(log(XTe), center = TRUE, scale = FALSE)
  # lrl.Yhat.test = as.numeric(WTe.c %*% lrl.betahat)
  # 
  # lrl.metrics = c(
  #   mse = as.vector(crossprod(YTe - lrl.Yhat.test) / nrow(XTe)),
  #   percselected = sum(abs(lrl.betahat) > 10e-8) / p,
  #   time = lrl.timing
  # )
  # 
  # saveRDS(
  #   lrl.metrics,
  #   paste0(output_dir, "/lrlasso_metrics", file.end))
  
}
