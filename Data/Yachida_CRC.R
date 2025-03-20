# Purpose: compare methods on data sets
# Date: 03/20/2025

jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

output_dir = "../Outputs/data_metrics/"

library(mvtnorm)
library(tidyverse)

library(Matrix)
library(glmnet)

library(selbal)
library(compositions)

source("Functions/slrs.R")
source("Functions/util.R")

# tuning parameter settings ----
run_CODALASSO <- T
run_SELBAL <- T
run_CODACORE <- T
run_slr_PC <- T

K = 5
alpha <- 1
ncomp <- 1
family = 'binomial'
hparam = "1se"
split.perc = 0.7 # partition into train and test data set
covar.train <- NULL
covar.test <- NULL

type.measure <- if (family=='gaussian'){
  "mse"
} else if (family=='binomial'){
  "auc"
} else if (family=='cox'){
  "deviance"
}
filename <- 'CRC'
load("../Data/Yachida_CRC_2019.rda")

W = data.list$X # n by p
p <- ncol(W)
n <- nrow(W)
Y2 = as.factor(ifelse(data.list$y == "Control", 0, 1))

###0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])####
W = as.matrix(selbal::cmultRepl2(W, zero.rep = "bayes"))

file.end = paste0(
  "/",filename,
  "_split", split.perc, 
  "_hparam", hparam, 
  "_",jid)

if (file.exists(paste0("../Outputs/datasets_new/",file.end,".rds"))){
  split.data <- readRDS(paste0("../Outputs/datasets_new/",file.end,".rds"))
  X <- split.data$X
  X.test <- split.data$X.test
  Y <- split.data$Y
  Y.test <- split.data$Y.test
} else{
  train.idx <- sample(1:n,split.perc*n) 
  
  X <- W[train.idx,]
  X.test <- W[-train.idx,]
  Y <- Y2[train.idx]
  Y.test <- Y2[-train.idx]
  split.data <- list(X=X,X.test = X.test,Y=Y,Y.test=Y.test)
  
  saveRDS(
    split.data,
    file=paste0(
      "../Outputs/datasets_new/",file.end,
      ".rds"))  
}



# codalasso #######################################################################
if (run_CODALASSO){
  
  source("Functions/codalasso.R")
  start.time = Sys.time()
  if(hparam == "min"){
    classo = codalasso(
      X, Y, numFolds = K, gamma = 0, type.measure = type.measure, stratify = FALSE)
  } else if(hparam == "1se"){
    classo = codalasso(
      X, Y, numFolds = K, gamma = 1, type.measure = type.measure, stratify = FALSE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  cl1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  saveRDS(
    classo,
    paste0(
      output_dir, file.end,
      "_classo",
      ".rds"))
  
  # get prediction error on test set
  classo1.Yhat.test = predict(classo, X.test) # before sigmoid
  
  cl1.metrics = c(
    acc = mean((classo1.Yhat.test > 0) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, classo1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum(abs(classo$cll$betas[-1]) > 10e-8) ,
    f1 = getF1(Y.test, classo1.Yhat.test > 0),
    time = cl1.timing
  )
  
  saveRDS(
    cl1.metrics,
    paste0(
      output_dir, file.end,
      "_classo_metrics",
      ".rds"))
  cat('\n Test set performance by classo ...', cl1.metrics, '\n')
  
}

# slr ----
if (run_slr_PC){
  cat('slr with constrained PC ...','\n')
  start.time = Sys.time()
  slrPCcv = cv.slr(x = X, y = Y, covar=NULL, method = "C", 
                   family=family, zeta = 0, ncomp = ncomp,
                   nfolds = K, type.measure = type.measure,
                   trace.it = FALSE)
  if(hparam == "min"){
    slrPC = slr(x = X, y = Y, covar=NULL, method = "C", 
                family=family, zeta = 0, ncomp = slrPCcv$nc["min",],
                threshold = slrPCcv$threshold[slrPCcv$index["min",]])
  } else if(hparam == "1se"){
    slrPC = slr(x = X, y = Y, covar=NULL, method = "C", 
                family=family, zeta = 0, ncomp = slrPCcv$nc["1se",],
                threshold = slrPCcv$threshold[slrPCcv$index["1se",]])
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slrPC.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  saveRDS(
    slrPC,
    paste0(
      output_dir, file.end,
      "_slr_constrainedPC",
      ".rds"))
  
  # compute metrics on the selected model #
  # get prediction error on test set
  slrPC.Yhat.test = predict.slr(slrPC,
                                X.test,
                                family = family)
  
  slrPC.metrics = c(
    acc = mean((slrPC.Yhat.test > 0.5) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, slrPC.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum( abs(slrPC$`log-contrast coefficients`) > 0) ,
    f1 = getF1(Y.test, slrPC.Yhat.test > 0.5),
    time = slrPC.timing
  )
  
  cat('\n Test set performance ...', slrPC.metrics, '\n')
  
  saveRDS(
    slrPC.metrics,
    paste0(
      output_dir, file.end,
      "_slr_constrainedPC_metrics",
      ".rds"))
}

# codacore - 1 balance #########################################################
if (run_CODACORE){
  cat('codacore1 ...','\n')
  
  library(codacore)
  start.time = Sys.time()
  if(hparam == "min"){
    codacore1 = codacore::codacore(
      x = X, y = Y, logRatioType = "ILR",
      objective = "binary classification", cvParams = list(numFolds = K),
      maxBaseLearners = 1,
      lambda = 0, fast = FALSE)
  } else if(hparam == "1se"){
    codacore1 = codacore::codacore(
      x = X, y = Y, logRatioType = "ILR",
      objective = "binary classification", cvParams = list(numFolds = K),
      maxBaseLearners = 1,
      lambda = 1, fast = FALSE)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  codacore1.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  saveRDS(
    codacore1,
    paste0(
      output_dir, file.end,
      "_codacore1",
      ".rds"))
  
  length(codacore1$ensemble) # one balance selected
  # numerator (I+) / denominator (I-) of selected balance
  codacore1$ensemble[[1]]$slope # positive (if negative, num -> den & vice versa)
  colnames(X)[codacore1$ensemble[[1]]$hard$numerator]
  colnames(X)[codacore1$ensemble[[1]]$hard$denominator]
  sum(c(codacore1$ensemble[[1]]$hard$numerator, codacore1$ensemble[[1]]$hard$denominator))
  
  # for overleaf
  cat('\n codacore1: Variables in numerator:\n')
  cat(paste(str_replace_all(
    colnames(X)[codacore1$ensemble[[1]]$hard$numerator],
    fixed("_"), "\\_"), collapse = ", "))
  cat('\n codacore1: Variables in denominator:\n')
  cat(paste(str_replace_all(
    colnames(X)[codacore1$ensemble[[1]]$hard$denominator],
    fixed("_"), "\\_"), collapse = ", "))
  
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
    codacore1.Yhat.test = predict(codacore1, X.test)
    # adjust codacore_SBP to correspond to positive theta-hats
    if(codacore1_coeffs[1] < 0){
      codacore1_SBP[, 1] = -codacore1_SBP[, 1]
    }
  } else{
    print(paste0("sim ", b, " -- codacore has no log-ratios"))
    codacore1_coeffs = c()
    codacore1_SBP = matrix(0, nrow = p, ncol = 1)
    codacore1model = stats::glm(Y ~ 1, family = family)
    codacore1.betahat = rep(0, p)
    codacore1.Yhat.test = predict(codacore1model, X.test)
  }
  rownames(codacore1_SBP) = colnames(X)
  
  codacore1.metrics = c(
    acc = mean((codacore1.Yhat.test > 0.5) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, codacore1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum(abs(codacore1.betahat) > 10e-8) ,
    f1 = getF1(Y.test, codacore1.Yhat.test > 0.5),
    time = codacore1.timing
  )
  
  cat('\n Test set performance by codacore1 ...', codacore1.metrics, '\n')
  
  saveRDS(
    codacore1.metrics,
    paste0(
      output_dir, file.end,
      "_codacore1_metrics",
      ".rds"))
}


# selbal #######################################################################
if (run_SELBAL){
  cat('selbal ...','\n')
  
  start.time = Sys.time()
  if(hparam == "min"){
    slbl = selbal::selbal.cv(x = X, y = Y, n.fold = K, opt.cri = "min",
                             seed = jid, maxV = 100)
  } else if(hparam == "1se"){
    slbl = selbal::selbal.cv(x = X, y = Y, n.fold = K, opt.cri = "1se",
                             seed = jid, maxV = 50)
  } else{
    stop("invalid hparam setting (method for selecting hyperparameter(s)).")
  }
  end.time = Sys.time()
  slbl0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  saveRDS(
    slbl,
    paste0(
      output_dir, file.end,
      "_selbal",
      ".rds"))
  
  # get theta-hat and gamma-hat
  slbl.sbp = getSBPSelbalcv(X,slbl)
  
  # get prediction error on test set
  slbl.Yhat.test = predict.glm(
    slbl$glm,
    newdata = data.frame(V1 = balance::balance.fromSBP(X.test, slbl.sbp)),
    type = "response")
  
  slbl.metrics = c(
    acc = mean((slbl.Yhat.test > 0.5) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, slbl.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum(slbl.sbp != 0),
    f1 = getF1(Y.test, slbl.Yhat.test > 0.5),
    time = slbl0.timing
  )
  cat('\n Test set performance by selbal ...', slbl.metrics, '\n')
  
  saveRDS(
    slbl.metrics,
    paste0(
      output_dir, file.end,
      "_selbal_metrics",
      ".rds"))
}
