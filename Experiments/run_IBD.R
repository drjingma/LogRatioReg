# update: 3/15/25

jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

# library(mvtnorm)
library(tidyverse)

library(Matrix)
library(glmnet)

library(selbal)
library(compositions)

source("Functions/slrs.R")
source("Functions/util.R")

# tuning parameter settings ----
K = 5
alpha <- 1
ncomp <- 1
family = 'binomial'
hparam = "1se"
type.measure <- if (family=='gaussian'){
  "mse"
} else if (family=='binomial'){
  "auc"
} else if (family=='cox'){
  "deviance"
}

run_SELBAL <- T
run_codaLasso <- T
run_CODACORE <- T
run_slr_PC <- T

output_dir = "../Outputs/data_final/"
filename <- 'IBD_highCor_sparsity1_ratio1_dim50'
model.params <- readRDS(paste0("../Data/",filename,".rds"))
X <- model.params$X
X.test <- model.params$X.test
sbp <- model.params$sbp
active.features <- names(model.params$sbp[model.params$sbp!=0,1])

datafile.name <- paste0("../Outputs/datasets_simulation","/",filename,
                        "_", jid,".rds")

if (file.exists(datafile.name)){
  data.list <- readRDS(datafile.name)
  Y <- data.list$Y
  Y.test <- data.list$Y.test
} else {
  # balance in train and test data
  z1 <- balance::balance.fromSBP(X,sbp)
  z1.test <- balance::balance.fromSBP(X.test,sbp)
  
  # generate response variable
  b1 <- 1 # coefficient for balance
  b0 <- log(53/34) - b1 * mean(z1)
  prob <- 1/(1+exp(-(b0 + b1 * z1)))
  prob.test <- 1/(1+exp(-(b0 + b1 * z1.test)))
  
  Y <- rbinom(n = nrow(X), size = 1, p = as.vector(prob))
  cat('High Cor ... training data', table(Y), "...\n")
  Y.test <- rbinom(n = nrow(X.test), size = 1, p = as.vector(prob.test))
  cat('High Cor ... test data', table(Y.test), "...\n")
  
  data.list <- list(Y = as.factor(Y), Y.test = as.factor(Y.test), b0 = b0, b1 = b1)
  
  saveRDS(data.list, file = datafile.name)
}
covar.test <- NULL
covar.train <- NULL

p = ncol(X)

# file.end for saving results
file.end <- paste0(
  "/",filename,
  "_hparam", hparam,
  "_", jid)

# codacore - 1 balance #########################################################
if (run_CODACORE){
  cat('codacore1 ...','\n')
  
  library(codacore)
  library(tensorflow)
  library(reticulate)
  start.time = Sys.time()
  if(hparam == "min"){
    codacore1 = codacore::codacore(
      x = X, y = Y, logRatioType = "balances",
      objective = "binary classification", cvParams = list(numFolds = K),
      maxBaseLearners = ncomp,
      lambda = 0, fast = FALSE)
  } else if(hparam == "1se"){
    codacore1 = codacore::codacore(
      x = X, y = Y, logRatioType = "balances",
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
  
  error <- getErrorRate(rownames(codacore1_SBP)[codacore1_SBP!=0],active.features,colnames(X))
  
  codacore1.metrics = c(
    acc = mean((codacore1.Yhat.test > 0.5) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, codacore1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum(abs(codacore1.betahat) > 10e-8) ,
    sensitivity = error$sensitivity,
    specificity = error$specificity,
    f1 = getF1(Y.test, codacore1.Yhat.test > 0.5),
    time = codacore1.timing
  )
  
  cat('\n Test set performance ...', codacore1.metrics, '\n')
  
  saveRDS(
    codacore1.metrics,
    paste0(
      output_dir, file.end,
      "_codacore1_metrics",
      ".rds"))
}


# codaLasso #######################################################################
if (run_codaLasso){
  cat('codaLasso ...','\n')
  
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
  
  error <- getErrorRate(classo$cll$`taxa with non-zero coeff`[-1],active.features,colnames(X))
  
  cl1.metrics = c(
    acc = mean((classo1.Yhat.test > 0) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, classo1.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum(abs(classo$cll$betas[-1]) > 10e-8) ,
    sensitivity = error$sensitivity,
    specificity = error$specificity,
    f1 = getF1(Y.test, classo1.Yhat.test > 0),
    time = cl1.timing
  )
  
  
  saveRDS(
    cl1.metrics,
    paste0(
      output_dir, file.end,
      "_classo_metrics",
      ".rds"))
  cat('\n Test set performance ...', cl1.metrics, '\n')
}

# slr - constrained PC  ----
if (run_slr_PC){
  source("lib/PLS-PBs/R: functions/fBalChip.r")
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
  
  error <- getErrorRate(names(slrPC$sbp[slrPC$sbp!=0,1]),active.features,colnames(X))
  
  slrPC.metrics = c(
    acc = mean((slrPC.Yhat.test > 0.5) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, slrPC.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = sum( abs(slrPC$`log-contrast coefficients`) > 0) ,
    sensitivity = error$sensitivity,
    specificity = error$specificity,
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
  slbl <- selbal::selbal(x = X, y = Y, th.imp = 0, 
                         maxV = slbl$opt.nvar, draw = F)
  end.time = Sys.time()
  slbl0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  saveRDS(
    slbl,
    paste0(
      output_dir, file.end,
      "_selbal",
      ".rds"))
  
  error <- getErrorRate(slbl$balance,active.features,colnames(X))
  
  slbl.sbp <- getSBPSelbal(X,slbl)
  
  # get prediction error on test set
  slbl.Yhat.test = predict.glm(
    slbl$glm,
    newdata = data.frame(V1 = balance::balance.fromSBP(X.test, slbl.sbp)),
    type = "response")
  
  slbl.metrics = c(
    acc = mean((slbl.Yhat.test > 0.5) == as.numeric(Y.test)-1),
    auc = pROC::roc(
      Y.test, slbl.Yhat.test, levels = c(0, 1), direction = "<")$auc,
    numselected = length(slbl$balance),
    sensitivity = error$sensitivity,
    specificity = error$specificity,
    f1 = getF1(Y.test, slbl.Yhat.test > 0.5),
    time = slbl0.timing
  )
  cat('\n Test set performance ...', slbl.metrics, '\n')
  
  saveRDS(
    slbl.metrics,
    paste0(
      output_dir, file.end,
      "_selbal_metrics",
      ".rds"))
}

