# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 3/7/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs"

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
  library(propr)
  
  library(pROC)
  
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  source("Kristyn/Functions/supervisedlogratioseta.R")
  source("Kristyn/Functions/HSClust.R")
  source("Kristyn/Functions/slrnew.R")
  source("Kristyn/Functions/codalasso.R")
  
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  source("Kristyn/Functions/helper_functions.R")
  
  # for plots
  library(ggraph) # make dendrogram
  library(igraph) # transform dataframe to graph object: graph_from_data_frame()
  library(tidygraph)
  
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
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 1 # 1, 2
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
  is0.beta = !non0.beta
  bspars = sum(non0.beta)
  # solve for beta
  c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
  beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
    as.vector(ilrtrans.true$ilr.trans)

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
  classo.AUC.train = roc(Y, classo.Yhat.train)$auc
  # get prediction error on test set
  classo.Yhat.test = predict(classo, X.test)
  classo.AUC.test = roc(Y, classo.Yhat.test)$auc
  # beta estimation accuracy, selection accuracy #
  cl.metrics = getMetricsLLC(
    betahat = cl.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta, 
    true.beta = beta.true, 
    metrics = c("betaestimation", "selection"), classification = TRUE)
  cl.metrics = c(
    AUCtr = classo.AUC.train, AUCte = classo.AUC.test, cl.metrics)

  saveRDS(c(
    cl.metrics,
    "betaSparsity" = bspars,
     
    "time" = cl.timing
  ),
  paste0(output_dir, "/metrics", "/classo_metrics", file.end))

  ##############################################################################
  # new slr method (a balance regression method)
  #   -- hierarchical spectral clustering (with rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slrnew = slr(x = X, y = Y, rank1approx = TRUE, classification = TRUE)
  end.time = Sys.time()
  slrnew.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrnew_activevars = names(slrnew$index)
  slrnew_SBP = matrix(slrnew$index)
  rownames(slrnew_SBP) = slrnew_activevars
  
  slrnew.coefs = coefficients(slrnew$model)
  slrnew.a0 = slrnew.coefs[1]
  slrnew.thetahat = slrnew.coefs[-1]

  slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
  slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrnew.betahat) = colnames(X)
  slrnew.betahat[slrnew_activevars, ] =
    as.numeric(slrnew.betahat.nonzero)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slrnew.Yhat.train = predict.glm(
    slrnew$model, newdata = data.frame(X), type = "response")
  slrnew.AUC.train = roc(Y, slrnew.Yhat.train)$auc
  # get prediction error on test set
  slrnew.Yhat.test = predict.glm(
    slrnew$model, newdata = data.frame(X.test), type = "response")
  slrnew.AUC.test = roc(Y, slrnew.Yhat.test)$auc
  # beta estimation accuracy, selection accuracy #
  slrnew.metrics = getMetricsBalanceReg(
    thetahat = slrnew.thetahat, betahat = slrnew.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slrnew.metrics = c(
    AUCtr = slrnew.AUC.train, AUCte = slrnew.AUC.test, slrnew.metrics)

  saveRDS(c(
    slrnew.metrics,
    "betaSparsity" = bspars,

    "time" = slrnew.timing
  ),
  paste0(output_dir, "/metrics", "/slr_approx_metrics", file.end))

  ##############################################################################
  # new slr method (a balance regression method)
  #   -- hierarchical spectral clustering (without rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slrnew2 = slr(x = X, y = Y, rank1approx = FALSE, classification = TRUE)
  end.time = Sys.time()
  slrnew2.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrnew2_activevars = names(slrnew2$index)
  slrnew2_SBP = matrix(slrnew2$index)
  rownames(slrnew2_SBP) = slrnew2_activevars
  
  slrnew2.coefs = coefficients(slrnew2$model)
  slrnew2.a0 = slrnew2.coefs[1]
  slrnew2.thetahat = slrnew2.coefs[-1]

  slrnew2.betahat.nonzero = getBetaFromTheta(slrnew2.thetahat, sbp = slrnew2_SBP)
  slrnew2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrnew2.betahat) = colnames(X)
  slrnew2.betahat[slrnew2_activevars, ] =
    as.numeric(slrnew2.betahat.nonzero)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slrnew2.Yhat.train = predict.glm(
    slrnew2$model, newdata = data.frame(X), type = "response")
  slrnew2.AUC.train = roc(Y, slrnew2.Yhat.train)$auc
  # get prediction error on test set
  slrnew2.Yhat.test = predict.glm(
    slrnew2$model, newdata = data.frame(X.test), type = "response")
  slrnew2.AUC.test = roc(Y.test, slrnew2.Yhat.test)$auc
  # beta estimation accuracy, selection accuracy #
  slrnew2.metrics = getMetricsBalanceReg(
    thetahat = slrnew2.thetahat, betahat = slrnew2.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slrnew2.metrics = c(
    AUCtr = slrnew2.AUC.train, AUCte = slrnew2.AUC.test, slrnew2.metrics)

  saveRDS(c(
    slrnew2.metrics,
    "betaSparsity" = bspars,

    "time" = slrnew2.timing
  ),
  paste0(output_dir, "/metrics", "/slr_metrics", file.end))

  ##############################################################################
  # DBA method (a balance regression method)
  ##############################################################################
  Y.dba = factor(Y, levels = c(0, 1))
  
  start.time = Sys.time()
  dba_SBPfull = sbp.fromADBA(x = X, group = Y.dba)
  dba_SBP = sbp.subset(dba_SBPfull)
  dba = cvBMLasso(
    y = Y, X = X, sbp = dba_SBP, nlam = nlam, nfolds = K, 
    intercept = intercept, standardize = scaling, classification = TRUE)
  end.time = Sys.time()
  dba.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  dba.lam.min.idx = which.min(dba$cvm)
  dba.a0 = dba$theta0[dba.lam.min.idx]
  dba.thetahat = dba$theta[, dba.lam.min.idx]
  dba.betahat = getBetaFromTheta(dba.thetahat, sbp = dba$sbp)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  dba.Yhat.train = as.vector(predict(
    dba$cv.glmnet, newx = balance.fromSBP(X, dba_SBP), s = "lambda.min"))
  dba.AUC.train = roc(Y, dba.Yhat.train)$auc
  # get prediction error on test set
  dba.Yhat.test = as.vector(predict(
    dba$cv.glmnet, newx = balance.fromSBP(X.test, dba_SBP), s = "lambda.min"))
  dba.AUC.test = roc(Y.test, dba.Yhat.test)$auc
  # beta estimation accuracy, selection accuracy #
  dba.metrics = getMetricsBalanceReg(
    thetahat = dba.thetahat, betahat = dba.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  dba.metrics = c(
    AUCtr = dba.AUC.train, AUCte = dba.AUC.test, dba.metrics)
  
  # # plot the tree given by slr-hsc, indicating significant covariates
  # dba_leaf_types = rep("covariate", nrow(dba$sbp))
  # dba_balance_types = rep("balance", ncol(dba$sbp))
  # dba_nodes_types = data.frame(
  #   name = c(colnames(dba$sbp), rownames(dba$sbp)),
  #   type = c(dba_balance_types, dba_leaf_types)
  # )
  # plotSBP(dba$sbp, title = "DBA", nodes_types = dba_nodes_types)

  saveRDS(c(
    dba.metrics,
    "betaSparsity" = bspars,

    "time" = dba.timing
  ),
  paste0(output_dir, "/metrics", "/dba_metrics", file.end))

  ##############################################################################
  # selbal method (a balance regression method)
  ##############################################################################
  library(selbal) # masks stats::cor()
  
  X.slbl = X
  Y.slbl = factor(Y, levels = c(0, 1))
  rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
  colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
  
  start.time = Sys.time()
  slbl = selbal.cv(x = X.slbl, y = Y.slbl, n.fold = K)
  end.time = Sys.time()
  slbl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  # U (transformation) matrix
  U.slbl = rep(0, p)
  names(U.slbl) = colnames(X.slbl)
  pba.pos = unlist(subset(
    slbl$global.balance, subset = Group == "NUM", select = Taxa))
  num.pos = length(pba.pos)
  pba.neg = unlist(subset(
    slbl$global.balance, subset = Group == "DEN", select = Taxa))
  num.neg = length(pba.neg)
  U.slbl[pba.pos] = 1 / num.pos
  U.slbl[pba.neg] = -1 / num.neg
  norm.const = sqrt((num.pos * num.neg) / (num.pos + num.neg))
  U.slbl = norm.const * U.slbl
  # check: these are equal
  # glm(Y.slbl ~ log(X.slbl) %*% as.matrix(U.slbl), family = binomial(link = "logit"))
  # slbl$glm
  slbl.thetahat = coefficients(slbl$glm)[-1]
  slbl.betahat = U.slbl %*% as.matrix(slbl.thetahat)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slbl.Yhat.train = predict.glm(
    slbl$glm, newdata = data.frame(X.slbl), type = "response")
  # roc_obj = roc(Y, slbl.Yhat.train)
  # plot(roc_obj)
  # roc_obj$auc
  slbl.AUC.train = roc(Y, slbl.Yhat.train)$auc
  # get prediction error on test set
  X.slbl.test = X.test
  rownames(X.slbl.test) = paste("Sample", 1:nrow(X.slbl.test), sep = "_")
  colnames(X.slbl.test) = paste("V", 1:ncol(X.slbl.test), sep = "")
  slbl.Yhat.test = predict.glm(
    slbl$glm, newdata = data.frame(X.slbl.test), type = "response")
  slbl.AUC.test = roc(Y.test, slbl.Yhat.test)$auc
  # beta estimation accuracy, selection accuracy #
  slbl.metrics = getMetricsBalanceReg(
    thetahat = slbl.thetahat, betahat = slbl.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slbl.metrics = c(AUCtr = slbl.AUC.train, AUCte = slbl.AUC.test, slbl.metrics)
  
  saveRDS(c(
    slbl.metrics,
    "betaSparsity" = bspars,

    "time" = slbl.timing
  ),
  paste0(output_dir, "/metrics", "/selbal_metrics", file.end))
  
  ##############################################################################
  # CoDaCoRe (a balance regression method)
  ##############################################################################
  library(codacore)
  
  if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
    reticulate::use_condaenv("anaconda3")
  }
  # tensorflow::install_tensorflow()
  # keras::install_keras()
  
  start.time = Sys.time()
  codacore0 = codacore(
    x = X, y = Y, logRatioType = "ILR", # instead of "balance" ?
    objective = "binary classification") 
  end.time = Sys.time()
  codacore0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  if(length(codacore0$ensemble) > 0){
    
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
    codacore0.AUC.train = roc(Y, codacore0.Yhat.train)$auc
    # get prediction error on test set
    codacore0.Yhat.test = predict(codacore0, X.test)
    codacore0.AUC.test = roc(Y.test, codacore0.Yhat.test)$auc
    # beta estimation accuracy, selection accuracy #
    codacore0.metrics = getMetricsBalanceReg(
      betahat = codacore0.betahat,
      true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
      true.beta = beta.true, metrics = c("betaestimation", "selection"))
    codacore0.metrics = c(
      AUCtr = codacore0.AUC.train, AUCte = codacore0.AUC.test, codacore0.metrics)
  } else{
    codacore0.metrics = c(
      AUCtr = NA, AUCte = NA, 
      EA1 = NA, EA2 = NA, EAInfty = NA, 
      EA1Active = NA, EA2Active = NA, EAInftyActive = NA, 
      EA1Inactive = NA, EA2Inactive = NA, EAInftyInactive = NA, 
      FP = 0, FN = p, TPR = 0, precision = NA, Fscore = NA, 
      "FP+" = 0, "FN+" = sum(SBP.true != 0), "TPR+" = 0, 
      "FP-" = p - sum(SBP.true != 0), "FN-" = 0, "TPR-" = 0
    )
  }
  
  saveRDS(c(
    codacore0.metrics,
    "betaSparsity" = bspars,
    
    "time" = codacore0.timing
  ),
  paste0(output_dir, "/metrics", "/codacore_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


