rm(list=ls())
# Purpose: observe different R^2 values for specified rho and sigma_eps values
#   for the one-balance model
#   this time, more generalizeable to different balances
# Date: 3/7/2022

################################################################################
# libraries and settings

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
source("Kristyn/Functions/simulatedata.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

# expdecay Sigma #############################################################

# Settings to toggle with
sigma.settings = "latentVarModel_binary" # !!!
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_eps2 = 0.01
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0
b1 = 1 # 1, 0.5, 0.25
a0 = 0 # 0
theta.value = 2 # weight on a1: 1

##############################################################################
# generate data
start.seed = 123
for(i in 1:100){
  print(paste0("sim ", i))
  seed = start.seed + i
  set.seed(seed)
  # get latent variable
  U.all = matrix(runif(2 * n), ncol = 1)
  # simulate y from latent variable
  sigmoid = function(x) 1 / (1 + exp(-x))
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
  # solve for beta
  c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
  beta.true = (b1 / c1plusc2) * theta.value * as.vector(ilrtrans.true$ilr.trans)
  
  ##############################################################################
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  classo = codalasso(X, Y, numFolds = K)
  
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
  
  ##############################################################################
  # CoDaCoRe (a balance regression method)
  ##############################################################################
  library(codacore)
  if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
    reticulate::use_condaenv("anaconda3")
  }
  
  codacore0 = codacore(
    x = X, y = Y, logRatioType = "ILR", # instead of "balance" ?
    objective = "binary classification") 
  if(length(codacore0$ensemble) > 0){
    codacore0_SBP = matrix(0, nrow = p, ncol = length(codacore0$ensemble))
    for(col.idx in 1:ncol(codacore0_SBP)){
      codacore0_SBP[codacore0$ensemble[[col.idx]]$hard$numerator, col.idx] = 1
      codacore0_SBP[codacore0$ensemble[[col.idx]]$hard$denominator, col.idx] = -1
    }
    
    # getSlopes(codacore0)
    # codacore0$ensemble[[1]]$intercept
    # codacore0$ensemble[[1]]$slope
    # codacore0 # the printed slope doesn't always match the slopes given above...
    codacore0_fit = glm(
      Y ~ getIlrX(X = X, sbp = codacore0_SBP), family = binomial(link = "logit"))
    # codacore0_fit
    # glm(Y ~ getLogRatios(codacore0, X), family = binomial(link = "logit"))
    
    
    codacore0.thetahat = coefficients(codacore0_fit)[-1] #########################
    U.codacore0 = getIlrTrans(sbp = codacore0_SBP)################################
    codacore0.betahat = U.codacore0 %*% as.matrix(codacore0.thetahat)
    
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
      thetahat = codacore0.thetahat, betahat = codacore0.betahat,
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
      "FP+" = 0, "FN+" = p, "TPR+" = 0, 
      "FP-" = p, "FN-" = 0, "TPR-" = 0
    )
  }
  
}

