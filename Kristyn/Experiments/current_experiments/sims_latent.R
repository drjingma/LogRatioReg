# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 2/21/2021

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
  
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  source("Kristyn/Functions/supervisedlogratioseta.R")
  source("Kristyn/Functions/HSClust.R")
  source("Kristyn/Functions/slrnew.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  source("Kristyn/Functions/simulatedata.R")
  
  # for plots
  library(ggraph) # make dendrogram
  library(igraph) # transform dataframe to graph object: graph_from_data_frame()
  library(tidygraph)
  
  # Tuning parameters###########################################################
  
  # Settings to toggle with
  sigma.settings = "diagSigma"
  theta.value = 1
  intercept = TRUE
  K = 10
  n = 100
  p = 30
  scaling = TRUE
  linkage = "average"
  tol = 1e-4
  nlam = 100
  neta = p
  #################
  # SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  #################
  rho = 0.2 #
  desired_Rsquared = 0.6 #
  # sigma_eps1 = get_sigma_eps(
  #   sbp = SBP.true, ilr.trans.constant = ilrtrans.true$const, theta = theta.value, 
  #   Rsq = desired_Rsquared, rho = rho)
  sigma_eps1 = 0.5
  sigma_eps2 = 0.1
  #################
  
  # Population parameters
  # SigmaW = rgExpDecay(p, rho)$Sigma
  # muW = c(rep(log(p), 5), rep(0, p - 5))
  # names(muW) = paste0('s', 1:p)
  
  file.end = paste0(
    "_", sigma.settings,
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_dim", n, "x", p, 
    # "_rho", rho, 
    "_noisey", sigma_eps1, 
    "_noisex", sigma_eps2, 
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # generate data
  # get info from sbp
  num.pos = ilrtrans.true$pos.vec
  num.neg = ilrtrans.true$neg.vec
  # get latent variable
  latent.var.all = runif(2 * n)
  # simulate y from latent variable
  b0 = 1
  y.all = b0 + theta.value * latent.var.all + rnorm(2 * n) * sigma_eps1
  # simulate X
  alrX.all.noises = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
  # alrX.covariates = as.matrix(
  #   c(rep(1, num.pos) / num.pos, 
  #     -rep(1, num.neg) / num.neg, 
  #     rep(0, p - 1 - num.pos - num.neg)))
  alrX.all.covariates = ilrtrans.true$ilr.trans[-p]
  coeffX1 = 1 # 10
  alrX.all = as.matrix(latent.var.all) %*% (coeffX1 * t(alrX.all.covariates)) + 
    alrX.all.noises
  X.all <- alrinv(alrX.all)
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
  
  ##############################################################################
  # supervised log-ratios (a balance regression method)
  #   -- hierarchical spectral clustering + thresholding with lasso
  ##############################################################################
  start.time = Sys.time()
  # apply hierarchical spectral clustering to the SLR similarity matrix
  slrSimMat = getSlrMatrix(
    y = Y, X = X, type = "similarity")
  slrhsc_btree = HSClust(
    W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  slrhsc_SBP = sbp.fromHSClust(
    levels_matrix = slrhsc_btree$allLevels, row_names = colnames(X))
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrhsc2 = cvBMLassoThresh(
    y = Y, X = X,
    W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
    hsc_method = "kmeans", # "shimalik", "kmeans"
    force_levelMax = TRUE,
    sbp = slrhsc_SBP,
    lambda = NULL, nlam = nlam,
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL,
    intercept = intercept,
    standardize = scaling
  )
  end.time = Sys.time()
  slrhsc2.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrhsc2.eta.min.idx = slrhsc2$min.idx[2]
  slrhsc2.lam.min.idx = slrhsc2$min.idx[1]
  slrhsc2.a0 = slrhsc2$theta0[[slrhsc2.eta.min.idx]][slrhsc2.lam.min.idx]
  slrhsc2.thetahat = slrhsc2$theta[[slrhsc2.eta.min.idx]][, slrhsc2.lam.min.idx]
  slrhsc2.SBP = slrhsc2$sbp_thresh[[slrhsc2.eta.min.idx]]
  slrhsc2.betahat.nonzero = getBetaFromTheta(slrhsc2.thetahat, sbp = slrhsc2.SBP)
  slrhsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrhsc2.betahat) = colnames(X)
  slrhsc2.betahat[slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], ] =
    as.numeric(slrhsc2.betahat.nonzero)

  # compute metrics on the selected model #
  slrhsc2.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(
      X[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
      sbp = slrhsc2.SBP),
    ilrX.test = getIlrX(
      X.test[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
      sbp = slrhsc2.SBP),
    n.train = n, n.test = n,
    thetahat0 = slrhsc2.a0, thetahat = slrhsc2.thetahat,
    betahat = slrhsc2.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    metrics = c("prediction", "selection"))

  # # plot the tree given by slr-hsc, indicating significant covariates
  # slrhsc2_leaf_types = rep("covariate", nrow(slrhsc2.SBP))
  # slrhsc2_balance_types = rep("balance", ncol(slrhsc2.SBP))
  # slrhsc2_nodes_types = data.frame(
  #   name = c(colnames(slrhsc2.SBP), rownames(slrhsc2.SBP)),
  #   type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
  # )
  # plotSBP(slrhsc2.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_nodes_types)
  # # fields::image.plot(slrSimMat)

  saveRDS(c(
    slrhsc2.metrics,
    "betaSparsity" = bspars,

    "time" = slrhsc2.timing
  ),
  paste0(output_dir, "/metrics", "/slr_hsc_thresh_lasso_metrics", file.end))

  ##############################################################################
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  start.time = Sys.time()
  classo = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  end.time = Sys.time()
  cl.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  cl.lam.min.idx = which.min(classo$cvm)
  cl.a0 = classo$int[cl.lam.min.idx]
  cl.betahat = classo$bet[, cl.lam.min.idx]

  # compute metrics on the selected model #
  cl.metrics = getMetricsLLC(
    y.train = Y, y.test = Y.test,
    logX.train = log(X),
    logX.test = log(X.test),
    n.train = n, n.test = n,
    betahat0 = cl.a0, betahat = cl.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta, 
    metrics = c("prediction", "selection"))

  saveRDS(c(
    cl.metrics,
    "betaSparsity" = bspars,
     
    "time" = cl.timing
  ),
  paste0(output_dir, "/metrics", "/classo_metrics", file.end))

  ##############################################################################
  # new slr method (a balance regression method)
  #   -- hierarchical spectral clustering
  ##############################################################################
  start.time = Sys.time()
  slrnew = slr(x = X, y = Y)
  slrnew_activevars = names(slrnew$index)
  slrnew_SBP = matrix(slrnew$index)
  rownames(slrnew_SBP) = slrnew_activevars
  end.time = Sys.time()
  slrnew.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrnew.coefs = coefficients(slrnew$model)
  slrnew.a0 = slrnew.coefs[1]
  slrnew.thetahat = slrnew.coefs[-1]

  slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
  slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrnew.betahat) = colnames(X)
  slrnew.betahat[slrnew_activevars, ] =
    as.numeric(slrnew.betahat.nonzero)

  # compute metrics on the selected model #
  slrnew.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
    ilrX.test = getIlrX(X.test[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrnew.a0, thetahat = slrnew.thetahat, betahat = slrnew.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    metrics = c("prediction", "selection"))

  saveRDS(c(
    slrnew.metrics,
    "betaSparsity" = bspars,

    "time" = slrnew.timing
  ),
  paste0(output_dir, "/metrics", "/slr_new_metrics", file.end))

  ##############################################################################
  # new slr method (a balance regression method)
  #   -- hierarchical spectral clustering (without rank 1 approximation)
  ##############################################################################
  start.time = Sys.time()
  slrnew = slr(x = X, y = Y, rank1approx = FALSE)
  slrnew_activevars = names(slrnew$index)
  slrnew_SBP = matrix(slrnew$index)
  rownames(slrnew_SBP) = slrnew_activevars
  end.time = Sys.time()
  slrnew.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrnew.coefs = coefficients(slrnew$model)
  slrnew.a0 = slrnew.coefs[1]
  slrnew.thetahat = slrnew.coefs[-1]

  slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
  slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrnew.betahat) = colnames(X)
  slrnew.betahat[slrnew_activevars, ] =
    as.numeric(slrnew.betahat.nonzero)

  # compute metrics on the selected model #
  slrnew.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
    ilrX.test = getIlrX(X.test[, slrnew_activevars, drop = FALSE], sbp = slrnew_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrnew.a0, thetahat = slrnew.thetahat, betahat = slrnew.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    metrics = c("prediction", "selection"))

  saveRDS(c(
    slrnew.metrics,
    "betaSparsity" = bspars,

    "time" = slrnew.timing
  ),
  paste0(output_dir, "/metrics", "/slr_new_noapprox_metrics", file.end))

  ##############################################################################
  # propr method (a balance regression method)
  ##############################################################################
  start.time = Sys.time()
  pr_res <- suppressMessages(propr(X, metric = "phs"))
  pr_btree = hclust(as.dist(pr_res@matrix),method = linkage)
  pr_SBP = sbp.fromHclust(pr_btree)
  pr = cvBMLasso(y = Y, X = X, sbp = pr_SBP, nlam = nlam,
                 nfolds = K, intercept = intercept, standardize = scaling)
  end.time = Sys.time()
  pr.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  pr.lam.min.idx = which.min(pr$cvm)
  pr.a0 = pr$theta0[pr.lam.min.idx]
  pr.thetahat = pr$theta[, pr.lam.min.idx]
  pr.betahat = getBetaFromTheta(pr.thetahat, sbp = pr$sbp)

  # compute metrics on the selected model #
  pr.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X, sbp = pr$sbp),
    ilrX.test = getIlrX(X.test, sbp = pr$sbp),
    n.train = n, n.test = n,
    thetahat0 = pr.a0, thetahat = pr.thetahat, betahat = pr.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    metrics = c("prediction", "selection"))

  # # plot the tree given by slr-hsc, indicating significant covariates
  # pr_leaf_types = rep("covariate", nrow(pr$sbp))
  # pr_balance_types = rep("balance", ncol(pr$sbp))
  # pr_nodes_types = data.frame(
  #   name = c(colnames(pr$sbp), rownames(pr$sbp)),
  #   type = c(pr_balance_types, pr_leaf_types)
  # )
  # plotSBP(pr$sbp, title = "propr", nodes_types = pr_nodes_types)
  # fields::image.plot(pr_res@matrix)

  saveRDS(c(
    pr.metrics,
    "betaSparsity" = bspars,

    "time" = pr.timing
  ),
  paste0(output_dir, "/metrics", "/propr_metrics", file.end))

  # ##############################################################################
  # # selbal method (a balance regression method)
  # ##############################################################################
  # library(selbal) # masks stats::cor()
  # 
  # start.time = Sys.time()
  # X.slbl = X
  # rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
  # colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
  # slbl = selbal.cv(x = X.slbl, y = as.vector(Y), n.fold = K)
  # end.time = Sys.time()
  # slbl.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # U (transformation) matrix
  # U.slbl = rep(0, p)
  # names(U.slbl) = colnames(X.slbl)
  # pba.pos = unlist(subset(
  #   slbl$global.balance, subset = Group == "NUM", select = Taxa))
  # num.pos = length(pba.pos)
  # pba.neg = unlist(subset(
  #   slbl$global.balance, subset = Group == "DEN", select = Taxa))
  # num.neg = length(pba.neg)
  # U.slbl[pba.pos] = 1 / num.pos
  # U.slbl[pba.neg] = -1 / num.neg
  # norm.const = sqrt((num.pos * num.neg) / (num.pos + num.neg))
  # U.slbl = norm.const * U.slbl
  # # check: these are equal
  # # lm(as.vector(Y) ~ log(X) %*% as.matrix(U.slbl))
  # # slbl$glm
  # slbl.thetahat = coefficients(slbl$glm)[2]
  # slbl.betahat = U.slbl %*% as.matrix(slbl.thetahat)
  # 
  # # compute metrics on the selected model #
  # # prediction errors
  # # get prediction error on training set
  # slbl.Yhat.train = predict.glm(
  #   slbl$glm, newdata = data.frame(X.slbl), type = "response")
  # slbl.PE.train = crossprod(Y - slbl.Yhat.train) / n
  # # get prediction error on test set
  # X.slbl.test = X.test
  # rownames(X.slbl.test) = paste("Sample", 1:nrow(X.slbl.test), sep = "_")
  # colnames(X.slbl.test) = paste("V", 1:ncol(X.slbl.test), sep = "")
  # slbl.Yhat.test = predict.glm(
  #   slbl$glm, newdata = data.frame(X.slbl.test), type = "response")
  # slbl.PE.test = crossprod(Y.test - slbl.Yhat.test) / n
  # # selection accuracy
  # slbl.non0.betahat = abs(slbl.betahat) > 10e-8
  # pseudo.true.beta = as.numeric(apply(true.sbp, 1, function(row) any(row != 0)))
  # slbl.SA = getSelectionAccuracy(
  #   is0.true.beta = is0.beta, non0.true.beta = non0.beta,
  #   non0.betahat = slbl.non0.betahat)
  # 
  # slbl.metrics = c(
  #   "PEtr" = slbl.PE.train,
  #   "PEte" = slbl.PE.test,
  #   "FP" = slbl.SA$FP,
  #   "FN" = slbl.SA$FN,
  #   "TPR" = slbl.SA$TPR,
  #   "precision" = slbl.SA$precision,
  #   "Fscore" = slbl.SA$Fscore
  # )
  # saveRDS(c(
  #   slbl.metrics,
  #   "betaSparsity" = bspars,
  #    
  #   "time" = slbl.timing
  # ),
  # paste0(output_dir, "/metrics", "/selbal_metrics", file.end))
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ### fin ###
}


