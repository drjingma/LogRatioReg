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
  sigma.settings = "latentVarModel"
  n = 100
  p = 30
  K = 10
  nlam = 100
  neta = p
  intercept = TRUE
  scaling = TRUE
  tol = 1e-4
  sigma_eps1 = 0.1
  sigma_eps2 = 0.1
  SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
  ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
  # ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
  #   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
  b0 = 0 # 0
  b1 = 1 # 1, 0.5, 0.25
  theta.value = 1 # weight on a1 -- 1
  a0 = 0 # 0
  
  file.end = paste0(
    "_", sigma.settings,
    "_", paste0(
      paste(which(SBP.true == 1), collapse = ""), "v", 
      paste(which(SBP.true == -1), collapse = "")),
    "_dim", n, "x", p, 
    "_noisey", sigma_eps1, 
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
  U.all = matrix(runif(2 * n), ncol = 1)
  # simulate y from latent variable
  y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
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
  
  # ##############################################################################
  # # supervised log-ratios (a balance regression method)
  # #   -- hierarchical spectral clustering + thresholding with lasso
  # ##############################################################################
  # start.time = Sys.time()
  # # apply hierarchical spectral clustering to the SLR similarity matrix
  # slrSimMat = getSlrMatrix(
  #   y = Y, X = X, type = "similarity")
  # slrhsc_btree = HSClust(
  #   W = slrSimMat, force_levelMax = TRUE, method = "kmeans")
  # slrhsc_SBP = sbp.fromHSClust(
  #   levels_matrix = slrhsc_btree$allLevels, row_names = colnames(X))
  # # apply supervised log-ratios, using CV to select threshold and also lambda
  # slrhsc2 = cvBMLassoThresh(
  #   y = Y, X = X,
  #   W = slrSimMat, # normalized similarity matrix (all values between 0 & 1)
  #   hsc_method = "kmeans", # "shimalik", "kmeans"
  #   force_levelMax = TRUE,
  #   sbp = slrhsc_SBP,
  #   lambda = NULL, nlam = nlam,
  #   eta = NULL, neta = neta,
  #   nfolds = K, foldid = NULL,
  #   intercept = intercept,
  #   standardize = scaling
  # )
  # end.time = Sys.time()
  # slrhsc2.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # slrhsc2.eta.min.idx = slrhsc2$min.idx[2]
  # slrhsc2.lam.min.idx = slrhsc2$min.idx[1]
  # slrhsc2.a0 = slrhsc2$theta0[[slrhsc2.eta.min.idx]][slrhsc2.lam.min.idx]
  # slrhsc2.thetahat = slrhsc2$theta[[slrhsc2.eta.min.idx]][, slrhsc2.lam.min.idx]
  # slrhsc2.SBP = slrhsc2$sbp_thresh[[slrhsc2.eta.min.idx]]
  # slrhsc2.betahat.nonzero = getBetaFromTheta(slrhsc2.thetahat, sbp = slrhsc2.SBP)
  # slrhsc2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  # rownames(slrhsc2.betahat) = colnames(X)
  # slrhsc2.betahat[slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], ] =
  #   as.numeric(slrhsc2.betahat.nonzero)
  # 
  # # compute metrics on the selected model #
  # slrhsc2.metrics = getMetricsBalanceReg(
  #   y.train = Y, y.test = Y.test,
  #   ilrX.train = getIlrX(
  #     X[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
  #     sbp = slrhsc2.SBP),
  #   ilrX.test = getIlrX(
  #     X.test[, slrhsc2$meets_threshold[[slrhsc2.eta.min.idx]], drop = FALSE],
  #     sbp = slrhsc2.SBP),
  #   n.train = n, n.test = n,
  #   thetahat0 = slrhsc2.a0, thetahat = slrhsc2.thetahat,
  #   betahat = slrhsc2.betahat,
  #   true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
  #   true.beta = beta.true)
  # 
  # # # plot the tree given by slr-hsc, indicating significant covariates
  # # slrhsc2_leaf_types = rep("covariate", nrow(slrhsc2.SBP))
  # # slrhsc2_balance_types = rep("balance", ncol(slrhsc2.SBP))
  # # slrhsc2_nodes_types = data.frame(
  # #   name = c(colnames(slrhsc2.SBP), rownames(slrhsc2.SBP)),
  # #   type = c(slrhsc2_balance_types, slrhsc2_leaf_types)
  # # )
  # # plotSBP(slrhsc2.SBP, title = "slr-hsc-eta", nodes_types = slrhsc2_nodes_types)
  # # # fields::image.plot(slrSimMat)
  # 
  # saveRDS(c(
  #   slrhsc2.metrics,
  #   "betaSparsity" = bspars,
  # 
  #   "time" = slrhsc2.timing
  # ),
  # paste0(output_dir, "/metrics", "/slr_thresh_lasso_metrics", file.end))

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
    true.beta = beta.true)

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
  slrnew = slr(x = X, y = Y, rank1approx = TRUE)
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
    true.beta = beta.true)

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
  slrnew2 = slr(x = X, y = Y, rank1approx = FALSE)
  slrnew2_activevars = names(slrnew2$index)
  slrnew2_SBP = matrix(slrnew2$index)
  rownames(slrnew2_SBP) = slrnew2_activevars
  end.time = Sys.time()
  slrnew2.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")

  slrnew2.coefs = coefficients(slrnew2$model)
  slrnew2.a0 = slrnew2.coefs[1]
  slrnew2.thetahat = slrnew2.coefs[-1]

  slrnew2.betahat.nonzero = getBetaFromTheta(slrnew2.thetahat, sbp = slrnew2_SBP)
  slrnew2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
  rownames(slrnew2.betahat) = colnames(X)
  slrnew2.betahat[slrnew2_activevars, ] =
    as.numeric(slrnew2.betahat.nonzero)

  # compute metrics on the selected model #
  slrnew2.metrics = getMetricsBalanceReg(
    y.train = Y, y.test = Y.test,
    ilrX.train = getIlrX(X[, slrnew2_activevars, drop = FALSE], sbp = slrnew2_SBP),
    ilrX.test = getIlrX(X.test[, slrnew2_activevars, drop = FALSE], sbp = slrnew2_SBP),
    n.train = n, n.test = n,
    thetahat0 = slrnew2.a0, thetahat = slrnew2.thetahat, betahat = slrnew2.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true)

  saveRDS(c(
    slrnew2.metrics,
    "betaSparsity" = bspars,

    "time" = slrnew2.timing
  ),
  paste0(output_dir, "/metrics", "/slr_metrics", file.end))

  # ##############################################################################
  # # propr method (a balance regression method)
  # ##############################################################################
  # start.time = Sys.time()
  # pr_res <- suppressMessages(propr(X, metric = "phs"))
  # pr_btree = hclust(as.dist(pr_res@matrix),method = linkage)
  # pr_SBP = sbp.fromHclust(pr_btree)
  # pr = cvBMLasso(y = Y, X = X, sbp = pr_SBP, nlam = nlam,
  #                nfolds = K, intercept = intercept, standardize = scaling)
  # end.time = Sys.time()
  # pr.timing = difftime(
  #   time1 = end.time, time2 = start.time, units = "secs")
  # 
  # pr.lam.min.idx = which.min(pr$cvm)
  # pr.a0 = pr$theta0[pr.lam.min.idx]
  # pr.thetahat = pr$theta[, pr.lam.min.idx]
  # pr.betahat = getBetaFromTheta(pr.thetahat, sbp = pr$sbp)
  # 
  # # compute metrics on the selected model #
  # pr.metrics = getMetricsBalanceReg(
  #   y.train = Y, y.test = Y.test,
  #   ilrX.train = getIlrX(X, sbp = pr$sbp),
  #   ilrX.test = getIlrX(X.test, sbp = pr$sbp),
  #   n.train = n, n.test = n,
  #   thetahat0 = pr.a0, thetahat = pr.thetahat, betahat = pr.betahat,
  #   true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
  #   true.beta = beta.true)
  # 
  # # # plot the tree given by slr-hsc, indicating significant covariates
  # # pr_leaf_types = rep("covariate", nrow(pr$sbp))
  # # pr_balance_types = rep("balance", ncol(pr$sbp))
  # # pr_nodes_types = data.frame(
  # #   name = c(colnames(pr$sbp), rownames(pr$sbp)),
  # #   type = c(pr_balance_types, pr_leaf_types)
  # # )
  # # plotSBP(pr$sbp, title = "propr", nodes_types = pr_nodes_types)
  # # fields::image.plot(pr_res@matrix)
  # 
  # saveRDS(c(
  #   pr.metrics,
  #   "betaSparsity" = bspars,
  # 
  #   "time" = pr.timing
  # ),
  # paste0(output_dir, "/metrics", "/propr_metrics", file.end))

  ##############################################################################
  # selbal method (a balance regression method)
  ##############################################################################
  library(selbal) # masks stats::cor()

  start.time = Sys.time()
  X.slbl = X
  rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
  colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
  slbl = selbal.cv(x = X.slbl, y = as.vector(Y), n.fold = K)
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
  # lm(as.vector(Y) ~ log(X) %*% as.matrix(U.slbl))
  # slbl$glm
  slbl.thetahat = coefficients(slbl$glm)[-1]
  slbl.betahat = U.slbl %*% as.matrix(slbl.thetahat)

  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  slbl.Yhat.train = predict.glm(
    slbl$glm, newdata = data.frame(X.slbl), type = "response")
  slbl.PE.train = crossprod(Y - slbl.Yhat.train) / n
  # get prediction error on test set
  X.slbl.test = X.test
  rownames(X.slbl.test) = paste("Sample", 1:nrow(X.slbl.test), sep = "_")
  colnames(X.slbl.test) = paste("V", 1:ncol(X.slbl.test), sep = "")
  slbl.Yhat.test = predict.glm(
    slbl$glm, newdata = data.frame(X.slbl.test), type = "response")
  slbl.PE.test = crossprod(Y.test - slbl.Yhat.test) / n
  # beta estimation accuracy, selection accuracy #
  slbl.metrics = getMetricsBalanceReg(
    thetahat = slbl.thetahat, betahat = slbl.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  slbl.metrics = c(PEtr = slbl.PE.train, PEte = slbl.PE.test, slbl.metrics)
  
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
    objective = "regression")
  codacore0_SBP = matrix(0, nrow = p, ncol = length(codacore0$ensemble))
  for(col.idx in 1:ncol(codacore0_SBP)){
    codacore0_SBP[codacore0$ensemble[[col.idx]]$hard$numerator, col.idx] = 1
    codacore0_SBP[codacore0$ensemble[[col.idx]]$hard$denominator, col.idx] = -1
  }
  end.time = Sys.time()
  codacore0.timing = difftime(
    time1 = end.time, time2 = start.time, units = "secs")
  
  # getSlopes(codacore0)
  # codacore0$ensemble[[1]]$intercept
  # codacore0$ensemble[[1]]$slope
  # codacore0 # the printed slope doesn't match the slopes given above...
  codacore0_fit = lm(Y ~ getIlrX(X = X, sbp = codacore0_SBP))
  # codacore0_fit
  # lm(Y ~ getLogRatios(codacore0, X))
  # lm(Y ~ -1 + getLogRatios(codacore0, X))
  
  
  codacore0.thetahat = coefficients(codacore0_fit)[-1] #########################
  U.codacore0 = getIlrTrans(sbp = codacore0_SBP)################################
  codacore0.betahat = U.codacore0 %*% as.matrix(codacore0.thetahat)
  
  # compute metrics on the selected model #
  # prediction errors
  # get prediction error on training set
  codacore0.Yhat.train = predict(codacore0, X)
  codacore0.PE.train = crossprod(Y - codacore0.Yhat.train) / n
  # get prediction error on test set
  codacore0.Yhat.test = predict(codacore0, X.test)
  codacore0.PE.test = crossprod(Y.test - codacore0.Yhat.test) / n
  # beta estimation accuracy, selection accuracy #
  codacore0.metrics = getMetricsBalanceReg(
    thetahat = codacore0.thetahat, betahat = codacore0.betahat,
    true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
    true.beta = beta.true, metrics = c("betaestimation", "selection"))
  codacore0.metrics = c(PEtr = codacore0.PE.train, PEte = codacore0.PE.test, codacore0.metrics)
  
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


