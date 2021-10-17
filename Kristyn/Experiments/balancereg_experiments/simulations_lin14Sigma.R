# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Note: Here, we simulate from the Sigma specified in Lin et al.2014 and intend 
#   to:
#   - observe the behaviors of slr under different measurement noise, sigma_eps
#   - compare to the compositional lasso model proposed in lin et al. 2014
# Date: 10/12/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_experiments/outputs"

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
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
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  
  library(balance) # for sbp.fromHclust()
  library(selbal)
  library(propr)
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  source("Kristyn/Functions/simulatedata.R")
  source("Kristyn/Functions/bic.R")
  
  # Settings to toggle with
  sigma.settings = "lin14Sigma"
  theta.settings = "dense" 
  # "dense" => j = 1 (of theta = theta_1, ..., theta_j, ..., theta_{p-1})
  # "sparse" = c(p - 3:1)
  values.theta = 1
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  linkage = "average"
  tol = 1e-4
  nlam = 100
  intercept = TRUE
  K = 10
  n = 100
  p = 200
  rho = 0.2 # 0.2, 0.5
  scaling = TRUE
  
  # Population parameters
  sigma_eps = 0.1 # 0.01, 0.1
  SigmaW <- rgExpDecay(p,rho)$Sigma
  SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
  U = getU(btree = SigmaWtree) # transformation matrix
  
  # theta settings
  if(theta.settings == "dense"){
    indices.theta = 1
    # } else if(theta.settings == "sparse"){
    #   indices.theta = p - 1
    # } else if(theta.settings == "both"){
    #   indices.theta = c(1, p - 1)
  } else{ # if(theta.settings == "multsparse")
    indices.theta = p - 3:1
  }
  # error checking indices.theta found based on theta.settings argument
  if(is.null(indices.theta)){
    stop("invalid indices.theta")
  }
  
  # get theta
  theta = rep(0, p - 1)
  if(is.null(values.theta)){
    # assume values.theta = 1
    values.theta = 1
  }
  if(length(values.theta) == 1){
    # if values.theta = 1, assume it's the same value for all nonzero indices
    values.theta = rep(values.theta, length(indices.theta))
  } else if(length(values.theta) != length(indices.theta)){
    # when 1 < length(values.theta) < total # of nonzero values
    stop("indices.theta does not have same length as values.theta")
  }
  theta[indices.theta] = values.theta
  theta = as.matrix(theta)
  
  # about beta
  beta = as.vector(getBeta(theta, U = U))
  names(beta) <- paste0('s', 1:p)
  non0.beta = (beta != 0)
  slr.non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
  # Population parameters, continued
  muW = c(rep(log(p), 5), rep(0, p - 5))
  names(muW) = names(beta)
  
  file.end = paste0(
    "_", sigma.settings,
    "_", theta.settings, 
    "_dim", n, "x", p, 
    "_noise", sigma_eps,
    "_rho", rho, 
    "_int", intercept,
    "_scale", scaling,
    "_sim", b,
    ".rds")
  
  ##############################################################################
  # simulate data
  fake.data = simulateBalanceReg(
    mu = muW, Sigma = SigmaW, U = U, n = 2 * n, theta = theta, 
    sigma.noise = sigma_eps)
  colnames(fake.data$X) = names(beta)
  
  # subset out training and test sets
  X = fake.data$X[1:n, ]
  X.test = fake.data$X[-(1:n), ]
  Y <- fake.data$y[1:n, , drop = TRUE]
  Y.test <- fake.data$y[-(1:n), , drop = TRUE]
  
  ##############################################################################
  # supervised log-ratios (a balance regression method)
  ##############################################################################
  
  # fit model ##################################################################
  if(!file.exists(paste0(output_dir, "/models", "/slr_model", file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/slr_timing", file.end))){
    slr.model.already.existed = FALSE
    # apply supervised log-ratios, using CV to select lambda
    start.time = Sys.time()
    slr = cvSLR(
      y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
      rho.type = rho.type, linkage = linkage, standardize = scaling)
    end.time = Sys.time()
    saveRDS(slr, paste0(output_dir, "/models", "/slr_model", file.end))
    
    # timing metric
    slr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      slr.timing, 
      paste0(output_dir, "/timing", "/slr_timing", file.end))
  } else{
    slr.model.already.existed = TRUE
  }
  
  # select tuning parameter and calculate metrics ##############################
  if(!file.exists(paste0(
    output_dir, "/metrics", "/slr_metrics", file.end)) | 
    !file.exists(paste0(
      output_dir, "/metrics", "/slr_bic_metrics", file.end)) | 
    slr.model.already.existed == FALSE){
    
    # import model and timing metric
    slr = readRDS(paste0(output_dir, "/models", "/slr_model", file.end))
    slr.timing = readRDS(paste0(
      output_dir, "/timing", "/slr_timing", file.end))
    slr.btree = slr$btree
    slr.btree = slr$btree
    slr.SBP = sbp.fromHclust(slr.btree)
    
    # choose lambda using cross-validated mse ##################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/slr_metrics", file.end))){
      slr.lam.min.idx = which.min(slr$cvm)
      slr.a0 = slr$int[slr.lam.min.idx]
      slr.thetahat = slr$bet[, slr.lam.min.idx]
      slr.Uhat = getU(btree = slr.btree)
      slr.betahat = getBeta(slr.thetahat, U = slr.Uhat)
      
      # compute metrics on the selected model #
      slr.metrics = getMetricsBalanceReg(
        y.train = Y, y.test = Y.test, 
        ilrX.train = computeBalances(X, slr.btree), 
        ilrX.test = computeBalances(X.test, slr.btree), 
        n.train = n, n.test = n, 
        thetahat0 = slr.a0, thetahat = slr.thetahat, betahat = slr.betahat, 
        sbp = slr.SBP, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        slr.metrics, 
        "timing" = slr.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/slr_metrics", file.end))
    }
    # choose lambda using bic ##################################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/slr_bic_metrics", file.end))){
      slr.bic = getBICseq(
        y = Y, predMat = computeBalances(X, slr.btree), 
        betahat0.vec = slr$int, betahat.mat = slr$bet)
      slr.lam.min.idx = which.min(slr.bic)
      slr.a0 = slr$int[slr.lam.min.idx]
      slr.thetahat = slr$bet[, slr.lam.min.idx]
      slr.Uhat = getU(btree = slr.btree)
      slr.betahat = getBeta(slr.thetahat, U = slr.Uhat)
      
      # compute metrics on the selected model #
      slr.metrics = getMetricsBalanceReg(
        y.train = Y, y.test = Y.test, 
        ilrX.train = computeBalances(X, slr.btree), 
        ilrX.test = computeBalances(X.test, slr.btree), 
        n.train = n, n.test = n, 
        thetahat0 = slr.a0, thetahat = slr.thetahat, betahat = slr.betahat, 
        sbp = slr.SBP, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        slr.metrics, 
        "timing" = slr.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/slr_bic_metrics", file.end))
    }
  }
  
  # roc curves #################################################################
  if(!file.exists(paste0(output_dir, "/roccurves", "/slr_roc", file.end)) | 
     slr.model.already.existed == FALSE){
    
    # import model
    slr = readRDS(paste0(output_dir, "/models", "/slr_model", file.end))
    slr.btree = slr$btree
    slr.btree = slr$btree
    slr.SBP = sbp.fromHclust(slr.btree)
    
    # roc
    slr.roc <- apply(slr$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, slr.SBP))
    
    saveRDS(slr.roc, paste0(output_dir, "/roccurves", "/slr_roc", file.end))
  }
  
  ##############################################################################
  # compositional lasso (a linear log contrast method)
  ##############################################################################
  
  # fit model ##################################################################
  if(!file.exists(paste0(output_dir, "/models", "/classo_model", file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/classo_timing", file.end))){
    cl.model.already.existed = FALSE
    # apply compositional lasso, using CV to select lambda
    start.time = Sys.time()
    classo = cv.func(
      method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
      nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
    end.time = Sys.time()
    saveRDS(classo, paste0(output_dir, "/models", "/classo_model", file.end))
    
    # timing metric
    cl.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      cl.timing, 
      paste0(output_dir, "/timing", "/classo_timing", file.end))
  } else{
    cl.model.already.existed = TRUE
  }
  
  # select tuning parameter and calculate metrics ##############################
  if(!file.exists(paste0(
    output_dir, "/metrics", "/classo_metrics", file.end)) | 
    !file.exists(paste0(
      output_dir, "/metrics", "/classo_bic_metrics", file.end)) |
    cl.model.already.existed == FALSE){
    
    # import model and timing metric
    classo = readRDS(paste0(output_dir, "/models", "/classo_model", file.end))
    cl.timing = readRDS(paste0(
      output_dir, "/timing", "/classo_timing", file.end))
    
    # choose lambda using cross-validated mse ##################################
    if(!file.exists(paste0(output_dir, "/metrics", "/classo_metrics", file.end))){
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
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        cl.metrics, 
        "timing" = cl.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/classo_metrics", file.end))
    }
    # choose lambda using bic ##################################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/classo_bic_metrics", file.end))){
      cl.bic = getBICseq(
        y = Y, predMat = log(X), 
        betahat0.vec = classo$int, betahat.mat = classo$bet)
      cl.lam.min.idx = which.min(cl.bic)
      cl.a0 = classo$int[cl.lam.min.idx]
      cl.betahat = classo$bet[, cl.lam.min.idx]
      
      # compute metrics on the selected model #
      cl.metrics = getMetricsLLC(
        y.train = Y, y.test = Y.test, 
        logX.train = log(X), 
        logX.test = log(X.test), 
        n.train = n, n.test = n, 
        betahat0 = cl.a0, betahat = cl.betahat, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        cl.metrics, 
        "timing" = cl.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/classo_bic_metrics", file.end))
    }
  }
  
  # roc curves #################################################################
  if(!file.exists(paste0(output_dir, "/roccurves", "/classo_roc", file.end)) | 
     cl.model.already.existed == FALSE){
    
    # import model
    classo = readRDS(paste0(output_dir, "/models", "/classo_model", file.end))
    
    # roc
    cl.roc <- apply(classo$bet, 2, function(a) 
      roc.for.coef(a, beta))
    
    saveRDS(cl.roc, paste0(output_dir, "/roccurves", "/classo_roc", file.end))
  }
  
  ##############################################################################
  # oracle method (a balance regression method)
  ##############################################################################
  
  # fit model ##################################################################
  if(!file.exists(paste0(output_dir, "/models", "/oracle_model", file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/oracle_timing", file.end))){
    or.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    start.time = Sys.time()
    oracle = cvILR(y = Y, X = X, btree = SigmaWtree, U = U, nlam = nlam, 
                   nfolds = K, intercept = intercept, standardize = scaling)
    end.time = Sys.time()
    saveRDS(oracle, paste0(output_dir, "/models", "/oracle_model", file.end))
    
    # timing metric
    or.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      or.timing, 
      paste0(output_dir, "/timing", "/oracle_timing", file.end))
  } else{
    or.model.already.existed = TRUE
  }
  
  # select tuning parameter and calculate metrics ##############################
  if(!file.exists(paste0(
    output_dir, "/metrics", "/oracle_metrics", file.end)) | 
    !file.exists(paste0(
      output_dir, "/metrics", "/oracle_bic_metrics", file.end)) |
    or.model.already.existed == FALSE){
    
    # import model and timing metric
    oracle = readRDS(paste0(output_dir, "/models", "/oracle_model", file.end))
    or.timing = readRDS(paste0(
      output_dir, "/timing", "/oracle_timing", file.end))
    or.btree = SigmaWtree
    or.SBP = sbp.fromHclust(or.btree)
    rownames(or.SBP) = names(beta)
    
    # choose lambda using cross-validated mse ##################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/oracle_metrics", file.end))){
      or.lam.min.idx = which.min(oracle$cvm)
      or.a0 = oracle$int[or.lam.min.idx]
      or.thetahat = oracle$bet[, or.lam.min.idx]
      or.Uhat = getU(btree = or.btree)
      or.betahat = getBeta(or.thetahat, U = or.Uhat)
      
      # compute metrics on the selected model #
      or.metrics = getMetricsBalanceReg(
        y.train = Y, y.test = Y.test, 
        ilrX.train = computeBalances(X, or.btree), 
        ilrX.test = computeBalances(X.test, or.btree), 
        n.train = n, n.test = n, 
        thetahat0 = or.a0, thetahat = or.thetahat, betahat = or.betahat, 
        sbp = or.SBP, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        or.metrics, 
        "timing" = or.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/oracle_metrics", file.end))
    }
    # choose lambda using bic ##################################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/oracle_bic_metrics", file.end))){
      or.bic = getBICseq(
        y = Y, predMat = computeBalances(X, or.btree), 
        betahat0.vec = oracle$int, betahat.mat = oracle$bet)
      or.lam.min.idx = which.min(or.bic)
      or.a0 = oracle$int[or.lam.min.idx]
      or.thetahat = oracle$bet[, or.lam.min.idx]
      or.Uhat = getU(btree = or.btree)
      or.betahat = getBeta(or.thetahat, U = or.Uhat)
      
      # compute metrics on the selected model #
      or.metrics = getMetricsBalanceReg(
        y.train = Y, y.test = Y.test, 
        ilrX.train = computeBalances(X, or.btree), 
        ilrX.test = computeBalances(X.test, or.btree), 
        n.train = n, n.test = n, 
        thetahat0 = or.a0, thetahat = or.thetahat, betahat = or.betahat, 
        sbp = or.SBP, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        or.metrics, 
        "timing" = or.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/oracle_bic_metrics", file.end))
    }
  } 
  
  # roc curves #################################################################
  if(!file.exists(paste0(output_dir, "/roccurves", "/oracle_roc", file.end)) | 
     or.model.already.existed == FALSE){
    
    # import model
    oracle = readRDS(paste0(output_dir, "/models", "/oracle_model", file.end))
    or.btree = SigmaWtree
    or.SBP = sbp.fromHclust(or.btree)
    rownames(or.SBP) = names(beta)
    
    # roc
    or.roc <- apply(oracle$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, or.SBP))
    
    saveRDS(or.roc, paste0(output_dir, "/roccurves", "/oracle_roc", file.end))
  }
  
  ##############################################################################
  # propr method (a balance regression method)
  ##############################################################################
  
  # fit model ##################################################################
  if(!file.exists(paste0(output_dir, "/models", "/propr_model", file.end)) |
     !file.exists(paste0(output_dir, "/timing", "/propr_timing", file.end))){
    pr.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    start.time = Sys.time()
    pr <- suppressMessages(propr(X, metric = "phs"))
    pr.tree = hclust(as.dist(pr@matrix),method = linkage)
    pr = cvILR(y = Y, X = X, btree = pr.tree, nlam = nlam, 
               nfolds = K, intercept = intercept, standardize = scaling)
    end.time = Sys.time()
    saveRDS(pr, paste0(output_dir, "/models", "/propr_model", file.end))
    
    # timing metric
    pr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
    saveRDS(
      pr.timing, 
      paste0(output_dir, "/timing", "/propr_timing", file.end))
  } else{
    pr.model.already.existed = TRUE
  }
  
  # select tuning parameter and calculate metrics ##############################
  if(!file.exists(paste0(
    output_dir, "/metrics", "/propr_metrics", file.end)) | 
    !file.exists(paste0(
      output_dir, "/metrics", "/propr_bic_metrics", file.end)) |
    pr.model.already.existed == FALSE){
    
    # import model and timing metric
    pr = readRDS(paste0(output_dir, "/models", "/propr_model", file.end))
    pr.timing = readRDS(paste0(
      output_dir, "/timing", "/propr_timing", file.end))
    pr.btree = pr$btree
    pr.SBP = sbp.fromHclust(pr.btree)
    
    # choose lambda using cross-validated mse ##################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/propr_metrics", file.end))){
      pr.lam.min.idx = which.min(pr$cvm)
      pr.a0 = pr$int[pr.lam.min.idx]
      pr.thetahat = pr$bet[, pr.lam.min.idx]
      pr.Uhat = getU(btree = pr.btree)
      pr.betahat = getBeta(pr.thetahat, U = pr.Uhat)
      
      # compute metrics on the selected model #
      pr.metrics = getMetricsBalanceReg(
        y.train = Y, y.test = Y.test, 
        ilrX.train = computeBalances(X, pr.btree), 
        ilrX.test = computeBalances(X.test, pr.btree), 
        n.train = n, n.test = n, 
        thetahat0 = pr.a0, thetahat = pr.thetahat, betahat = pr.betahat, 
        sbp = pr.SBP, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        pr.metrics, 
        "timing" = pr.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/propr_metrics", file.end))
    }
    # choose lambda using bic ##################################################
    if(!file.exists(paste0(
      output_dir, "/metrics", "/propr_bic_metrics", file.end))){
      pr.bic = getBICseq(
        y = Y, predMat = computeBalances(X, pr.btree), 
        betahat0.vec = pr$int, betahat.mat = pr$bet)
      pr.lam.min.idx = which.min(pr.bic)
      pr.a0 = pr$int[pr.lam.min.idx]
      pr.thetahat = pr$bet[, pr.lam.min.idx]
      pr.Uhat = getU(btree = pr.btree)
      pr.betahat = getBeta(pr.thetahat, U = pr.Uhat)
      
      # compute metrics on the selected model #
      pr.metrics = getMetricsBalanceReg(
        y.train = Y, y.test = Y.test, 
        ilrX.train = computeBalances(X, pr.btree), 
        ilrX.test = computeBalances(X.test, pr.btree), 
        n.train = n, n.test = n, 
        thetahat0 = pr.a0, thetahat = pr.thetahat, betahat = pr.betahat, 
        sbp = pr.SBP, 
        true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
      
      saveRDS(c(
        pr.metrics, 
        "timing" = pr.timing,
        "betaSparsity" = bspars
      ), 
      paste0(output_dir, "/metrics", "/propr_bic_metrics", file.end))
    }
  }
  
  # roc curves #################################################################
  if(!file.exists(paste0(output_dir, "/roccurves", "/propr_roc", file.end)) | 
     pr.model.already.existed == FALSE){
    
    # import model
    pr = readRDS(paste0(output_dir, "/models", "/propr_model", file.end))
    pr.btree = pr$btree
    pr.SBP = sbp.fromHclust(pr.btree)
    
    # roc
    pr.roc <- apply(pr$bet, 2, function(a) 
      roc.for.coef.LR(a, beta, pr.SBP))
    
    saveRDS(pr.roc, paste0(output_dir, "/roccurves", "/propr_roc", file.end))
  }
  
  ##############################################################################
  # full linear log-contrast model
  ##############################################################################
  # went with linear log-contrast model, because then I don't have to decide
  #   which U transformation matrix to use for a balance regression model
  
  if(!file.exists(paste0(
    output_dir, "/metrics", "/full_metrics", file.end))){
    Q = as.matrix(rep(1, p))
    Q2 = rbind(0, Q)
    full = lsei(A = cbind(1, log(X)), B = Y, E = t(Q2), F = 0)
    full.a0 = full$X[1]
    full.betahat = full$X[-1]
    
    # compute metrics on the selected model #
    full.metrics = getMetricsLLC(
      y.train = Y, y.test = Y.test, 
      logX.train = log(X), 
      logX.test = log(X.test), 
      n.train = n, n.test = n, 
      betahat0 = full.a0, betahat = full.betahat, 
      true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)
    
    saveRDS(c(
      full.metrics, 
      "timing" = cl.timing,
      "betaSparsity" = bspars
    ), 
    paste0(output_dir, "/metrics", "/full_metrics", file.end))
  }
  
  # ##############################################################################
  # # selbal
  # ##############################################################################
  # rownames(X) = paste("Sample", 1:nrow(X), sep = "_")
  # colnames(X) = paste("V", 1:ncol(X), sep = "")
  # 
  # # apply selbal, using CV to select the optimal number of variables
  # start.time = Sys.time()
  # selbal.fit = selbal.cv(x = X, y = as.vector(Y), n.fold = K)
  # end.time = Sys.time()
  # selbal.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")
  # 
  # # U (transformation) matrix
  # u.selbal = rep(0, p)
  # names(u.selbal) = colnames(X)
  # pba.pos = unlist(subset(
  #   selbal.fit$global.balance, subset = Group == "NUM", select = Taxa))
  # r = length(pba.pos)
  # pba.neg = unlist(subset(
  #   selbal.fit$global.balance, subset = Group == "DEN", select = Taxa))
  # s = length(pba.neg)
  # u.selbal[pba.pos] = 1 / r
  # u.selbal[pba.neg] = -1 / s
  # norm.const = sqrt((r * s) / (r + s))
  # u.selbal = norm.const * u.selbal
  # # check: these are equal
  # # lm(as.vector(Y) ~ log(X) %*% as.matrix(u.selbal))
  # # selbal.fit$glm
  # selbal.thetahat = coefficients(selbal.fit$glm)[2]
  # selbal.betahat = u.selbal %*% as.matrix(selbal.thetahat)
  # 
  # # evaluate model #
  # 
  # # 1. prediction error #
  # # 1a. on training set #
  # # get prediction error on training set
  # selbal.Yhat.train = predict.glm(selbal.fit$glm, 
  #                          newdata = data.frame(X), 
  #                          type = "response")
  # selbal.PE.train = crossprod(Y - selbal.Yhat.train) / n
  # # 1b. on test set #
  # # get prediction error on test set
  # rownames(X.test) = paste("Sample", 1:nrow(X), sep = "_")
  # colnames(X.test) = paste("V", 1:ncol(X), sep = "")
  # selbal.Yhat.test = predict.glm(selbal.fit$glm, newdata = data.frame(X.test), 
  #                         type = "response")
  # selbal.PE.test = crossprod(Y.test - selbal.Yhat.test) / n
  # 
  # # 2. estimation accuracy (i.t.o. beta) #
  # selbal.EA1 = sum(abs(selbal.betahat - beta))
  # selbal.EA2 = sqrt(crossprod(selbal.betahat - beta))
  # selbal.EAInfty = max(abs(selbal.betahat - beta))
  # 
  # # 3. selection accuracy (i.t.o. beta) #
  # selbal.non0.betahat = abs(selbal.betahat) > 10e-8
  # selbal.is0.betahat = selbal.betahat <= 10e-8
  # # FP
  # selbal.FP = sum(is0.beta & selbal.non0.betahat)
  # # FN
  # selbal.FN = sum((non0.beta != cl.non0.betahat) & non0.beta)
  # # TPR
  # selbal.TPR = sum((non0.beta == selbal.non0.betahat) & selbal.non0.betahat) / 
  #   sum(non0.beta)
  # 
  # saveRDS(c(
  #   "PEtr" = selbal.PE.train, 
  #   "PEte" = selbal.PE.test, 
  #   "EA1" = selbal.EA1, 
  #   "EA2" = selbal.EA2, 
  #   "EAInfty" = selbal.EAInfty, 
  #   "FP" = selbal.FP, 
  #   "FN" = selbal.FN, 
  #   "TPR" = selbal.TPR, 
  #   "timing" = selbal.timing,
  #   "betaSparsity" = bspars
  # ), 
  # file = paste0(output_dir, "/selbal_metrics", b, file.end))
  
  # roc
  # selbal.roc <- apply(slr$bet, 2, function(a) 
  #   roc.for.coef.LR(a, beta, ...))
}


