# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Note: Here, we simulate from a 10-block-diagonal Sigma and intend to:
#   - have one active contrast pair within each block
#   - compare slr to coat, where the latter captures correlated structure
#       among covariates, but not whether the covariates are predicted
# Date: 10/19/2021

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
  source("Kristyn/Functions/supervisedlogratiosalpha.R")
  source("COAT-master/coat.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  source("Kristyn/Functions/simulatedata.R")
  
  # Settings to toggle with
  sigma.settings = "10blockSigmaExpDecay"
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "1blockpair4halves"  
  # "pairperblock" => choose j corresp. to one pair of covariates for each block
  # "2blockpairs4halves" => 
  #   2 contrasts corresponding to 2 blocks each (accounts for 4 blocks so far), 
  #   4 contrasts, each corresponding to half (or approx. half) of the variables 
  #     in 4 different blocks (accounts for 8 blocks so far), and 
  #   the other two blocks with inactive variables (i.e. not in any of the 
  #     selected contrasts).
  values.theta = 1
  linkage = "average"
  tol = 1e-4
  nlam = 100
  intercept = TRUE
  K = 10
  n = 100
  p = 200
  rho = 0.2
  scaling = TRUE
  
  # Population parameters
  sigma_eps = 0.1 # 0.1, 0.5
  num.blocks = 10
  SigmaWblock = rgExpDecay(p / num.blocks, rho)$Sigma
  SigmaW = as.matrix(bdiag(
    SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, 
    SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock))
  SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
  U = getU(btree = SigmaWtree) # transformation matrix
  # plot(SigmaWtree)
  # sbp.fromHclust(SigmaWtree)
  
  # theta settings
  SBP = sbp.fromHclust(SigmaWtree)
  # for each column (contrast), find which variables are included (1 or -1)
  contrast.vars = apply(SBP, 2, FUN = function(col) which(col != 0))
  if(theta.settings == "1blockpair4halves"){
    # "1blockpair4halves" => 
    #   1 contrast corresponding to 2 blocks (accounts for 2 blocks so far), 
    #   4 contrasts, each corresponding to half (or approx. half) of the vars 
    #     in 4 different blocks (accounts for 8 blocks so far), and 
    #   the other 4 blocks with inactive vars (i.e. not in any of the 
    #     selected contrasts).
    # get the 1 contrast corresponding to 2 blocks
    block.contrasts.1blockpair = which(
      sapply(contrast.vars, length) == 2 * (p / num.blocks))
    indices.theta1 = unname(block.contrasts.1blockpair)
    # get the 4 contrasts, each corresponding to half (or approx. half) of the 
    #   vars in 4 different blocks
    block.contrasts.halves = which(
      sapply(contrast.vars, length) == 9) # 0.5 * (p / num.blocks))
    indices.theta2 = unname(block.contrasts.halves[1:4])
    indices.theta = c(indices.theta1, indices.theta2)
  } else{
    stop("invalid theta.settings")
  }
  print(indices.theta)
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
  
  file.end0 = paste0(
    "_", sigma.settings,
    "_", theta.settings, 
    "_dim", n, "x", p, 
    "_noise", sigma_eps,
    "_rho", rho, 
    "_int", intercept,
    "_scale", scaling,
    "_sim")
  
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
  # get lamda sequences
  
  # import lambdas
  cl.lams.mat = matrix(NA, nrow = nlam, ncol = numSims)
  slr.lams.mat = matrix(NA, nrow = nlam, ncol = numSims)
  # selbal.lams.mat = matrix(NA, nrow = nlam, ncol = numSims)
  or.lams.mat = matrix(NA, nrow = nlam, ncol = numSims)
  # coat.lams.mat = matrix(NA, nrow = nlam, ncol = numSims)
  pr.lams.mat = matrix(NA, nrow = nlam, ncol = numSims)
  for(i in 1:numSims){
    # classo
    cl.lams.mat[, i] = readRDS(
      paste0(output_dir, "/models", "/classo_model", file.end0, i, ".rds"
      ))$lambda
    # slr
    slr.lams.mat[, i] = readRDS(
      paste0(output_dir, "/models", "/slr_model", file.end0, i, ".rds"
      ))$lambda
    # # selbal
    # selbal.lams.mat[, i] = readRDS(
    #   paste0(output_dir, "/models", "/selbal_model", file.end0, i, ".rds"
    #   ))$lambda
    # oracle
    or.lams.mat[, i] = readRDS(
      paste0(output_dir, "/models", "/oracle_model", file.end0, i, ".rds"
      ))$lambda
    # # coat
    # coat.lams.mat[, i] = readRDS(
    #   paste0(output_dir, "/models", "/coat_model", file.end0, i, ".rds"
    #   ))$lambda
    # propr
    pr.lams.mat[, i] = readRDS(
      paste0(output_dir, "/models", "/propr_model", file.end0, i, ".rds"
      ))$lambda
  }
  
  # lambda min and max values
  cl.bounds = log(c(min(cl.lams.mat), max(cl.lams.mat)))
  slr.bounds = log(c(min(slr.lams.mat), max(slr.lams.mat)))
  # selbal.bounds = log(c(min(has.selbal), max(has.selbal)))
  or.bounds = log(c(min(or.lams.mat), max(or.lams.mat)))
  # coat.bounds = log(c(min(coat.lams), max(coat.lams)))
  pr.bounds = log(c(min(pr.lams.mat), max(pr.lams.mat)))
  
  # new lambda sequences, one for each method
  cl.lambda.seq = exp(seq(max(cl.bounds), min(cl.bounds),length.out = nlam))
  slr.lambda.seq = exp(seq(max(slr.bounds), min(slr.bounds),length.out = nlam))
  # selbal.lambda.seq = exp(seq(max(selbal.bounds), min(selbal.bounds),length.out = nlam))
  or.lambda.seq = exp(seq(max(or.bounds), min(or.bounds),length.out = nlam))
  # coat.lambda.seq = exp(seq(max(coat.bounds), min(coat.bounds),length.out = nlam))
  pr.lambda.seq = exp(seq(max(pr.bounds), min(pr.bounds),length.out = nlam))
  
  ##############################################################################
  # supervised log-ratios
  ##############################################################################
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/models_samelam", "/slr_model", file.end))){
    slr.model.already.existed = FALSE
    # apply supervised log-ratios, using CV to select lambda
    slr = cvSLR(
      y = Y, X = X, lambda = slr.lambda.seq, nfolds = K, intercept = intercept, 
      rho.type = rho.type, linkage = linkage, standardize = scaling)
    saveRDS(
      slr, 
      paste0(
        output_dir, "/roc_samelam/models_samelam", "/slr_model", file.end))
  } else{
    slr.model.already.existed = TRUE
  }
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/roccurves_samelam", "/slr_roc", file.end)) |
    slr.model.already.existed == FALSE){
    
    # import model
    slr = readRDS(paste0(
      output_dir, "/roc_samelam/models_samelam", "/slr_model", file.end))
    slr.btree = slr$btree
    slr.SBP = sbp.fromHclust(slr.btree)
    
    # roc
    slr.roc <- apply(slr$bet, 2, function(a)
      roc.for.coef.LR(a, beta, slr.SBP))
    slr.roc = as.data.frame(t(slr.roc))
    slr.roc$lambda = slr.lambda.seq
    
    saveRDS(slr.roc, paste0(
      output_dir, "/roc_samelam/roccurves_samelam", "/slr_roc", file.end))
  }
  
  ##############################################################################
  # compositional lasso
  ##############################################################################
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/models_samelam", "/classo_model", file.end))){
    cl.model.already.existed = FALSE
    # apply compositional lasso, using CV to select lambda
    classo = cv.func(
      method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), lambda = cl.lambda.seq, 
      nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
    saveRDS(classo, paste0(
      output_dir, "/roc_samelam/models_samelam", "/classo_model", file.end))
  } else{
    cl.model.already.existed = TRUE
  }
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/roccurves_samelam", "/classo_roc", file.end)) |
    cl.model.already.existed == FALSE){
    
    # import model
    classo = readRDS(paste0(
      output_dir, "/roc_samelam/models_samelam", "/classo_model", file.end))
    
    # roc
    cl.roc <- apply(classo$bet, 2, function(a)
      roc.for.coef(a, beta))
    cl.roc = as.data.frame(t(cl.roc))
    cl.roc$lambda = cl.lambda.seq
    
    saveRDS(cl.roc, paste0(
      output_dir, "/roc_samelam/roccurves_samelam", "/classo_roc", file.end))
  }
  
  ##############################################################################
  # oracle method
  ##############################################################################
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/models_samelam", "/oracle_model", file.end))){
    or.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    oracle = cvILR(y = Y, X = X, btree = SigmaWtree, U = U, lambda = or.lambda.seq, 
                   nfolds = K, intercept = intercept, standardize = scaling)
    saveRDS(oracle, paste0(
      output_dir, "/roc_samelam/models_samelam", "/oracle_model", file.end))
  } else{
    or.model.already.existed = TRUE
  }
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/roccurves_samelam", "/oracle_roc", file.end)) |
    or.model.already.existed == FALSE){
    
    # import model
    oracle = readRDS(paste0(
      output_dir, "/roc_samelam/models_samelam", "/oracle_model", file.end))
    or.btree = SigmaWtree
    or.SBP = sbp.fromHclust(or.btree)
    rownames(or.SBP) = names(beta)
    
    # roc
    or.roc <- apply(oracle$bet, 2, function(a)
      roc.for.coef.LR(a, beta, or.SBP))
    or.roc = as.data.frame(t(or.roc))
    or.roc$lambda = or.lambda.seq
    
    saveRDS(or.roc, paste0(
      output_dir, "/roc_samelam/roccurves_samelam", "/oracle_roc", file.end))
  }
  
  ##############################################################################
  # propr method
  ##############################################################################
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/models_samelam", "/propr_model", file.end))){
    pr.model.already.existed = FALSE
    # apply oracle method, using CV to select lambda
    pr <- propr(X, metric = "phs")
    pr.tree = hclust(as.dist(pr@matrix),method = linkage)
    pr = cvILR(y = Y, X = X, btree = pr.tree, lambda = pr.lambda.seq, 
               nfolds = K, intercept = intercept, standardize = scaling)
    saveRDS(pr, paste0(
      output_dir, "/roc_samelam/models_samelam", "/propr_model", file.end))
  } else{
    pr.model.already.existed = TRUE
  }
  
  if(!file.exists(paste0(
    output_dir, "/roc_samelam/roccurves_samelam", "/propr_roc", file.end)) |
    pr.model.already.existed == FALSE){
    
    # import model
    pr = readRDS(paste0(
      output_dir, "/roc_samelam/models_samelam", "/propr_model", file.end))
    pr.btree = pr$btree
    pr.SBP = sbp.fromHclust(pr.btree)
    
    # roc
    pr.roc <- apply(pr$bet, 2, function(a)
      roc.for.coef.LR(a, beta, pr.SBP))
    pr.roc = as.data.frame(t(pr.roc))
    pr.roc$lambda = pr.lambda.seq
    
    saveRDS(pr.roc, paste0(
      output_dir, "/roc_samelam/roccurves_samelam", "/propr_roc", file.end))
  }
  
  # ##############################################################################
  # # supervised log-ratios alpha = 0.5
  # ##############################################################################
  # 
  # if(!file.exists(paste0(
  #   output_dir, "/roc_samelam/models_samelam", "/slralpha0.5_model", file.end))){
  #   slr0.5.model.already.existed = FALSE
  #   # apply supervised log-ratios, using CV to select lambda
  #   alpha = 0.5
  #   slr0.5 = cvSLRalpha(
  #     y = Y, X = X, lambda = slr0.5.lambda.seq, nfolds = K, alpha = alpha, 
  #     intercept = intercept, rho.type = rho.type, linkage = linkage, 
  #     scaling = scaling)
  #   saveRDS(slr0.5, paste0(
  #     output_dir, "/roc_samelam/models_samelam", "/slralpha0.5_model", file.end))
  # } else{
  #   slr0.5.model.already.existed = TRUE
  #   slr0.5 = readRDS(paste0(
  #     output_dir, "/roc_samelam/models_samelam", "/slralpha0.5_model", file.end))
  # }
  # 
  #   slr0.5.btree = slr0.5$btree
  #   slr0.5.SBP = sbp.fromHclust(slr0.5.btree)
  # 
  # if(!file.exists(paste0(
  #   output_dir, "/roc_samelam/roccurves_samelam", "/slralpha0.5_roc", file.end)) | 
  #    slr0.5.model.already.existed == FALSE){
  #   # roc
  #   slr0.5.roc <- apply(slr0.5$bet, 2, function(a) 
  #     roc.for.coef.LR(a, beta, slr0.5.SBP))
  #   
  #   saveRDS(slr0.5.roc, paste0(
  #     output_dir, "/roc_samelam/roccurves_samelam", "/slralpha0.5_roc", file.end))
  # }
}


