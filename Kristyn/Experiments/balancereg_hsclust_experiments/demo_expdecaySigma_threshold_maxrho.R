# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 12/19/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_hsclust_experiments/outputs"

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
# numSims = 100

# > max(abs(beta)) / 2 # 4.472136
rho_seq = seq(0, 1, by = 0.1)
rho_seq = rho_seq[-length(rho_seq)]

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:length(rho_seq)
) %dorng% {
  rho_seq = seq(0, 1, by = 0.1) ############# adjust according to settings above
  rho_seq = rho_seq[-length(rho_seq)]
  
  library(mvtnorm)
  
  library(Matrix)
  library(glmnet)
  
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  source("Kristyn/Functions/HSClust.R")
  
  # helper functions
  source("Kristyn/Functions/metrics.R")
  source("Kristyn/Functions/simulatedata.R")
  
  # for plots
  library(ggraph) # make dendrogram
  library(igraph) # transform dataframe to graph object: graph_from_data_frame()
  library(tidygraph)
  
  # expdecay Sigma #############################################################
  
  # Settings to toggle with
  sigma.settings = "expdecaySigma"
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  theta.settings = "pminus4"
  values.theta = 10
  linkage = "average"
  tol = 1e-4
  nlam = 100
  neta = 50
  intercept = TRUE
  K = 10
  n = 100
  p = 30
  rho = rho_seq[b] 
  scaling = TRUE
  
  # Population parameters
  sigma_eps = 0 ###################
  SigmaW = rgExpDecay(p, rho)$Sigma
  # fields::image.plot(SigmaW)
  
  # theta settings
  SigmaW_hsclust = HSClust(
    W = getSimilarityMatrix(unnormalized_similarity_matrix = SigmaW),
    levelMax = p, force_levelMax = TRUE)
  SBP = sbp.fromHSClust(levels_matrix = SigmaW_hsclust$allLevels)
  
  # a preliminary plot of the tree given by covariance matrix SigmaW
  # nodes_types = data.frame(
  #   name = c(colnames(SBP), rownames(SBP)),
  #   type = c(rep("balance", ncol(SBP)), rep("covariate", nrow(SBP)))
  # )
  # plotSBP(SBP, title = "Sigma", nodes_types = nodes_types) 
  
  # for each column (contrast), find which variables are included (1 or -1)
  indices.theta = p - 4
  
  # print(indices.theta)
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
  beta = as.vector(getBeta(theta, sbp = SBP))
  names(beta) <- paste0('s', 1:p)
  non0.beta = (beta != 0)
  is0.beta = abs(beta) <= 10e-8
  bspars = sum(non0.beta)
  
  # Population parameters, continued
  muW = c(rep(log(p), 5), rep(0, p - 5))
  names(muW) = names(beta)
  
  # plot the tree given by covariance matrix SigmaW, indicating 
  #   significant covariates and balances (theta's)
  leaf_types = rep("insignif cov", nrow(SBP))
  leaf_types[non0.beta] = "signif cov"
  balance_types = rep("insignif bal", ncol(SBP))
  balance_types[theta[, 1] != 0] = "signif bal"
  nodes_types = data.frame(
    name = c(colnames(SBP), rownames(SBP)),
    type = c(balance_types, leaf_types)
  )
  plt1 = plotSBP(
    SBP, title = "Sigma", nodes_types = nodes_types) # ...
  
  ##############################################################################
  # generate data
  fake.data = simulateBalanceReg(
    mu = muW, Sigma = SigmaW, U = getU(sbp = SBP), n = 2 * n, theta = theta, 
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
  
  # apply supervised log-ratios, using CV to select threshold and also lambda
  slrMat = getSupervisedMatrix(y = Y, X = X, rho.type = rho.type)
  # fields::image.plot(slrMat)
  slr.hsclust = HSClust(
    W = slrMat, 
    force_levelMax = TRUE, method = "kmeans")
  slr.SBP = sbp.fromHSClust(
    levels_matrix = slr.hsclust$allLevels, row_names = names(beta))
  slr = cvILReta(
    y = Y, X = X, 
    W = slrMat, # normalized similarity matrix (all values between 0 & 1)
    hsc_method = "kmeans", # "shimalik", "kmeans"
    force_levelMax = TRUE, 
    sbp = slr.SBP,
    lambda = NULL, nlam = nlam, 
    eta = NULL, neta = neta,
    nfolds = K, foldid = NULL, 
    intercept = intercept, 
    standardize = scaling
  )
  
  slr_sbp = slr$sbp_thresh[[slr$min.idx[2]]]
  # plt2 = plotSBP(
  #   slr_sbp, 
  #   title = paste0("rho = ", rho, ", sigma_eps = ", sigma_eps))
  slr_thetahat = slr$theta[[slr$min.idx[2]]][, slr$min.idx[[1]]]
  slr_betahat = getBeta(theta = slr_thetahat, sbp = slr_sbp)[, 1]
  slr_leaf_types = rep("not-selected cov", nrow(slr_sbp))
  slr_leaf_types[slr_betahat != 0] = "selected cov"
  slr_balance_types = rep("not-selected bal", ncol(slr_sbp))
  slr_balance_types[slr_thetahat != 0] = "selected bal"
  slr_nodes_types = data.frame(
    name = c(colnames(slr_sbp), rownames(slr_sbp)),
    type = c(slr_balance_types, slr_leaf_types)
  )
  plt2 = plotSBP(
    slr_sbp, title = "Sigma", nodes_types = slr_nodes_types) # ...
  
  plt1and2 = list(
    true_tree = plt1, 
    slr_tree = plt2
  )
  return(plt1and2)
}

saveRDS(
  res, 
  file = paste0(
    output_dir, 
    "/expdecaySigma_threshold_maxrho",
    ".rds"
  ))
