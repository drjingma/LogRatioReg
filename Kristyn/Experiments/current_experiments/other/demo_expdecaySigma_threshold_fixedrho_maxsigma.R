# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 12/13/2021

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
rho = 0.2
sigma_eps_seq = seq(0, 4.5, by = 0.1)

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:length(sigma_eps_seq)
) %dorng% {
  rho = 0.2 ################################# adjust according to settings above
  sigma_eps_seq = seq(0, 4.5, by = 0.1) ##### adjust according to settings above
  
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
  # rho = 0.2 # 0.2, 0.5
  scaling = TRUE
  
  # Population parameters
  sigma_eps =sigma_eps_seq[b]
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
  # leaf_types = rep("insignif cov", nrow(SBP))
  # leaf_types[non0.beta] = "signif cov"
  # balance_types = rep("insignif bal", ncol(SBP))
  # balance_types[theta[, 1] != 0] = "signif bal"
  # nodes_types = data.frame(
  #   name = c(colnames(SBP), rownames(SBP)),
  #   type = c(balance_types, leaf_types)
  # )
  # plotSBP(SBP, title = "Sigma", nodes_types = nodes_types) 
  
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
  print("checkpoint 1: SLR matrix ########################################################")
  # fields::image.plot(slrMat)
  slr.hsclust = HSClust(
    W = slrMat, 
    force_levelMax = TRUE, method = "kmeans")
  print("checkpoint 2: SLR hsclust ########################################################")
  slr.SBP = sbp.fromHSClust(
    levels_matrix = slr.hsclust$allLevels, row_names = names(beta))
  print("checkpoint 3: SLR sbp matrix ########################################################")
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
  print("checkpoint 4: SLR threshold ########################################################")
  
  plt = plotSBP(
    slr$sbp_thresh[[slr$min.idx[2]]], 
    title = paste0("sigma_eps = ", sigma_eps))
  return(plt)
}

saveRDS(
  res, 
  file = paste0(
    output_dir, 
    "/diagSigma_threshold_fixedrho_maxsigma",
    ".rds"
  ))
