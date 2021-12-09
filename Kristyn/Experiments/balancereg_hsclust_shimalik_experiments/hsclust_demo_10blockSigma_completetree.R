# Purpose: demonstrate hierarchical spectral clustering 
# Date: 12/02/2021

################################################################################
# libraries and settings

library(limSolve)
library(mvtnorm)
library(Matrix)
library(glmnet)
library(compositions)
library(stats)

library(balance) # for sbp.fromHclust()
library(selbal)
library(propr)
library(sClust)
source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/HSClust.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/simulatedata.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidyverse)
require(tidygraph)

# Sigma with 10 blocks #########################################################

# Settings to toggle with
sigma.settings = "10blockSigma"
rho.type = "square" # 1 = "absolute value", 2 = "square"
theta.settings = "1blockpair4halves"  
# "pairperblock" => choose j corresp. to one pair of covariates for each block
# "2blockpairs4halves" => 
#   2 contrasts corresponding to 2 blocks each (accounts for 4 blocks so far), 
#   4 contrasts, each corresponding to half (or approx. half) of the variables 
#     in 4 different blocks (accounts for 8 blocks so far), and 
#   the other two blocks with inactive variables (i.e. not in any of the 
#     selected contrasts).
values.theta = 10
linkage = "average"
tol = 1e-4
nlam = 100
intercept = TRUE
K = 10
n = 100
p = 30
rho = 0.2 # 0.2, 0.5
scaling = TRUE

# Population parameters
sigma_eps = 0.0 # 0, 0.01 0.1
num.blocks = 10
SigmaWblock = matrix(rho, p / num.blocks, p / num.blocks)
for(i in 1:nrow(SigmaWblock)) SigmaWblock[i, i] = 1
SigmaW = as.matrix(bdiag(
  SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, 
  SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock, SigmaWblock))
heatmap(SigmaW, Rowv = NA, Colv = NA)
# SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
# U = getU(btree = SigmaWtree) # transformation matrix
# plot(SigmaWtree)

# theta settings
SigmaW_hsclust = HSClust(
  W = getSimilarityMatrix(unnormalized_similarity_matrix = SigmaW), 
  levelMax = p - 1, force_levelMax = TRUE)
SBP = sbp.fromHSClust(levels_matrix = SigmaW_hsclust$allLevels)
U = getU(sbp = SBP)

# a preliminary plot of the tree given by covariance matrix SigmaW
nodes_types = data.frame(
  name = c(colnames(SBP), rownames(SBP)),
  type = c(rep("balance", ncol(SBP)), rep("covariate", nrow(SBP)))
) 
plotSBP(SBP, title = "Sigma", nodes_types = nodes_types) 

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
  contrast_lengths = sapply(contrast.vars, length)
  contrasts.blockpair = which(
    contrast_lengths == 2 * (p / num.blocks))
  indices.theta1 = unname(contrasts.blockpair[1])
  blockpair.vars = contrast.vars[[contrasts.blockpair[1]]]
  # get the 4 contrasts, each corresponding to half (or approx. half) of the 
  #   vars in 4 different blocks
  num.closest.to.half = contrast_lengths[
    which.min(contrast_lengths - (p / num.blocks / 2))]
  contrasts.halves = which(
    contrast_lengths == num.closest.to.half) # 0.5 * (p / num.blocks))
  # identify the contrasts corresponding to half-blocks that aren't in
  #   the already-selected two blocks and that aren't in the same blocks
  contrasts.10blocks = which(contrast_lengths == (p / num.blocks))
  contrasts.halves.selected = c()
  for(i in 1:length(contrasts.halves)){
    # current half and which of the 10 blocks it is in
    contrast.half.tmp = contrasts.halves[i]
    contrast.half.vars.tmp = contrast.vars[[contrast.half.tmp]]
    contrast.half.block.tmp = NA
    for(l in 1:length(contrasts.10blocks)){
      contrast.block.tmp = contrasts.10blocks[l]
      contrast.block.vars.tmp = contrast.vars[[contrast.block.tmp]]
      if(any(contrast.half.vars.tmp %in% contrast.block.vars.tmp)){
        contrast.half.block.tmp = contrasts.10blocks[l]
        break
      }
    }
    contrast.half.block.vars.tmp = contrast.vars[[contrast.half.block.tmp]]
    # first check if the current half has any variables in the block-pair's set
    if(!any(contrast.half.vars.tmp %in% blockpair.vars)){
      # then check if it is any of the prev. selected halves' blocks
      if(length(contrasts.halves.selected) > 0){
        is.in.block = rep(FALSE, length(contrasts.halves.selected))
        for(j in 1:length(contrasts.halves.selected)){
          contrast.half.sel.tmp = contrasts.halves.selected[j]
          contrast.half.sel.vars.tmp = contrast.vars[[contrast.half.sel.tmp]]
          contrast.half.sel.block.tmp = NA
          for(l in 1:length(contrasts.10blocks)){
            contrast.block.tmp = contrasts.10blocks[l]
            contrast.block.vars.tmp = contrast.vars[[contrast.block.tmp]]
            if(any(contrast.half.sel.vars.tmp %in% contrast.block.vars.tmp)){
              contrast.half.sel.block.tmp = contrasts.10blocks[l]
              break
            }
          }
          if(contrast.half.block.tmp == contrast.half.sel.block.tmp){
            is.in.block[j] = TRUE
          }
        }
        if(all(!is.in.block)){
          contrasts.halves.selected = c(
            contrasts.halves.selected, contrast.half.tmp)
          if(length(contrasts.halves.selected) >= 4) break
        }
      } else{
        contrasts.halves.selected = c(
          contrasts.halves.selected, contrast.half.tmp)
      }
    }
  }
  indices.theta2 = unname(contrasts.halves.selected[1:4])
  indices.theta = c(indices.theta1, indices.theta2)
} else{
  stop("invalid theta.settings")
}
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
beta = as.vector(getBeta(theta, U = U))
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
plotSBP(SBP, title = "Sigma", nodes_types = nodes_types) 
# ggsave(
#   filename = paste0(
#     "20211202_",
#     sigma.settings, "_noise", sigma_eps,
#     "_", theta.settings,
#     "_val", values.theta[1],
#     "_SigmaWtree3.pdf"),
#   plot = last_plot(),
#   width = 8, height = 5, units = c("in")
# )

##############################################################################
# simulate data
set.seed(123)
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

# apply supervised log-ratios, using CV to select lambda
start.time = Sys.time()
# slr = cvSLR(
#   y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
#   rho.type = rho.type, linkage = linkage, standardize = scaling)
slrMat = getSupervisedDistanceMatrix(y = Y, X = X, rho.type = rho.type)
slr.hsclust = HSClust(
  W = getSimilarityMatrix(unnormalized_distance_matrix = slrMat), 
  levelMax = p - 1, force_levelMax = TRUE)
slr.SBP = sbp.fromHSClust(
  levels_matrix = slr.hsclust$allLevels, row_names = names(beta))
slr = cvILR(y = Y, X = X, sbp = slr.SBP, nlam = nlam, 
            nfolds = K, intercept = intercept, standardize = scaling)
end.time = Sys.time()

# timing metric
slr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

# choose lambda using cross-validated mse ##################################
slr.lam.min.idx = which.min(slr$cvm)
slr.a0 = slr$int[slr.lam.min.idx]
slr.thetahat = slr$bet[, slr.lam.min.idx]
slr.Uhat = getU(sbp = slr$sbp)
slr.betahat = getBeta(slr.thetahat, U = slr.Uhat)

# compute metrics on the selected model #
slr.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test, 
  ilrX.train = computeBalances(X, sbp = slr$sbp), 
  ilrX.test = computeBalances(X.test, sbp = slr$sbp), 
  n.train = n, n.test = n, 
  thetahat0 = slr.a0, thetahat = slr.thetahat, betahat = slr.betahat, 
  sbp = slr$sbp, 
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

slr.metrics = c(
  slr.metrics, 
  "timing" = slr.timing,
  "betaSparsity" = bspars
)

# plot the tree given by slr, this time also coloring the selected variables
slr.is0.betahat = abs(slr.betahat[, 1]) <= 1e-8
slr.non0.betahat = abs(slr.betahat[, 1]) > 1e-8
slr.is0.thetahat = slr.thetahat == 0
slr.non0.thetahat = !slr.is0.thetahat
leaf_types = rep(NA, nrow(slr.SBP))
leaf_types[non0.beta & slr.non0.betahat] = "selected, signif cov"
leaf_types[non0.beta & slr.is0.betahat] = "not-selected, signif cov"
leaf_types[is0.beta & slr.non0.betahat] = "selected, insignif cov"
leaf_types[is0.beta & slr.is0.betahat] = "not-selected, insignif cov"
balance_types = rep(NA, ncol(slr.SBP))
balance_types[slr.non0.thetahat] = "selected bal"
balance_types[slr.is0.thetahat] = "not-selected bal"
slr.nodestypes = data.frame(
  name = c(colnames(slr.SBP), rownames(slr.SBP)),
  type = c(balance_types, leaf_types)
)
plotSBP(slr.SBP, title = "supervised log-ratios", nodes_types = slr.nodestypes) 
# ggsave(
#   filename = paste0(
#     "20211202_",
#     sigma.settings, "_noise", sigma_eps,
#     "_", theta.settings,
#     "_val", values.theta[1],
#     "_slrtree2.pdf"),
#   plot = last_plot(),
#   width = 8, height = 5, units = c("in")
# )

# roc curves #################################################################
slr.roc <- apply(slr$bet, 2, function(a) 
  roc.for.coef.LR(a, beta, slr$sbp))

##############################################################################
# compositional lasso (a linear log contrast method)
##############################################################################

# fit model ##################################################################

# apply compositional lasso, using CV to select lambda
start.time = Sys.time()
classo = cv.func(
  method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam, 
  nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
end.time = Sys.time()

# timing metric
cl.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

# select tuning parameter and calculate metrics ##############################
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

cl.metrics = c(
  cl.metrics, 
  "timing" = cl.timing,
  "betaSparsity" = bspars
)

# roc curves #################################################################
cl.roc <- apply(classo$bet, 2, function(a) 
  roc.for.coef(a, beta))

##############################################################################
# oracle method (a balance regression method)
##############################################################################

# fit model ##################################################################

# apply oracle method, using CV to select lambda
start.time = Sys.time()
or.hsclust = HSClust(
  W = getSimilarityMatrix(unnormalized_similarity_matrix = SigmaW), 
  levelMax = p - 1, force_levelMax = TRUE)
or.SBP = sbp.fromHSClust(
  levels_matrix = or.hsclust$allLevels, row_names = names(beta))
oracle = cvILR(y = Y, X = X, sbp = or.SBP, nlam = nlam, 
               nfolds = K, intercept = intercept, standardize = scaling)
end.time = Sys.time()

# select tuning parameter and calculate metrics ##############################
or.lam.min.idx = which.min(oracle$cvm)
or.a0 = oracle$int[or.lam.min.idx]
or.thetahat = oracle$bet[, or.lam.min.idx]
or.Uhat = getU(sbp = oracle$sbp)
or.betahat = getBeta(or.thetahat, U = or.Uhat)

# compute metrics on the selected model #
or.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test, 
  ilrX.train = computeBalances(X, sbp = oracle$sbp), 
  ilrX.test = computeBalances(X.test, sbp = oracle$sbp), 
  n.train = n, n.test = n, 
  thetahat0 = or.a0, thetahat = or.thetahat, betahat = or.betahat, 
  sbp = oracle$sbp, 
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

or.metrics = c(
  or.metrics, 
  "timing" = or.timing,
  "betaSparsity" = bspars
)

# plot the tree given by oracle, this time also coloring the selected variables
or.is0.betahat = abs(or.betahat[, 1]) <= 1e-8
or.non0.betahat = abs(or.betahat[, 1]) > 1e-8
or.is0.thetahat = or.thetahat == 0
or.non0.thetahat = !or.is0.thetahat
leaf_types = rep(NA, nrow(or.SBP))
leaf_types[non0.beta & or.non0.betahat] = "selected, signif cov"
leaf_types[non0.beta & or.is0.betahat] = "not-selected, signif cov"
leaf_types[is0.beta & or.non0.betahat] = "selected, insignif cov"
leaf_types[is0.beta & or.is0.betahat] = "not-selected, insignif cov"
balance_types = rep(NA, ncol(or.SBP))
balance_types[or.non0.thetahat] = "selected bal"
balance_types[or.is0.thetahat] = "not-selected bal"
or.nodestypes = data.frame(
  name = c(colnames(or.SBP), rownames(or.SBP)),
  type = c(balance_types, leaf_types)
)
plotSBP(or.SBP, title = "oracle", nodes_types = or.nodestypes) 
# ggsave(
#   filename = paste0(
#     "20211202_",
#     sigma.settings, "_noise", sigma_eps,
#     "_", theta.settings,
#     "_val", values.theta[1],
#     "_oracletree2.pdf"),
#   plot = last_plot(),
#   width = 8, height = 5, units = c("in")
# )

# roc curves #################################################################
or.roc <- apply(oracle$bet, 2, function(a) 
  roc.for.coef.LR(a, beta, oracle$sbp))

##############################################################################
# propr method (a balance regression method)
##############################################################################

# fit model ##################################################################

# apply propr method, using CV to select lambda
start.time = Sys.time()
pr_res <- suppressMessages(propr(X, metric = "phs"))
pr.hsclust = HSClust(
  W = getSimilarityMatrix(unnormalized_distance_matrix = pr_res@matrix), 
  levelMax = p - 1, force_levelMax = TRUE)
pr.SBP = sbp.fromHSClust(
  levels_matrix = pr.hsclust$allLevels, row_names = names(beta))
pr = cvILR(y = Y, X = X, sbp = pr.SBP, nlam = nlam, 
           nfolds = K, intercept = intercept, standardize = scaling)
end.time = Sys.time()

# timing metric
pr.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

# select tuning parameter and calculate metrics ##############################
pr.lam.min.idx = which.min(pr$cvm)
pr.a0 = pr$int[pr.lam.min.idx]
pr.thetahat = pr$bet[, pr.lam.min.idx]
pr.Uhat = getU(sbp = pr$sbp)
pr.betahat = getBeta(pr.thetahat, U = pr.Uhat)

# compute metrics on the selected model #
pr.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test, 
  ilrX.train = computeBalances(X, sbp = pr$sbp), 
  ilrX.test = computeBalances(X.test, sbp = pr$sbp), 
  n.train = n, n.test = n, 
  thetahat0 = pr.a0, thetahat = pr.thetahat, betahat = pr.betahat, 
  sbp = pr.SBP, 
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

pr.metrics = c(
  pr.metrics, 
  "timing" = pr.timing,
  "betaSparsity" = bspars
)

# plot the tree given by oracle, this time also coloring the selected variables
pr.is0.betahat = abs(pr.betahat[, 1]) <= 1e-8
pr.non0.betahat = abs(pr.betahat[, 1]) > 1e-8
pr.is0.thetahat = pr.thetahat == 0
pr.non0.thetahat = !pr.is0.thetahat
leaf_types = rep(NA, nrow(pr.SBP))
leaf_types[non0.beta & pr.non0.betahat] = "selected, signif cov"
leaf_types[non0.beta & pr.is0.betahat] = "not-selected, signif cov"
leaf_types[is0.beta & pr.non0.betahat] = "selected, insignif cov"
leaf_types[is0.beta & pr.is0.betahat] = "not-selected, insignif cov"
balance_types = rep(NA, ncol(pr.SBP))
balance_types[pr.non0.thetahat] = "selected bal"
balance_types[pr.is0.thetahat] = "not-selected bal"
pr.nodestypes = data.frame(
  name = c(colnames(pr.SBP), rownames(pr.SBP)),
  type = c(balance_types, leaf_types)
)
plotSBP(pr.SBP, title = "propr", nodes_types = pr.nodestypes) 
# ggsave(
#   filename = paste0(
#     "20211202_",
#     sigma.settings, "_noise", sigma_eps,
#     "_", theta.settings,
#     "_val", values.theta[1],
#     "_proprtree2.pdf"),
#   plot = last_plot(),
#   width = 8, height = 5, units = c("in")
# )

# roc curves #################################################################
pr.roc <- apply(pr$bet, 2, function(a) 
  roc.for.coef.LR(a, beta, pr$sbp))

################################################################################
# compare beta-hats for slr and propr

betahats = data.frame(
  index = 1:p, beta = beta, slr = slr.betahat, propr = pr.betahat)
betahats.mlt = melt(betahats, id.vars = "index")
ggplot(
  betahats.mlt, aes(x = index, y = value, color = variable, shape = variable)) +
  geom_point(size = 3, alpha = 0.5)
# ggsave(
#   filename = paste0(
#     "20211202_",
#     sigma.settings, "_noise", sigma_eps,
#     "_", theta.settings,
#     "_val", values.theta[1],
#     "_betahats.pdf"),
#   plot = last_plot(),
#   width = 8, height = 5, units = c("in")
# )
