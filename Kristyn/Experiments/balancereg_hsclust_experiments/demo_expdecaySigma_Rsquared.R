# Purpose: observe different R^2 values for specified rho and sigma_eps values
# Date: 1/2/2021

################################################################################
# libraries and settings

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
values.theta = 1
linkage = "average"
tol = 1e-4
nlam = 100
neta = 50
intercept = TRUE
K = 10
n = 100
p = 30
scaling = TRUE
#################
# if rho = 0, 
#   sigma_eps = sqrt(2/3) => R^2 = 0.6
#   sigma_eps = sqrt(1/4) => R^2 = 0.8
# if rho = 0.2, 
#   sigma_eps = sqrt(0.7125333) => R^2 = 0.6
#   sigma_eps = sqrt(0.2672) => R^2 = 0.8
# if rho = 0.5, 
#   sigma_eps = sqrt(0.808333) => R^2 = 0.6
#   sigma_eps = sqrt(0.303125) => R^2 = 0.8
rho = 0.5
sigma_eps = sqrt(0.303125)
#################

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
seed = 123
# set.seed(seed)
# generate X
logW.all <- mvrnorm(n = 2 * n, mu = muW, Sigma = SigmaW) 
W.all <- exp(logW.all)
X.all <- sweep(W.all, 1, rowSums(W.all), FUN='/')
colnames(X.all) = names(beta)
# ilrX = computeBalances(X = X.all, U = getU(sbp = SBP))
# create the ilr(X.all) covariate by hand to
# generate y
SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
U.true.details = getUdetailed(sbp = SBP.true)
U.true = U.true.details$U
y.all = as.numeric(log(X.all) %*% U.true * values.theta) + rnorm(n) * sigma_eps

# subset out training and test sets
X = X.all[1:n, ]
X.test = X.all[-(1:n), ]
Y <- y.all[1:n]
Y.test <- y.all[-(1:n)]

# ##############################################################################
# # supervised log-ratios (a balance regression method)
# ##############################################################################
# 
# # apply supervised log-ratios, using CV to select threshold and also lambda
# slrMat = getSupervisedMatrix(y = Y, X = X, rho.type = rho.type)
# # fields::image.plot(slrMat)
# slr.hsclust = HSClust(
#   W = slrMat, 
#   force_levelMax = TRUE, method = "kmeans")
# slr.SBP = sbp.fromHSClust(
#   levels_matrix = slr.hsclust$allLevels, row_names = names(beta))
# slr = cvILReta(
#   y = Y, X = X, 
#   W = slrMat, # normalized similarity matrix (all values between 0 & 1)
#   hsc_method = "kmeans", # "shimalik", "kmeans"
#   force_levelMax = TRUE, 
#   sbp = slr.SBP,
#   lambda = NULL, nlam = nlam, 
#   eta = NULL, neta = neta,
#   nfolds = K, foldid = NULL, 
#   intercept = intercept, 
#   standardize = scaling
# )
# 
# slr_sbp = slr$sbp_thresh[[slr$min.idx[2]]]
# # plt2 = plotSBP(
# #   slr_sbp, 
# #   title = paste0("rho = ", rho, ", sigma_eps = ", sigma_eps))
# slr_thetahat = slr$theta[[slr$min.idx[2]]][, slr$min.idx[[1]]]
# slr_betahat = getBeta(theta = slr_thetahat, sbp = slr_sbp)[, 1]
# slr_leaf_types = rep("not-selected cov", nrow(slr_sbp))
# slr_leaf_types[slr_betahat != 0] = "selected cov"
# slr_balance_types = rep("not-selected bal", ncol(slr_sbp))
# slr_balance_types[slr_thetahat != 0] = "selected bal"
# slr_nodes_types = data.frame(
#   name = c(colnames(slr_sbp), rownames(slr_sbp)),
#   type = c(slr_balance_types, slr_leaf_types)
# )
# plt2 = plotSBP(
#   slr_sbp, title = "SLR tree", nodes_types = slr_nodes_types) # ...
# 
# plt1and2 = list(
#   true_tree = plt1, 
#   slr_tree = plt2
# )

# R-squared (in terms of the true model)
SSres = sum((y.all - as.numeric(log(X.all) %*% U.true * values.theta))^2)
SStot = sum((y.all - mean(y.all))^2)
Rsq = 1 - SSres/SStot
Rsq

