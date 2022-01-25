rm(list=ls())
# Purpose: try using Data/Q.txt phylogenetic tree to generate data from
# Date: 1/24/2022

################################################################################
# libraries and settings

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(propr)

source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/supervisedlogratioseta.R")
source("Kristyn/Functions/HSClust.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/simulatedata.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

################################################################################
# data
Q = read.table(file = "Data/Q.txt")
Q = as.matrix(Q)

# Sigma with exponential decay ###############################################

# Settings to toggle with
sigma.settings = "Qmatrix"
values.theta = 1
linkage = "average"
tol = 1e-4
nlam = 100
neta = 50
intercept = TRUE
K = 10
n = 100
p = ncol(Q)
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
get_sigma_eps = function(theta_val, Rsq_val, rho_val){
  sigma_eps_sq.tmp = theta_val^2 * (1 - Rsq_val) / Rsq_val + 
    theta_val^2 * (1 - Rsq_val) * (rho_val^3 + 2 * rho_val^2 + 3 * rho_val) / 
    (10 * Rsq_val)
  return(sqrt(sigma_eps_sq.tmp))
}
rho = 0.2 #
desired_Rsquared = 0.6 #
sigma_eps = get_sigma_eps(
  theta_val = values.theta, Rsq_val = desired_Rsquared, rho_val = rho)

# Population parameters

### apply hierarchical spectral clustering to similarity matrix Q?
# hsclust_Q = HSClust(
#   W = Q, # assuming Q is a similarity matrix
#   force_levelMax = TRUE, method = "kmeans")
# sbp_Q = sbp.fromHSClust(levels_matrix = hsclust_Q$allLevels)
### apply hierarchical clustering to similarity matrix Q?
hclust_Q = hclust(
  as.dist(getSimilarityMatrix(unnormalized_similarity_matrix = Q)), 
  method = linkage)
sbp_Q = sbp.fromHclust(hclust_Q)

# define covariate parameters
SigmaW = rgExpDecay(p, rho)$Sigma
muW = c(rep(log(p), 5), rep(0, p - 5))
names(muW) = paste0('s', 1:p)

# define theta (directly related to sparsity)
# randomly choose 5% of balances/internal nodes that at most 5 active variables
#   i.e. excluse balances/internal nodes with more than 5 active variables
contrast_vars = apply(sbp_Q, 2, FUN = function(col) which(col != 0))
contrast_lens = sapply(contrast_vars, length)
viable_contrasts = colnames(sbp_Q)[contrast_lens <= 5]
num_select_viable_contrasts = signif(length(viable_contrasts) * 0.05, 0)
# View(sbp_Q[, viable_contrasts])
theta = rep(0, p - 1)
names(theta) = colnames(sbp_Q)
selected_viable_contrasts = sample(
  x = viable_contrasts, size = num_select_viable_contrasts, 
  replace = FALSE)
theta[selected_viable_contrasts] = values.theta
beta = as.vector(getBeta(theta = theta, sbp = sbp_Q))

# define beta
names(beta) <- rownames(sbp_Q)
non0.beta = (beta != 0)
is0.beta = abs(beta) <= 10e-8
bspars = sum(non0.beta)

##############################################################################
# generate data

fake.data = simulateBalanceReg(
  mu = muW, Sigma = SigmaW, sbp = sbp_Q, n = 2 * n, theta = theta, 
  sigma.noise = sigma_eps)
colnames(fake.data$X) = names(beta)

# subset out training and test sets
X = fake.data$X[1:n, ]
X.test = fake.data$X[-(1:n), ]
Y <- fake.data$y[1:n, , drop = TRUE]
Y.test <- fake.data$y[-(1:n), , drop = TRUE]

##############################################################################
# supervised log-ratios (a balance regression method)
#   -- hierarchical clustering
##############################################################################
start.time = Sys.time()
# apply hierarchical clustering to the SLR distance matrix
slrDistMat = getSupervisedMatrix(
  y = Y, X = X, type = "distance")
slrhc_btree = hclust(as.dist(slrDistMat), method = linkage)
slrhc_SBP = sbp.fromHclust(slrhc_btree)
# apply supervised log-ratios, using CV to select threshold and also lambda
slrhc = cvILR(y = Y, X = X, sbp = slrhc_SBP, nlam = nlam,
              nfolds = K, intercept = intercept, standardize = scaling)
end.time = Sys.time()
slrhc.timing = difftime(time1 = end.time, time2 = start.time, units = "secs")

slrhc.lam.min.idx = which.min(slrhc$cvm)
slrhc.a0 = slrhc$int[slrhc.lam.min.idx]
slrhc.thetahat = slrhc$bet[, slrhc.lam.min.idx]
slrhc.betahat = getBeta(slrhc.thetahat, sbp = slrhc$sbp)

# compute metrics on the selected model #
slrhc.metrics = getMetricsBalanceReg(
  y.train = Y, y.test = Y.test,
  ilrX.train = computeBalances(X, sbp = slrhc$sbp),
  ilrX.test = computeBalances(X.test, sbp = slrhc$sbp),
  n.train = n, n.test = n,
  thetahat0 = slrhc.a0, thetahat = slrhc.thetahat, betahat = slrhc.betahat,
  sbp = slrhc$sbp,
  true.beta = beta, is0.true.beta = is0.beta, non0.true.beta = non0.beta)

# plot the tree given by slr-hc, indicating significant covariates
slrhc_leaf_types = rep("covariate", nrow(slrhc$sbp))
slrhc_balance_types = rep("balance", ncol(slrhc$sbp))
slrhc_nodes_types = data.frame(
  name = c(colnames(slrhc$sbp), rownames(slrhc$sbp)),
  type = c(slrhc_balance_types, slrhc_leaf_types)
)
plotSBP(slrhc$sbp, title = "slr-hc", nodes_types = slrhc_nodes_types)
fields::image.plot(slrDistMat)
