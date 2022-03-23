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

source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/HSClust.R")
source("Kristyn/Functions/slrnew.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/helper_functions.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

# expdecay Sigma #############################################################

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
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0
b1 = 1 # 1, 0.5, 0.25
theta.value = 1 # weight on a1: 1
a0 = 0 # 0

##############################################################################
# generate data
seed = 5
set.seed(seed)
# get latent variable
U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
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
# solve for beta
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
  as.vector(ilrtrans.true$ilr.trans)

##############################################################################
# slr
##############################################################################

############################
slrnew0 = slr_old(x = X, y = Y, rank1approx = FALSE)
slrnew0_activevars = names(slrnew0$index)
slrnew0_SBP = matrix(slrnew0$index)
rownames(slrnew0_SBP) = slrnew0_activevars

slrnew0.coefs = coefficients(slrnew0$model)
slrnew0.a0 = slrnew0.coefs[1]
slrnew0.thetahat = slrnew0.coefs[-1]

slrnew0.betahat.nonzero = getBetaFromTheta(slrnew0.thetahat, sbp = slrnew0_SBP)
slrnew0.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew0.betahat) = colnames(X)
slrnew0.betahat[slrnew0_activevars, ] = as.numeric(slrnew0.betahat.nonzero)

############################
slrnew1 = slr(x = X, y = Y, num.clusters = 2, approx = FALSE)
# x = X
# y = Y
# num.clusters = 2
# approx = FALSE
# classification = FALSE
slrnew1_activevars = names(slrnew1$index)
slrnew1_SBP = matrix(slrnew1$index)
rownames(slrnew1_SBP) = slrnew1_activevars

slrnew1.coefs = coefficients(slrnew1$model)
slrnew1.a0 = slrnew1.coefs[1]
slrnew1.thetahat = slrnew1.coefs[-1]

slrnew1.betahat.nonzero = getBetaFromTheta(slrnew1.thetahat, sbp = slrnew1_SBP)
slrnew1.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew1.betahat) = colnames(X)
slrnew1.betahat[slrnew1_activevars, ] = as.numeric(slrnew1.betahat.nonzero)

############################
slrnew2 = hslr(x = X, y = Y, max.levels = 1, approx = FALSE)
slrnew2fit = slrnew2$fits[[1]]

slrnew2_activevars = names(slrnew2fit$index)
slrnew2_SBP = matrix(slrnew2fit$index)
rownames(slrnew2_SBP) = slrnew2_activevars

slrnew2.coefs = coefficients(slrnew2fit$model)
slrnew2.a0 = slrnew2.coefs[1]
slrnew2.thetahat = slrnew2.coefs[-1]

slrnew2.betahat.nonzero = getBetaFromTheta(slrnew2.thetahat, sbp = slrnew2_SBP)
slrnew2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew2.betahat) = colnames(X)
slrnew2.betahat[slrnew2_activevars, ] = as.numeric(slrnew2.betahat.nonzero)

# compare
all.equal(as.numeric(slrnew0.betahat), as.numeric(slrnew1.betahat))
all.equal(as.numeric(slrnew0.betahat), as.numeric(slrnew2.betahat))









##############################################################################
# slr with max.clusters = 5
##############################################################################
slrmult0 = slrmult(x = X, y = Y, max.clusters = 5, approx = FALSE)

which_fit_slrmult0 = which.max(slrmult0$max.Rsqs)
hslr0fit = slrmult0$fits[[which_fit_slrmult0]]
# hslr0fit$sbp

##############################################################################
# slr with cv
##############################################################################
cvslr0 = cv.slr(x = X, y = Y, max.clusters = 5, approx = FALSE, nfolds = K)

cvslr0$nclusters_min_idx
cvslr0$nclusters_1se_idx
cvslr0fit = cvslr0$models[[cvslr0$nclusters_1se_idx]]
# cvslr0fit$sbp


##############################################################################
# hslr using cor to choose level
##############################################################################
hslr0 = hslr(x = X, y = Y, max.levels = 5, approx = FALSE)

# which_fit = which(slrnew5$Rsqs == max(slrnew5$Rsqs, na.rm = TRUE), arr.ind = TRUE)
which_fit_hslr0 = which.max(hslr0$max.Rsqs)
hslr0fit = hslr0$fits[[which_fit_hslr0]]

