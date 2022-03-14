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
source("Kristyn/Functions/supervisedlogratioseta.R")
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
sigma_eps1 = 0.01
sigma_eps2 = 0.01
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
# seed = 5
# set.seed(seed)
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
# Dr. Ma's new slr
##############################################################################

# apply supervised log-ratios, using CV to select threshold
slrnew = slr(x = X, y = Y, rank1approx = FALSE)
slrnew_activevars = names(slrnew$index)
slrnew_SBP = matrix(slrnew$index)
rownames(slrnew_SBP) = slrnew_activevars

slrnew.coefs = coefficients(slrnew$model)
slrnew.a0 = slrnew.coefs[1]
slrnew.thetahat = slrnew.coefs[-1]

slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew.betahat) = colnames(X)
slrnew.betahat[slrnew_activevars, ] = as.numeric(slrnew.betahat.nonzero)

beta.true
as.numeric(slrnew.betahat)
# slrnew$cors

# ##############################################################################
# # compositional lasso (a linear log contrast method)
# ##############################################################################
# classo = cv.func(
#   method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
#   nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
# 
# cl.lam.min.idx = which.min(classo$cvm)
# cl.a0 = classo$int[cl.lam.min.idx]
# cl.betahat = classo$bet[, cl.lam.min.idx]
# round(cl.betahat, 5)
