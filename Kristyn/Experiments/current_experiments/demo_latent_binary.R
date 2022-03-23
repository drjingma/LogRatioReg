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

library(pROC)

source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/HSClust.R")
source("Kristyn/Functions/slrnew.R")
source("Kristyn/Functions/codalasso.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/helper_functions.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

# expdecay Sigma #############################################################

# Settings to toggle with
sigma.settings = "latentVarModel_binary" # !!!
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_eps2 = 0.1
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0
b1 = 5 # 1, 0.5, 0.25
a0 = 0 # 0
theta.value = 1 # weight on a1: 1

##############################################################################
# generate data

# for(i in 1:500){
#   print(i)
# seed = i
seed = 465
set.seed(seed)
# get latent variable
# U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
U.all = matrix(runif(min = 0, max = 1, 2 * n), ncol = 1)
# simulate y from latent variable
y.all = rbinom(n = 2 * n, size = 1, p = as.vector(sigmoid(b0 + b1 * U.all)))
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
beta.true

# # plot the logistic regression curve
# x_seq = U.all
# y_seq = as.vector(sigmoid(b0 + b1 * U.all))
# plot(x_seq, y_seq)

##############################################################################
# new slr method (a balance regression method)
#   -- hierarchical spectral clustering (with rank 1 approximation)
##############################################################################
slrnew = slr(
  x = X, y = Y, num.clusters = 2, classification = TRUE, approx = FALSE)
slrnew_SBP = slrnew$sbp

slrnew.coefs = coefficients(slrnew$model)
slrnew.a0 = slrnew.coefs[1]
slrnew.thetahat = slrnew.coefs[-1]
slrnew.thetahat

slrnew.betahat.nonzero = getBetaFromTheta(slrnew.thetahat, sbp = slrnew_SBP)
slrnew.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew.betahat) = colnames(X)
slrnew.betahat[slrnew_activevars, ] = as.numeric(slrnew.betahat.nonzero)

# compute metrics on the selected model #
# prediction errors
# get prediction error on training set
slrnew.Yhat.train = predict.glm(
  slrnew$model, 
  newdata = data.frame(balance::balance.fromSBP(x = X, y = slrnew_SBP)), 
  type = "response")
slrnew.AUC.train = roc(
  Y, slrnew.Yhat.train, levels = c(0, 1), direction = "<")$auc
# get prediction error on test set
slrnew.Yhat.test = predict.glm(
  slrnew$model, newdata = data.frame(X.test), type = "response")
slrnew.AUC.test = roc(
  Y.test, slrnew.Yhat.test, levels = c(0, 1), direction = "<")$auc
# beta estimation accuracy, selection accuracy #
slrnew.metrics = getMetricsBalanceReg(
  thetahat = slrnew.thetahat, betahat = slrnew.betahat,
  true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
  true.beta = beta.true, metrics = c("betaestimation", "selection"))
slrnew.metrics = c(
  AUCtr = slrnew.AUC.train, AUCte = slrnew.AUC.test, slrnew.metrics)

slrnew2 = slrmult(x = X, y = Y, num.clusters = 2, approx = FALSE, classification = TRUE)
# x = X
# y = Y
# num.clusters = 2
# approx = FALSE
# classification = FALSE
slrnew2_activevars = names(slrnew2$index)
slrnew2_SBP = matrix(slrnew2$index)
rownames(slrnew2_SBP) = slrnew2_activevars

slrnew2.coefs = coefficients(slrnew2$model)
slrnew2.a0 = slrnew2.coefs[1]
slrnew2.thetahat = slrnew2.coefs[-1]

slrnew2.betahat.nonzero = getBetaFromTheta(slrnew2.thetahat, sbp = slrnew2_SBP)
slrnew2.betahat = matrix(0, nrow = ncol(X), ncol = 1)
rownames(slrnew2.betahat) = colnames(X)
slrnew2.betahat[slrnew2_activevars, ] = as.numeric(slrnew2.betahat.nonzero)

# compare
all.equal(as.numeric(slrnew.betahat), as.numeric(slrnew2.betahat))
