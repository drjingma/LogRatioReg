# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 4/21/2022

rm(list=ls())
library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/util.R")

source("Kristyn/Functions/slrtesting.R")

# Tuning parameters###########################################################

# Settings to toggle with
sigma.settings = "latentVarModel_corX"
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
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.25 # 1, 0.5, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0
rho_alrXj = 0.2

##############################################################################
# generate data
set.seed(123)

# get latent variable
U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
# simulate y from latent variable
y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
# simulate X: 
epsj.all = mvrnorm(
  n = 2 * n, mu = rep(0, p - 1), 
  Sigma = sigma_eps2 * rgExpDecay(p - 1, rho_alrXj)$Sigma)
# epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
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
bspars = sum(non0.beta)
# solve for beta
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
  as.vector(ilrtrans.true$ilr.trans)

##############################################################################
# slr method using population Gamma (subtractFrom1 = FALSE)
#   rank 1 approximation -- FALSE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
slr0 = slr_testing(
  subtractFrom1 = FALSE, 
  x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
  highdegree.regularization = FALSE)

slr0.coefs = getCoefsBM(
  coefs = coefficients(slr0$model), sbp = slr0$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = TRUE)
#   rank 1 approximation -- FALSE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
slr1 = slr_testing(
  subtractFrom1 = TRUE, 
  x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
  highdegree.regularization = FALSE)

slr1.coefs = getCoefsBM(
  coefs = coefficients(slr1$model), sbp = slr1$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = FALSE)
#   rank 1 approximation -- FALSE
#   amini regularization -- TRUE
#   high degree regularization -- FALSE
##############################################################################
slr0am = slr_testing(
  subtractFrom1 = FALSE, 
  x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
  highdegree.regularization = FALSE)

slr0am.coefs = getCoefsBM(
  coefs = coefficients(slr0am$model), sbp = slr0am$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = TRUE)
#   rank 1 approximation -- FALSE
#   amini regularization -- TRUE
#   high degree regularization -- FALSE
##############################################################################
slr1am = slr_testing(
  subtractFrom1 = TRUE, 
  x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
  highdegree.regularization = FALSE)

slr1am.coefs = getCoefsBM(
  coefs = coefficients(slr1am$model), sbp = slr1am$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = FALSE)
#   rank 1 approximation -- TRUE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
slr0ap = slr_testing(
  subtractFrom1 = FALSE, 
  x = X, y = Y, approx = TRUE, amini.regularization = FALSE, 
  highdegree.regularization = FALSE)

slr0ap.coefs = getCoefsBM(
  coefs = coefficients(slr0ap$model), sbp = slr0ap$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = TRUE)
#   rank 1 approximation -- TRUE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
slr1ap = slr_testing(
  subtractFrom1 = TRUE, 
  x = X, y = Y, approx = TRUE, amini.regularization = FALSE, 
  highdegree.regularization = FALSE)

slr1ap.coefs = getCoefsBM(
  coefs = coefficients(slr1ap$model), sbp = slr1ap$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = FALSE)
#   rank 1 approximation -- TRUE
#   amini regularization -- TRUE
#   high degree regularization -- FALSE
##############################################################################
slr0amap = slr_testing(
  subtractFrom1 = FALSE, 
  x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
  highdegree.regularization = FALSE)

slr0amap.coefs = getCoefsBM(
  coefs = coefficients(slr0amap$model), sbp = slr0amap$sbp)

##############################################################################
# slr method using population Gamma (subtractFrom1 = TRUE)
#   rank 1 approximation -- TRUE
#   amini regularization -- TRUE
#   high degree regularization -- FALSE
##############################################################################
slr1amap = slr_testing(
  subtractFrom1 = TRUE, 
  x = X, y = Y, approx = TRUE, amini.regularization = TRUE, 
  highdegree.regularization = FALSE)

slr1amap.coefs = getCoefsBM(
  coefs = coefficients(slr1amap$model), sbp = slr1amap$sbp)













