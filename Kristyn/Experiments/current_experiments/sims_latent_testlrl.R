# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 5/15/2022

rm(list=ls())

b = 1
set.seed(b)

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/util.R")

# Tuning parameters###########################################################

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
sigma_eps1 = 0.1 # sigma (for y)
sigma_eps2 = 0.1 # sigma_j (for x)
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.25 # 1, 0.5, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0

##############################################################################
# generate data
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
bspars = sum(non0.beta)
# solve for beta
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
  as.vector(ilrtrans.true$ilr.trans)

##############################################################################
# Log-Ratio Lasso
# -- regresses on pairwise log-ratios
##############################################################################
library(logratiolasso)
source("Kristyn/Functions/logratiolasso.R")
W = log(X)
Wc = scale(W, center = TRUE, scale = FALSE)
Yc = Y - mean(Y)
lrl_model <- cv_two_stage(z = Wc, y = Yc, n_folds = K)

# testing code
z = Wc
y = Yc
family = "gaussian"
lambda_1 = NULL
k_max = 10
n_folds = K
nlambda = 20
folds = NULL





