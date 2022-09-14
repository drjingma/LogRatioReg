# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 8/9/2022
rm(list=ls())

################################################################################
# libraries and settings

# Other simulation settings
numSims = 100

### in parallel loop ###

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("slr_analyses/Functions/slrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

library(ggplot2)

# Tuning parameters###########################################################

# Settings to toggle with
settings.name = "ContinuousResponseUnlabeled"
n = 100 # to make easier to threshold
p = 30
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_y = 0.1
sigma_x = 0.15
SBP.true = matrix(c(1, 1, 1, 1, 1, -1, -1, -1, -1, -1, rep(0, p - 10)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 3 # 0.5
c.value = 1 # Dr. Ma's formulation: c1 = c / k+
a0 = 0 # 0
n.unlabeled = n * 100
ulimit = 0.5

################################################################################
# Simulations #
################################################################################

b = 1
# for(b in 1:numSims){
set.seed(123 + b)
# rm(list=ls())

##############################################################################
# generate data
# get latent variable
U.all = matrix(runif(min = -ulimit, max = ulimit, 2 * n + n.unlabeled), ncol = 1)
# simulate y from latent variable
y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n + n.unlabeled) * sigma_y)
# simulate X: 
eps.all = matrix(
  rnorm((2 * n + n.unlabeled) * (p - 1)), 
  nrow = (2 * n + n.unlabeled)) * sigma_x
a1 = c.value * ilrtrans.true$ilr.trans.unscaled[-p] 
alrX.all = a0 + U.all %*% t(a1) + eps.all 
X.all <- alrinv(alrX.all)
colnames(X.all) = paste0('s', 1:p)

# subset out labeled and unlabeled sets
X.labeled = X.all[1:(2 * n), ]
Y.labeled = y.all[1:(2 * n)]
X2 = X.all[-(1:(2 * n)), ]
# subset out training and test sets from the labeled data
X = X.labeled[1:n, ]
X.test = X.labeled[-(1:n), ]
Y <- Y.labeled[1:n]
Y.test <- Y.labeled[-(1:n)]

# about beta
non0.beta = as.vector(SBP.true != 0)
bspars = sum(non0.beta)
# solve for beta
theta.value = c.value / ilrtrans.true$const
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
llc.coefs.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
  as.vector(ilrtrans.true$ilr.trans)

##############################################################################
# about the chosen settings

# Aitchison variation when j != k
# when j != k, aitchison var is Sjk = (c1 + c2)^2 Var[U] + 2 * sigma_eps2
varU = (2 * ulimit)^2 / 12
c1plusc2^2 * varU # term 1 # or want this to be small????
2 * sigma_x^2 # term 2 (want this term to dominate)
c1plusc2^2 * varU + 2 * sigma_x^2 # Sjk

# Correlation bt clr(Xj) & y
covclrXy = a1 * b1 * varU # covariance, in numerator
varclrX = a1^2 * varU + (1 - (1 / (p))) * sigma_x^2 # variance of clrX
vary = b1^2 * varU + sigma_y^2 # variance of y
# population correlations?
covclrXy / (sqrt(varclrX) * sqrt(vary))

# clustering
fields::image.plot(
  getAitchisonVar(X[, SBP.true != 0]))
spectral.clustering(getAitchisonVar(X[, SBP.true != 0]))

fields::image.plot(
  getAitchisonVar(rbind(X[, SBP.true != 0], X2[, SBP.true != 0])))
spectral.clustering(
  getAitchisonVar(rbind(X[, SBP.true != 0], X2[, SBP.true != 0]))
)

# how does theta affect Cor(clr(X)_j, y)?
func1 = function(theta){
  theta * abs(ilrtrans.true$ilr.trans[1, 1])
}
func1plusc2 = function(theta){
  theta * sum(abs(unique(ilrtrans.true$ilr.trans)))
}
funThetaCorXjY = function(theta, j = 1){
  a1j = (theta * ilrtrans.true$ilr.trans[j]) # j must be less than p
  varu = (2 * ulimit)^2 / 12
  covclrXjy = a1j * b1 * varu # covariance, in numerator
  varclrXj = a1j^2 * varu + (1 - (1 / (p))) * sigma_x^2 # variance of clrX
  vary = b1^2 * varu + sigma_y^2 # variance of y
  covclrXjy / (sqrt(varclrXj) * sqrt(vary))# population correlation
}
ggplot() + 
  geom_function(fun = func1, xlim = c(0, theta.value * 2), aes(color = "c1")) + 
  geom_function(fun = func1plusc2, xlim = c(0, theta.value * 2), aes(color = "c1+c2")) + 
  geom_function(fun = funThetaCorXjY, xlim = c(0, theta.value * 2), aes(color = "Cor(clr(Xj), y)")) + 
  geom_point(aes(x = theta.value, y = func1(theta.value), color = "c1")) +
  geom_point(aes(x = theta.value, y = func1plusc2(theta.value), color = "c1+c2")) +
  geom_point(aes(x = theta.value, y = funThetaCorXjY(theta.value), color = "Cor(clr(Xj), y)")) +
  labs(x = "x = value of theta in balance regression model", y = "")

# how does beta1 affect Cor(clr(X)_j, y)?
funBeta1CorXjY = function(beta1, j = 1){
  a1j = (theta.value * ilrtrans.true$ilr.trans[j]) # j must be less than p
  varu = (2 * ulimit)^2 / 12
  covclrXjy = a1j * beta1 * varu # covariance, in numerator
  varclrXj = a1j^2 * varu + (1 - (1 / (p))) * sigma_x^2 # variance of clrX
  vary = beta1^2 * varu + sigma_y^2 # variance of y
  covclrXjy / (sqrt(varclrXj) * sqrt(vary))# population correlation
}
ggplot() + 
  geom_function(fun = funBeta1CorXjY, xlim = c(0, b1 * 2), aes(color = "Cor(clr(Xj), y)")) + 
  geom_point(aes(x = b1, y = funBeta1CorXjY(b1), color = "Cor(clr(Xj), y)")) +
  labs(x = "x = value of beta1 in balance regression model", y = "")

# how is thresholding on labeled data?
clrX.labeled <- t(apply(X.labeled,1,clr))
plot(cor(clrX.labeled,as.matrix(Y.labeled)))
# # after CV splitting?
set.seed(1)
cvfoldid = sample(1:nrow(X.labeled), n / K)
clrX.labeled.cvfold <- t(apply(X.labeled[cvfoldid, ],1,clr))
plot(cor(clrX.labeled.cvfold,as.matrix(Y.labeled[cvfoldid])))

