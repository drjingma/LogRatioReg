# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 8/7/2022
rm(list=ls())

################################################################################
# libraries and settings

set.seed(1)
library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("slr_analyses/Functions/slrs.R")
source("slr_analyses/Functions/semislrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

# Tuning parameters###########################################################

# Settings to toggle with
sigma.settings = "latentVarModel_miss"
n = 5000
p = 30
sigma_eps1 = 0.1 # zeta (for y)
sigma_eps2 = 0.1 # sigma_j (for x)
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix, i.e.
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0
b1 = 0.75 
theta.value = 0.5 # weight on a1
a0 = 0 
prop.missing = 0.5
ulimit = 0.5 # limits of U ~ Uniform(-ulimit, ulimit)

# generate data ################################################################

# get latent variable
U.all = matrix(runif(min = -ulimit, max = ulimit, 2 * n), ncol = 1)
# simulate y from latent variable
y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
# simulate X: 
epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
a1 = theta.value * ilrtrans.true$ilr.trans[-p] 
#   alpha1j = {
#     c1 = theta*ilr.const/k+   if j \in I+
#     -c2 = -theta*ilr.const/k-  if j \in I-
#     0                       o/w
#   }
alrXj.all = a0 + U.all %*% t(a1) + epsj.all #log(Xj/Xp) =alpha0j+alpha1j*U+epsj
X.all <- alrinv(alrXj.all)
colnames(X.all) = paste0('s', 1:p)

# solve for the coefficients in the linear log contrast model
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
linearlogcontrast.coefficients = (b1 / (ilrtrans.true$const * c1plusc2)) * 
  as.vector(ilrtrans.true$ilr.trans)

##############################################################################
##############################################################################
##############################################################################
# about the chosen settings

# b1 / c * (c1 + c2)
b1
ilrtrans.true$const * c1plusc2
b1 / (ilrtrans.true$const * c1plusc2)

# Aitchison variation when j != k ############################################
# when j != k, aitchison var is Sjk = (c1 + c2)^2 Var[U] + 2 * sigma_eps2
varU = (2 * ulimit)^2 / 12
c1plusc2^2 * varU # term 1
2 * sigma_eps2^2 # term 2 (want this term to dominate)

# Correlation bt clr(Xj) & y #################################################
covclrXy = a1 * b1 * varU # covariance, in numerator
varclrX = a1^2 * varU + (1 - (1 / (p))) * sigma_eps2^2 # variance of clrX
vary = b1^2 * varU + sigma_eps1^2 # variance of y
# population correlations?
covclrXy / (sqrt(varclrX) * sqrt(vary))

# data correlations:
round(apply(
  t(apply(
    X.all,1,
    function(a) log(a) - mean(log(a)))), 2, 
  function(xj) cor(xj, y.all)
), 3)

# variance of y
var(y.all)
vary
