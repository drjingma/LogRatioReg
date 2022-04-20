# date: 4/19/2022
# purpose: apply spectral clustering to population Gamma to see if we can 
#   recover the active set
rm(list=ls())

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/util.R")

library(ggplot2)

library(pheatmap)

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
sigma_eps1 = 0.1
sigma_eps2 = 0.1
# SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.25 # 1, 0.5, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0

##############################################################################
# generate data
set.seed(123)

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
# apply spectral clustering to population Gamma

popGammajk1 = function(
  alpha1j, alpha1k, beta1, var_epsilon, var_epsilonj, var_epsilonk, U){
  varU = stats::var(U)
  corrjk = ((alpha1j - alpha1k) * beta1 * varU) / 
    ((beta1^2 * varU + var_epsilon) * 
       ((alpha1j - alpha1k)^2 * varU + var_epsilonj + var_epsilonk))
  return(abs(corrjk))
}
popGammajk2 = function(
  alpha1j, alpha1k, beta1, var_epsilon, var_epsilonj, var_epsilonk){
  varU = (1 / 12) * (0.5 - (-0.5))
  corrjk = ((alpha1j - alpha1k) * beta1 * varU) / 
    ((beta1^2 * varU + var_epsilon) * 
       ((alpha1j - alpha1k)^2 * varU + var_epsilonj + var_epsilonk))
  return(abs(corrjk))
}
popGamma = function(
  alpha1, beta1, var_epsilon, var_epsilon2, U = NULL, 
  type = 1
){
  p = length(alpha1)
  if(length(var_epsilon2) == 1) var_epsilon2 = rep(var_epsilon2, p)
  
  rhoMat = matrix(0, p, p)
  for (j in 1:p){
    for (k in 1:p){
      if (k==j){next}
      else {
        if(type == 1){
          if(is.null(U)){
            stop("popGamma(): need to provide U argument when type == 1")
          }
          rhoMat[j, k] = popGammajk1(
            alpha1j = alpha1[j], alpha1k = alpha1[k], beta1 = beta1, 
            var_epsilon = var_epsilon, var_epsilonj = var_epsilon2[j], 
            var_epsilonk = var_epsilon2[k], U = U)
        } else if(type == 2){
          rhoMat[j, k] = popGammajk2(
            alpha1j = alpha1[j], alpha1k = alpha1[k], beta1 = beta1, 
            var_epsilon = var_epsilon, var_epsilonj = var_epsilon2[j], 
            var_epsilonk = var_epsilon2[k])
        } else{
          stop("popGamma(): invalid type argument!!")
        }
      }
    }
  }
  return(rhoMat)
}

popGamma1 = popGamma(
  alpha1 = c(a1, 0), beta1 = b1, var_epsilon = sigma_eps1^2, 
  var_epsilon2 = sigma_eps2^2, U = U.all[1:n, ], type = 1)
pheatmap(popGamma1)
popGamma2 = popGamma(
  alpha1 = c(a1, 0), beta1 = b1, var_epsilon = sigma_eps1^2, 
  var_epsilon2 = sigma_eps2^2, type = 2)
pheatmap(popGamma2)
spectral.clustering2(
  W = popGamma1, n_eig = 2, amini.regularization = FALSE, 
  highdegree.regularization = FALSE, include.leading.eigenvector = FALSE)


