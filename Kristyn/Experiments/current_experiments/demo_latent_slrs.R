# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 4/27/2022

rm(list=ls())
library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/slr1sc.R")
source("Kristyn/Functions/slrtesting.R")
source("Kristyn/Functions/util.R")

source("Kristyn/Functions/slrtesting.R")
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
set.seed(1)

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
# slr method
#   similarity.matrix -- TRUE
#   maxGamma -- FALSE
#   spectral.clustering.algorithm == "kmeans"
#   rank 1 approximation -- FALSE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
slrk = slr_testing(
  similarity.matrix = TRUE, maxGamma = TRUE,
  x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
  highdegree.regularization = FALSE, spectral.clustering.method = "kmeans")

Gamma = slrmatrix(x = X, y = Y)
par(mfrow = c(1, 1))
fields::image.plot(max(Gamma) - Gamma)

slrk.coefs = getCoefsBM(
  coefs = coefficients(slrk$model), sbp = slrk$sbp)

par(mfrow = c(1,1))
fields::image.plot(slrk$kernel)
pheatmap(slrk$spectralclustering1$L$W)
pheatmap(slrk$spectralclustering1$L$W.tmp)
par(mfrow = c(1, 3))
# plot(slrk$spectralclustering1$ei$values, ylab = "eigenvalues")
plot(slrk$spectralclustering1$ei$vectors[, 1], ylab = "1st eigenvector")
plot(slrk$spectralclustering1$ei$vectors[, 2], ylab = "2nd eigenvector")
plot(slrk$spectralclustering1$ei$vectors[, 3], ylab = "3nd eigenvector")


##############################################################################
# slr method using k-means spectral clustering with K = 3
#   similarity.matrix -- TRUE
#   maxGamma -- TRUE
#   rank 1 approximation -- FALSE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
slr1sc = slr1sc(
  similarity.matrix = TRUE, maxGamma = TRUE,
  x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
  highdegree.regularization = FALSE)

slr1sc.coefs = getCoefsBM(
  coefs = coefficients(slr1sc$model), sbp = slr1sc$sbp)
slr1sc.coefs$llc.coefs
