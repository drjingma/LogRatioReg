# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 5/2/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_slrs"

# # set up parallelization
# library(foreach)
# library(future)
# library(doFuture)
# library(parallel)
# registerDoFuture()
# nworkers = detectCores() / 2
# plan(multisession, workers = nworkers)
# 
# library(rngtools)
# library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Other simulation settings
numSims = 100

################################################################################
# Simulation Libraries and Settings #
################################################################################

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/slrtesting.R")
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

################################################################################
# Simulations #
################################################################################
b = 3

file.end = paste0(
  "_", sigma.settings,
  "_", paste0(
    paste(which(SBP.true == 1), collapse = ""), "v", 
    paste(which(SBP.true == -1), collapse = "")),
  "_dim", n, "x", p, 
  "_noisey", sigma_eps1, 
  "_noisex", sigma_eps2, 
  "_b0", b0, 
  "_b1", b1, 
  "_a0", a0, 
  "_theta", theta.value,
  "_sim", b,
  ".rds")

# ##############################################################################
# # generate data
# # get latent variable
# U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
# # simulate y from latent variable
# y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
# # simulate X: 
# epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
# a1 = theta.value * ilrtrans.true$ilr.trans[-p] 
# #   alpha1j = {
# #     c1=theta*ilr.const/k+   if j \in I+
# #     -c2=-theta*ilr.const/k-  if j \in I-
# #     0                       o/w
# #   }
# alrXj.all = a0 + U.all %*% t(a1) + epsj.all #log(Xj/Xp) =alpha0j+alpha1j*U+epsj
# X.all <- alrinv(alrXj.all)
# colnames(X.all) = paste0('s', 1:p)
# 
# # subset out training and test sets
# X = X.all[1:n, ]
# X.test = X.all[-(1:n), ]
# Y <- y.all[1:n]
# Y.test <- y.all[-(1:n)]
# 
# # about beta
# non0.beta = as.vector(SBP.true != 0)
# is0.beta = !non0.beta
# bspars = sum(non0.beta)
# # solve for beta
# c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
# beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
#   as.vector(ilrtrans.true$ilr.trans)

data.tmp = readRDS(paste0(output_dir, "/data", file.end))
X = data.tmp$X
Y = data.tmp$Y
X.test = data.tmp$X.test
Y.test = data.tmp$Y.test
SBP.true = data.tmp$SBP.true
beta.true = data.tmp$beta.true
is0.beta = data.tmp$is0.beta
non0.beta = data.tmp$non0.beta

##############################################################################
# compositional lasso (a linear log contrast method)
##############################################################################
# classo = cv.func(
#   method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), nlam = nlam,
#   nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
# 
# cl.lam.min.idx = which.min(classo$cvm)
# cl.a0 = classo$int[cl.lam.min.idx]
# cl.betahat = classo$bet[, cl.lam.min.idx]
# 
# # compute metrics on the selected model #
# cl.metrics = getMetricsLLC(
#   y.train = Y, y.test = Y.test,
#   logX.train = log(X),
#   logX.test = log(X.test),
#   n.train = n, n.test = n,
#   betahat0 = cl.a0, betahat = cl.betahat,
#   true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta, 
#   true.beta = beta.true)

classo_metrics = readRDS(paste0(output_dir, "/classo_metrics", file.end))
classo_metrics["Fscore"]

##############################################################################
# slr method using k-means spectral clustering with K = 3
#   rank 1 approximation -- FALSE
#   amini regularization -- FALSE
#   high degree regularization -- FALSE
##############################################################################
# slr0 = slr(
#   x = X, y = Y, approx = FALSE, amini.regularization = FALSE, 
#   highdegree.regularization = FALSE)
# 
# slr0.coefs = getCoefsBM(
#   coefs = coefficients(slr0$model), sbp = slr0$sbp)
# 
# # compute metrics on the selected model #
# slr0.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = slr0$sbp),
#   ilrX.test = getIlrX(X.test, sbp = slr0$sbp),
#   n.train = n, n.test = n,
#   thetahat0 = slr0.coefs$a0, thetahat = slr0.coefs$bm.coefs,
#   betahat = slr0.coefs$llc.coefs,
#   true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
#   true.beta = beta.true)

slr0_metrics = readRDS(paste0(output_dir, "/slr_metrics", file.end))
slr0_metrics["Fscore"]

slr0test = slr_testing(
  x = X, y = Y, approx = FALSE, amini.regularization = FALSE,
  highdegree.regularization = FALSE)
fields::image.plot(slr0test$kernel)
fields::image.plot(slr0test$spectralclustering$L$W)
fields::image.plot(slr0test$spectralclustering$L$W.tmp)
slr0test$cors
slr0test$Rsqs
slr0test$spectralclustering$cl
slr0test$sbp

##############################################################################
# slr method using k-means spectral clustering with K = 3
#   rank 1 approximation -- FALSE
#   amini regularization -- TRUE
#   high degree regularization -- FALSE
##############################################################################
# slr0am = slr(
#   x = X, y = Y, approx = FALSE, amini.regularization = TRUE, 
#   highdegree.regularization = FALSE)
# 
# slr0am.coefs = getCoefsBM(
#   coefs = coefficients(slr0am$model), sbp = slr0am$sbp)
# 
# # compute metrics on the selected model #
# slr0am.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = slr0am$sbp),
#   ilrX.test = getIlrX(X.test, sbp = slr0am$sbp),
#   n.train = n, n.test = n,
#   thetahat0 = slr0am.coefs$a0, thetahat = slr0am.coefs$bm.coefs,
#   betahat = slr0am.coefs$llc.coefs,
#   true.sbp = SBP.true, is0.true.beta = is0.beta, non0.true.beta = non0.beta,
#   true.beta = beta.true)

slr0am_metrics = readRDS(paste0(output_dir, "/slr_amini_metrics", file.end))
slr0am_metrics["Fscore"]


