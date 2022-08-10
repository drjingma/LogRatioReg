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
source("slr_analyses/Functions/semislrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

# Tuning parameters###########################################################

# Settings to toggle with
sigma.settings = "latentVarModel_missing"
n = 100
p = 30
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_y = 0.01
sigma_x = 0.1
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.5 # 0.5
theta.value = 0.5 # weight on a1 -- 1, 0.75, 0.5
a0 = 0 # 0
prop.missing = 0.5 # 0.25, 0.5
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
U.all = matrix(runif(min = -ulimit, max = ulimit, 2 * n), ncol = 1)
# simulate y from latent variable
y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_y)
# simulate X: 
epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_x
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

# missingness in training set
train.keep.upto = floor(n * (1 - prop.missing))
Y = Y[1:train.keep.upto]
X2 = X[(train.keep.upto + 1):n, ]
X = X[1:train.keep.upto, ]

# about beta
non0.beta = as.vector(SBP.true != 0)
bspars = sum(non0.beta)
# solve for beta
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
abs(c1plusc2^2 * varU - 2 * sigma_x^2)

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

# 
# ##############################################################################
# # compositional lasso
# # -- fits a linear log contrast model
# ##############################################################################
# start.time = Sys.time()
# classo = cv.func(
#   method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1),
#   nlam = nlam, nfolds = K, tol = tol, intercept = intercept,
#   scaling = scaling)
# end.time = Sys.time()
# cl.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# # cl.lam.idx = which.min(classo$cvm)
# oneSErule = min(classo$cvm) + classo$cvsd[which.min(classo$cvm)] * 1
# cl.lam.idx = which(classo$cvm <= oneSErule)[1]
# cl.a0 = classo$int[cl.lam.idx]
# cl.betahat = classo$bet[, cl.lam.idx]
# 
# # compute metrics on the selected model #
# cl.metrics = getMetricsLLC(
#   y.train = Y, y.test = Y.test,
#   logX.train = log(X),
#   logX.test = log(X.test),
#   n.train = n, n.test = n,
#   betahat0 = cl.a0, betahat = cl.betahat,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # slr
# #   screen.method = "wald"
# #   cluster.method = "spectral"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# slrspec0cv = cv.slr(
#   x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", threshold = NULL, s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   scale = scaling)
# slrspec0 = slr(
#   x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = slrspec0cv$threshold[slrspec0cv$index["1se",]])
# end.time = Sys.time()
# slrspec0.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# slrspec0.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(slrspec0.fullSBP) = colnames(X)
# slrspec0.fullSBP[match(
#   names(slrspec0$sbp), rownames(slrspec0.fullSBP))] = slrspec0$sbp
# 
# slrspec0.coefs = getCoefsBM(
#   coefs = coefficients(slrspec0$fit), sbp = slrspec0.fullSBP)
# 
# # compute metrics on the selected model #
# slrspec0.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = slrspec0.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = slrspec0.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = slrspec0.coefs$a0, thetahat = slrspec0.coefs$bm.coefs,
#   betahat = slrspec0.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # slr
# #   screen.method = "wald"
# #   cluster.method = "hieararchical"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# slrhier0cv = cv.slr(
#   x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# slrhier0 = slr(
#   x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = slrhier0cv$threshold[slrhier0cv$index["1se",]])
# end.time = Sys.time()
# slrhier0.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# slrhier0.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(slrhier0.fullSBP) = colnames(X)
# slrhier0.fullSBP[match(
#   names(slrhier0$sbp), rownames(slrhier0.fullSBP))] = slrhier0$sbp
# 
# slrhier0.coefs = getCoefsBM(
#   coefs = coefficients(slrhier0$fit), sbp = slrhier0.fullSBP)
# 
# # compute metrics on the selected model #
# slrhier0.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = slrhier0.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = slrhier0.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = slrhier0.coefs$a0, thetahat = slrhier0.coefs$bm.coefs,
#   betahat = slrhier0.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # semislr
# #   fold.x2 = FALSE
# #   screen.method = "wald"
# #   cluster.method = "spectral"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# sslrspec0cv = cv.semislr(
#   x = X, x2 = X2, fold.x2 = FALSE,
#   y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# sslrspec0 = semislr(
#   x = X, x2 = X2, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = sslrspec0cv$threshold[sslrspec0cv$index["1se",]])
# end.time = Sys.time()
# sslrspec0.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# sslrspec0.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(sslrspec0.fullSBP) = colnames(X)
# sslrspec0.fullSBP[match(
#   names(sslrspec0$sbp), rownames(sslrspec0.fullSBP))] = sslrspec0$sbp
# 
# sslrspec0.coefs = getCoefsBM(
#   coefs = coefficients(sslrspec0$fit), sbp = sslrspec0.fullSBP)
# 
# # compute metrics on the selected model #
# sslrspec0.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = sslrspec0.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = sslrspec0.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = sslrspec0.coefs$a0, thetahat = sslrspec0.coefs$bm.coefs,
#   betahat = sslrspec0.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # semislr
# #   fold.x2 = TRUE
# #   screen.method = "wald"
# #   cluster.method = "spectral"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# sslrspec1cv = cv.semislr(
#   x = X, x2 = X2, fold.x2 = TRUE,
#   y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# sslrspec1 = semislr(
#   x = X, x2 = X2, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = sslrspec1cv$threshold[sslrspec1cv$index["1se",]])
# end.time = Sys.time()
# sslrspec1.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# sslrspec1.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(sslrspec1.fullSBP) = colnames(X)
# sslrspec1.fullSBP[match(
#   names(sslrspec1$sbp), rownames(sslrspec1.fullSBP))] = sslrspec1$sbp
# 
# sslrspec1.coefs = getCoefsBM(
#   coefs = coefficients(sslrspec1$fit), sbp = sslrspec1.fullSBP)
# 
# # compute metrics on the selected model #
# sslrspec1.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = sslrspec1.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = sslrspec1.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = sslrspec1.coefs$a0, thetahat = sslrspec1.coefs$bm.coefs,
#   betahat = sslrspec1.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # semislr with no CV on x2
# #   screen.method = "wald"
# #   cluster.method = "spectral"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# sslrspec2cv = cv.slr(
#   x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# sslrspec2 = semislr(
#   x = X, x2 = X2, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = sslrspec2cv$threshold[sslrspec2cv$index["1se",]])
# end.time = Sys.time()
# sslrspec2.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# sslrspec2.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(sslrspec2.fullSBP) = colnames(X)
# sslrspec2.fullSBP[match(
#   names(sslrspec2$sbp), rownames(sslrspec2.fullSBP))] = sslrspec2$sbp
# 
# sslrspec2.coefs = getCoefsBM(
#   coefs = coefficients(sslrspec2$fit), sbp = sslrspec2.fullSBP)
# 
# # compute metrics on the selected model #
# sslrspec2.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = sslrspec2.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = sslrspec2.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = sslrspec2.coefs$a0, thetahat = sslrspec2.coefs$bm.coefs,
#   betahat = sslrspec2.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # semislr
# #   fold.x2 = FALSE
# #   screen.method = "wald"
# #   cluster.method = "hieararchical"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# sslrhier0cv = cv.semislr(
#   x = X, x2 = X2, fold.x2 = FALSE,
#   y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# sslrhier0 = semislr(
#   x = X, x2 = X2, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = sslrhier0cv$threshold[sslrhier0cv$index["1se",]])
# end.time = Sys.time()
# sslrhier0.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# sslrhier0.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(sslrhier0.fullSBP) = colnames(X)
# sslrhier0.fullSBP[match(
#   names(sslrhier0$sbp), rownames(sslrhier0.fullSBP))] = sslrhier0$sbp
# 
# sslrhier0.coefs = getCoefsBM(
#   coefs = coefficients(sslrhier0$fit), sbp = sslrhier0.fullSBP)
# 
# # compute metrics on the selected model #
# sslrhier0.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = sslrhier0.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = sslrhier0.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = sslrhier0.coefs$a0, thetahat = sslrhier0.coefs$bm.coefs,
#   betahat = sslrhier0.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # semislr
# #   fold.x2 = TRUE
# #   screen.method = "wald"
# #   cluster.method = "hieararchical"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# sslrhier1cv = cv.semislr(
#   x = X, x2 = X2, fold.x2 = TRUE,
#   y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# sslrhier1 = semislr(
#   x = X, x2 = X2, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = sslrhier1cv$threshold[sslrhier1cv$index["1se",]])
# end.time = Sys.time()
# sslrhier1.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# sslrhier1.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(sslrhier1.fullSBP) = colnames(X)
# sslrhier1.fullSBP[match(
#   names(sslrhier1$sbp), rownames(sslrhier1.fullSBP))] = sslrhier1$sbp
# 
# sslrhier1.coefs = getCoefsBM(
#   coefs = coefficients(sslrhier1$fit), sbp = sslrhier1.fullSBP)
# 
# # compute metrics on the selected model #
# sslrhier1.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = sslrhier1.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = sslrhier1.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = sslrhier1.coefs$a0, thetahat = sslrhier1.coefs$bm.coefs,
#   betahat = sslrhier1.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# ##############################################################################
# # semislr with no CV on x2
# #   screen.method = "wald"
# #   cluster.method = "hieararchical"
# #   response.type = "continuous"
# #   s0.perc = 0
# #   zeta = 0
# #   type.measure = "mse"
# # -- fits a balance regression model with one balance
# ##############################################################################
# start.time = Sys.time()
# sslrhier2cv = cv.slr(
#   x = X, y = Y, screen.method = "wald", cluster.method = "hieararchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# sslrhier2 = semislr(
#   x = X, x2 = X2, y = Y, screen.method = "wald", cluster.method = "hieararchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0, 
#   threshold = sslrhier2cv$threshold[sslrhier2cv$index["1se",]])
# end.time = Sys.time()
# sslrhier2.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# sslrhier2.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(sslrhier2.fullSBP) = colnames(X)
# sslrhier2.fullSBP[match(
#   names(sslrhier2$sbp), rownames(sslrhier2.fullSBP))] = sslrhier2$sbp
# 
# sslrhier2.coefs = getCoefsBM(
#   coefs = coefficients(sslrhier2$fit), sbp = sslrhier2.fullSBP)
# 
# # compute metrics on the selected model #
# sslrhier2.metrics = getMetricsBM(
#   y.train = Y, y.test = Y.test,
#   ilrX.train = getIlrX(X, sbp = sslrhier2.fullSBP),
#   ilrX.test = getIlrX(X.test, sbp = sslrhier2.fullSBP),
#   n.train = n, n.test = n,
#   thetahat0 = sslrhier2.coefs$a0, thetahat = sslrhier2.coefs$bm.coefs,
#   betahat = sslrhier2.coefs$llc.coefs,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true)
# 
# # ##############################################################################
# # # selbal method (a balance regression method)
# # # -- fits a balance regression model with one balance
# # ##############################################################################
# # library(selbal) # masks stats::cor()
# # slbl.data = getSelbalData(X = X, y = Y, classification = FALSE)
# # 
# # start.time = Sys.time()
# # slbl = selbal.cv(x = slbl.data$X, y = slbl.data$y, n.fold = K)
# # end.time = Sys.time()
# # slbl.timing = difftime(
# #   time1 = end.time, time2 = start.time, units = "secs")
# # 
# # slbl.coefs = getCoefsSelbal(
# #   X = slbl.data$X, y = slbl.data$y, selbal.fit = slbl, classification = FALSE, 
# #   check = TRUE)
# # 
# # # compute metrics on the selected model #
# # # prediction errors
# # # get prediction error on training set
# # slbl.Yhat.train = predict.glm(
# #   slbl$glm, 
# #   newdata = data.frame(V1 = balance::balance.fromSBP(
# #     x = slbl.data$X, y = slbl.coefs$sbp)), 
# #   type = "response")
# # slbl.PE.train = crossprod(Y - slbl.Yhat.train) / n
# # # get prediction error on test set
# # slbl.test.data = getSelbalData(X = X.test, y = Y.test, classification = FALSE)
# # slbl.Yhat.test = predict.glm(
# #   slbl$glm, 
# #   newdata = data.frame(V1 = balance::balance.fromSBP(
# #     x = slbl.test.data$X, y = slbl.coefs$sbp)), 
# #   type = "response")
# # slbl.PE.test = crossprod(Y.test - slbl.Yhat.test) / n
# # # beta estimation accuracy, selection accuracy #
# # slbl.metrics = getMetricsBM(
# #   thetahat = slbl.coefs$bm.coefs, betahat = slbl.coefs$llc.coefs,
# #   true.sbp = SBP.true, non0.true.beta = non0.beta,
# #   true.beta = llc.coefs.true, metrics = c("betaestimation", "selection"))
# # slbl.metrics = c(PEtr = slbl.PE.train, PEte = slbl.PE.test, slbl.metrics)
# 
# ##############################################################################
# # CoDaCoRe
# # -- fits a balance regression model with possibly multiple balances
# ##############################################################################
# library(codacore)
# 
# start.time = Sys.time()
# codacore0 = codacore(
#   x = X, y = Y, logRatioType = "ILR",
#   objective = "regression", cvParams = list(numFolds = K)) 
# end.time = Sys.time()
# codacore0.timing = difftime(
#   time1 = end.time, time2 = start.time, units = "secs")
# 
# if(length(codacore0$ensemble) > 0){ # at least 1 log-ratio found
#   codacore0_SBP = matrix(0, nrow = p, ncol = length(codacore0$ensemble))
#   codacore0_coeffs = rep(NA, length(codacore0$ensemble))
#   for(col.idx in 1:ncol(codacore0_SBP)){
#     codacore0_SBP[
#       codacore0$ensemble[[col.idx]]$hard$numerator, col.idx] = 1
#     codacore0_SBP[
#       codacore0$ensemble[[col.idx]]$hard$denominator, col.idx] = -1
#     codacore0_coeffs[col.idx] = codacore0$ensemble[[col.idx]]$slope
#   }
#   
#   codacore0.betahat = getBetaFromCodacore(
#     SBP_codacore = codacore0_SBP, coeffs_codacore = codacore0_coeffs, p = p)
#   
#   # compute metrics on the selected model #
#   # prediction errors
#   # get prediction error on training set
#   codacore0.Yhat.train = predict(codacore0, X)
#   # get prediction error on test set
#   codacore0.Yhat.test = predict(codacore0, X.test)
#   
# } else{
#   print(paste0("sim ", i, " -- codacore has no log-ratios"))
#   codacore0_coeffs = c()
#   codacore0model = stats::glm(Y ~ 1, family = "gaussian")
#   codacore0.betahat = rep(0, p)
#   
#   # compute metrics on the selected model #
#   # prediction errors
#   # get prediction error on training set
#   codacore0.Yhat.train = predict(codacore0model, X)
#   # get prediction error on test set
#   codacore0.Yhat.test = predict(codacore0model, X.test)
# }
# codacore0.PE.train = crossprod(Y - codacore0.Yhat.train) / n
# codacore0.PE.test = crossprod(Y.test - codacore0.Yhat.test) / n
# 
# # beta estimation accuracy, selection accuracy #
# codacore0.metrics = getMetricsBM(
#   betahat = codacore0.betahat,
#   true.sbp = SBP.true, non0.true.beta = non0.beta,
#   true.beta = llc.coefs.true, metrics = c("betaestimation", "selection"))
# codacore0.metrics = c(
#   PEtr = codacore0.PE.train, PEte = codacore0.PE.test, codacore0.metrics)
# 
# # ##############################################################################
# # # Log-Ratio Lasso
# # # -- regresses on pairwise log-ratios
# # ##############################################################################
# # library(logratiolasso)
# # source("slr_analyses/Functions/logratiolasso.R")
# # Wc = scale(log(X), center = TRUE, scale = FALSE)
# # Yc = Y - mean(Y)
# # 
# # start.time = Sys.time()
# # lrl <- cv_two_stage(z = Wc, y = Yc, n_folds = K)
# # end.time = Sys.time()
# # lrl.timing = difftime(
# #   time1 = end.time, time2 = start.time, units = "secs")
# # 
# # # compute metrics on the selected model #
# # # prediction errors
# # # get prediction error on training set
# # lrl.Yhat.train = Wc %*% lrl$beta_min
# # lrl.PE.train = crossprod(Y - lrl.Yhat.train) / n
# # # get prediction error on test set
# # Wc.test = scale(log(X.test), center = TRUE, scale = FALSE)
# # Yc.test = Y.test - mean(Y.test)
# # lrl.Yhat.test = Wc.test %*% lrl$beta_min
# # lrl.PE.test = crossprod(Yc.test - lrl.Yhat.test) / n
# # 
# # # beta estimation accuracy, selection accuracy #
# # lrl.metrics = getMetricsBM(
# #   betahat = lrl$beta_min, # don't back-scale bc only centered X (didn't scale)
# #   true.sbp = SBP.true, non0.true.beta = non0.beta,
# #   true.beta = llc.coefs.true, metrics = c("betaestimation", "selection"))
# # lrl.metrics = c(
# #   PEtr = lrl.PE.train, PEte = lrl.PE.test, lrl.metrics)
# 
# ##############################################################################
# ##############################################################################
# ##############################################################################
# ### fin ###
# 
# cbind(
#   cl = cl.metrics, 
#   slr.s = slrspec0.metrics, 
#   slr.h = slrhier0.metrics, 
#   semislr.s = sslrspec0.metrics,
#   semislr.h = sslrhier0.metrics,
#   semislr.s.f = sslrspec1.metrics,
#   semislr.h.f = sslrhier1.metrics, 
#   cdcr = codacore0.metrics
#   )
# 
# # pheatmap::pheatmap(slrspec0$cluster.mat)
# pheatmap::pheatmap(slrspec0$cluster.mat)
# pheatmap::pheatmap(sslrspec0$cluster.mat)
# pheatmap::pheatmap(sslrspec1$cluster.mat)
# fields::image.plot(slrspec0$cluster.mat)
# fields::image.plot(sslrspec0$cluster.mat)
# fields::image.plot(sslrspec1$cluster.mat)
# 
# pheatmap::pheatmap(slrhier0$cluster.mat)
# pheatmap::pheatmap(sslrhier0$cluster.mat)
# pheatmap::pheatmap(sslrhier1$cluster.mat)
# fields::image.plot(slrhier0$cluster.mat)
# fields::image.plot(sslrhier0$cluster.mat)
# fields::image.plot(sslrhier1$cluster.mat)

##############################################################################
##############################################################################
##############################################################################
### fin ###
# }


# Aitchison variation when j != k
# when j != k, aitchison var is Sjk = (c1 + c2)^2 Var[U] + 2 * sigma_eps2
c1plusc2^2 * varU # term 1 # or wnant this to be small????
2 * sigma_x^2 # term 2 (want this term to dominate)
abs(c1plusc2^2 * varU - 2 * sigma_x^2)

# Correlation bt clr(Xj) & y
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
