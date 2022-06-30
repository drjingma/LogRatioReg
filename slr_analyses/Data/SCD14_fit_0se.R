# Purpose: compare slr to other methods on data sets
# Date: 6/29/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs"

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(selbal)

source("RCode/func_libs.R")
source("slr_analyses/Functions/slrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

# tuning parameter settings
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4

################################################################################
# sCD14: another HIV data set in selbal package
#   n = 151 samples (a subset from sCD14 data set), 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
W = selbal::sCD14[, 1:60]
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::sCD14[, 61]

################################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

################################################################################
# fit methods
################################################################################
p = ncol(X)

# classo #######################################################################
# classo = cv.func(
#   method="ConstrLasso", y = Y, x = log(X_gbm), Cmat = matrix(1, p, 1),
#   nlam = nlam, nfolds = K, tol = tol, intercept = intercept,
#   scaling = scaling)
# saveRDS(
#   classo,
#   paste0(
#     output_dir, "/sCD14_0se",
#     "_classo",
#     "_gbm",
#     ".rds"))
cl = readRDS(
  paste0(
    output_dir, "/sCD14_0se",
    "_classo",
    "_gbm",
    ".rds"))

# slr - spectral ###############################################################
# slrspeccv = cv.slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "mse",
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# slrspec = slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   threshold = slrspeccv$threshold[slrspeccv$index["min",]])
# saveRDS(
#   slrspec,
#   paste0(
#     output_dir, "/sCD14_0se",
#     "_slr_spectral",
#     "_gbm",
#     ".rds"))
slrspec = readRDS(
  paste0(
    output_dir, "/sCD14_0se",
    "_slr_spectral",
    "_gbm",
    ".rds"))

# slr - hierarchical ###########################################################
# slrhiercv = cv.slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "mse",
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# slrhier = slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   threshold = slrhiercv$threshold[slrhiercv$index["min",]])
# saveRDS(
#   slrhier,
#   paste0(
#     output_dir, "/sCD14_0se",
#     "_slr_hierarchical",
#     "_gbm",
#     ".rds"))
slrhier = readRDS(
  paste0(
    output_dir, "/sCD14_0se",
    "_slr_hierarchical",
    "_gbm",
    ".rds"))

# selbal #######################################################################
# slbl = selbal::selbal.cv(x = X_gbm, y = Y, n.fold = K, opt.cri = "max")
# saveRDS(
#   slbl,
#   paste0(
#     output_dir, "/sCD14_0se",
#     "_selbal",
#     "_gbm",
#     ".rds"))
slbl = readRDS(
  paste0(
    output_dir, "/sCD14_0se",
    "_selbal",
    "_gbm",
    ".rds"))

# codacore #####################################################################
# codacore0 = codacore::codacore(
#   x = X_gbm, y = Y, logRatioType = "ILR",
#   objective = "regression", cvParams = list(numFolds = K), lambda = 0)
# saveRDS(
#   codacore0,
#   paste0(
#     output_dir, "/sCD14_0se",
#     "_codacore",
#     "_gbm",
#     ".rds"))
cdcr = readRDS(
  paste0(
    output_dir, "/sCD14_0se",
    "_codacore",
    "_gbm",
    ".rds"))

# log-ratio lasso ##############################################################
# library(logratiolasso)
# source("slr_analyses/Functions/logratiolasso.R")
# W_gbm.c = scale(log(X_gbm), center = TRUE, scale = FALSE)
# lrl <- cv_two_stage(
#   z = W_gbm.c, y = Y, n_folds = K, family="gaussian")
# saveRDS(
#   lrl,
#   paste0(
#     output_dir, "/sCD14_0se",
#     "_lrlasso",
#     "_gbm",
#     ".rds"))
lrl = readRDS(
  paste0(
    output_dir, "/sCD14_0se",
    "_lrlasso",
    "_gbm",
    ".rds"))

################################################################################
# get active sets and selected balances (if applicable)
################################################################################

# classo #######################################################################
# selected variables
# oneSErule = min(cl$cvm) + cl$cvsd[which.min(cl$cvm)] * 1
# cl.lam.idx = which(cl$cvm <= oneSErule)[1]
cl.lam.idx = which.min(cl$cvm)
cl.a0 = cl$int[cl.lam.idx]
cl.betahat = cl$bet[, cl.lam.idx]
# positive/negative effect on response
colnames(X)[cl.betahat > 0 & abs(cl.betahat) > 1e-8] # positive effect
colnames(X)[cl.betahat < 0 & abs(cl.betahat) > 1e-8] # negative effect
sum(abs(cl.betahat) > 1e-8)

# slr - spectral ###############################################################
# SBP
slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrspec.fullSBP) = colnames(X)
slrspec.fullSBP[match(
  names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp
# thetahat 
slrspec.coefs = getCoefsBM(
  coefs = coefficients(slrspec$fit), sbp = slrspec.fullSBP)
# numerator (I+) / denominator (I-) of selected balance
rownames(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs > 0]
rownames(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs < 0]
sum(slrspec.fullSBP[, 1] != 0)

# slr - hierarchical ###########################################################
# SBP
slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrhier.fullSBP) = colnames(X)
slrhier.fullSBP[match(
  names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
# thetahat 
slrhier.coefs = getCoefsBM(
  coefs = coefficients(slrhier$fit), sbp = slrhier.fullSBP)
# numerator (I+) / denominator (I-) of selected balance
rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs > 0]
rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs < 0]
sum(slrhier.fullSBP[, 1] != 0)

# selbal #######################################################################
# numerator (I+) / denominator (I-) of selected balance
slbl$global.balance[slbl$global.balance$Group == "NUM", "Taxa"] # 4
slbl$global.balance[slbl$global.balance$Group == "DEN", "Taxa"] # 4

# codacore #####################################################################
length(cdcr$ensemble) # one balance selected
# numerator (I+) / denominator (I-) of selected balance
cdcr$ensemble[[1]]$slope # positive (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[1]]$hard$numerator]
colnames(X)[cdcr$ensemble[[1]]$hard$denominator]
# 2nd balance
cdcr$ensemble[[2]]$slope # negative (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[2]]$hard$denominator]
colnames(X)[cdcr$ensemble[[2]]$hard$numerator]
# 3rd balance
cdcr$ensemble[[3]]$slope # positive (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[3]]$hard$numerator]
colnames(X)[cdcr$ensemble[[3]]$hard$denominator]
# 4th balance
cdcr$ensemble[[4]]$slope # negative (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[4]]$hard$denominator]
colnames(X)[cdcr$ensemble[[4]]$hard$numerator]
# 4th balance
cdcr$ensemble[[5]]$slope # negative (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[5]]$hard$denominator]
colnames(X)[cdcr$ensemble[[5]]$hard$numerator]

# # log-ratio lasso ############################################################
# # positive/negative effect on response
# colnames(X)[lrl$beta_min > 0 & abs(lrl$beta_min) > 1e-8] # positive effect
# colnames(X)[lrl$beta_min < 0 & abs(lrl$beta_min) > 1e-8] # negative effect
# sum(abs(lrl$beta_min) > 1e-8)















