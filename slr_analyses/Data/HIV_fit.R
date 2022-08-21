# Purpose: compare slr to other methods on data sets
# Date: 8/15/2022
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
scaling = TRUE

##############################################################################
# HIV: one of the HIV data sets in selbal package -- has binary response
#   n = 155 samples, 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 covariate (MSM), 
#   1 response (HIV_Status - binary)
W = selbal::HIV[, 1:60]
covar = data.frame(MSM = selbal::HIV[, 61])
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::HIV[, 62]
# levels(Y) # (control, case)
Y2 = ifelse(Y == "Pos", 1, 0)

################################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

################################################################################
# fit methods
################################################################################

# classo #######################################################################
# cl_gbm = codalasso(X_gbm, Y2, numFolds = K)
# saveRDS(
#   cl_gbm,
#   paste0(
#     output_dir, "/HIV",
#     "_classo",
#     "_gbm",
#     ".rds"))
cl = readRDS(
  paste0(
    output_dir, "/HIV",
    "_classo",
    "_gbm",
    ".rds"))

# slr - spectral ###############################################################
# slrspeccv = cv.slr(
#   x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "spectral",
#   response.type = "binary", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "auc",
#   scale = scaling, trace.it = FALSE)
# slrspec = slr(
#   x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "spectral",
#   response.type = "binary", s0.perc = 0, zeta = 0,
#   threshold = slrspeccv$threshold[slrspeccv$index["1se",]])
# saveRDS(
#   slrspeccv,
#   paste0(
#     output_dir, "/HIV",
#     "_slrcv_spectral",
#     "_gbm",
#     ".rds"))
# saveRDS(
#   slrspec,
#   paste0(
#     output_dir, "/HIV",
#     "_slr_spectral",
#     "_gbm",
#     ".rds"))

slrspeccv = readRDS(
  paste0(
    output_dir, "/HIV",
    "_slrcv_spectral",
    "_gbm",
    ".rds"))
slrspec = readRDS(
  paste0(
    output_dir, "/HIV",
    "_slr_spectral",
    "_gbm",
    ".rds"))

# slr - hierarchical ###########################################################
# slrhiercv = cv.slr(
#   x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "binary", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "auc",
#   scale = scaling, trace.it = FALSE)
# slrhier = slr(
#   x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "binary", s0.perc = 0, zeta = 0,
#   threshold = slrhiercv$threshold[slrhiercv$index["1se",]])
# saveRDS(
#   slrhiercv,
#   paste0(
#     output_dir, "/HIV",
#     "_slrcv_hierarchical",
#     "_gbm",
#     ".rds"))
# saveRDS(
#   slrhier,
#   paste0(
#     output_dir, "/HIV",
#     "_slr_hierarchical",
#     "_gbm",
#     ".rds"))

slrhiercv = readRDS(
  paste0(
    output_dir, "/HIV",
    "_slrcv_hierarchical",
    "_gbm",
    ".rds"))
slrhier = readRDS(
  paste0(
    output_dir, "/HIV",
    "_slr_hierarchical",
    "_gbm",
    ".rds"))

# selbal #######################################################################
# slbl_gbm = selbal::selbal.cv(x = X_gbm, y = Y, n.fold = K)
# saveRDS(
#   slbl_gbm,
#   paste0(
#     output_dir, "/HIV",
#     "_selbal",
#     "_gbm",
#     ".rds"))
slbl = readRDS(
  paste0(
    output_dir, "/HIV",
    "_selbal",
    "_gbm",
    ".rds"))

# selbal - covariate ###########################################################
# slbl_covar_gbm = selbal::selbal.cv(x = X_gbm, y = Y, n.fold = K, covar = covar)
# saveRDS(
#   slbl_covar_gbm,
#   paste0(
#     output_dir, "/HIV",
#     "_selbal_covar",
#     "_gbm",
#     ".rds"))
slblc = readRDS(
  paste0(
    output_dir, "/HIV",
    "_selbal_covar",
    "_gbm",
    ".rds"))

# codacore #####################################################################
# codacore0 = codacore::codacore(
#   x = X_gbm, y = Y2, logRatioType = "ILR",
#   objective = "binary classification", cvParams = list(numFolds = K))
# saveRDS(
#   codacore0,
#   paste0(
#     output_dir, "/HIV",
#     "_codacore",
#     "_gbm",
#     ".rds"))
cdcr = readRDS(
  paste0(
    output_dir, "/HIV",
    "_codacore",
    "_gbm",
    ".rds"))

# log-ratio lasso ############################################################
# library(logratiolasso)
# source("slr_analyses/Functions/logratiolasso.R")
# W_gbm.c = scale(log(X_gbm), center = TRUE, scale = FALSE)
# lrl <- cv_two_stage(
#   z = W_gbm.c, y = Y2, n_folds = K, family="binomial")
# saveRDS(
#   lrl,
#   paste0(
#     output_dir, "/HIV",
#     "_lrlasso",
#     "_gbm",
#     ".rds"))
lrl = readRDS(
  paste0(
    output_dir, "/HIV",
    "_lrlasso",
    "_gbm",
    ".rds"))

################################################################################
# get active sets and selected balances (if applicable)
################################################################################
p = ncol(X)

# classo #######################################################################
# selected variables
cl.betahat = cl$cll$betas[-1]
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
slbl$global.balance[slbl$global.balance$Group == "NUM", "Taxa"] # 1
slbl$global.balance[slbl$global.balance$Group == "DEN", "Taxa"] # 1

# selbal - covariate ###########################################################
# numerator (I+) / denominator (I-) of selected balance
slblc$global.balance[slblc$global.balance$Group == "NUM", "Taxa"] # 1
slblc$global.balance[slblc$global.balance$Group == "DEN", "Taxa"] # 1

# codacore #####################################################################
length(cdcr$ensemble) # one balance selected
# numerator (I+) / denominator (I-) of 1st selected balance
cdcr$ensemble[[1]]$slope # positive (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[1]]$hard$numerator]
colnames(X)[cdcr$ensemble[[1]]$hard$denominator]
# numerator (I+) / denominator (I-) of 2nd selected balance
cdcr$ensemble[[2]]$slope # positive (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[2]]$hard$denominator]
colnames(X)[cdcr$ensemble[[2]]$hard$numerator]

# log-ratio lasso ############################################################
# positive/negative effect on response
colnames(X)[lrl$beta_min > 0 & abs(lrl$beta_min) > 1e-8] # positive effect
colnames(X)[lrl$beta_min < 0 & abs(lrl$beta_min) > 1e-8] # negative effect
sum(abs(lrl$beta_min) > 1e-8)














