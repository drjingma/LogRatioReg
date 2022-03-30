# Purpose: compare slr to selbal on data sets
# Date: 3/29/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "Kristyn/Data/outputs"

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(propr)
library(selbal)

library(pROC)

library(zCompositions)
library(pheatmap)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/codalasso.R")

# helper functions
source("Kristyn/Functions/util.R")

# tuning parameter settings
K = 10



################################################################################
# sCD14: another HIV data set in selbal package
#   n = 151 samples (a subset from sCD14 data set), 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
W = selbal::sCD14[, 1:60]
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::sCD14[, 61]
################################################################################
# need to handle 0's



################################################################################
# Strategy 1: Replace 0's with 0.5 (used in Lin et al. 2014 [classo])
W_0.5 = W
W_0.5[W_0.5 == 0] = 0.5
X_0.5 = sweep(W_0.5, 1, rowSums(W_0.5), FUN='/')

slrcor_0.5 = slrmatrix(x = X_0.5, y = Y)
pheatmap(slrcor_0.5)

##### slr with approximation step #####
# slr0approx_0.5 = slr(
#   x = X_0.5, y = Y, num.clusters = 2, classification = FALSE, approx = TRUE)
# saveRDS(
#   slr0approx_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_slr_approx",
#     "_0.5",
#     ".rds"))
slr0approx_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slr_approx",
    "_0.5",
    ".rds"))
slr0approx_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0approx_0.5$model), sbp = slr0approx_0.5$sbp)
rownames(slr0approx_0.5_coefs$llc.coefs)[slr0approx_0.5_coefs$llc.coefs != 0]

##### slr (no approximation step) #####
# slr0_0.5 = slr(
#   x = X_0.5, y = Y, num.clusters = 2, classification = FALSE, approx = TRUE)
# saveRDS(
#   slr0_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_slr",
#     "_0.5",
#     ".rds"))
slr0_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slr",
    "_0.5",
    ".rds"))
slr0_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0_0.5$model), sbp = slr0_0.5$sbp)
rownames(slr0_0.5_coefs$llc.coefs)[slr0_0.5_coefs$llc.coefs != 0]

##### selbal #####
# slbl_0.5 = selbal.cv(x = X_0.5, y = Y, n.fold = K)
# saveRDS(
#   slbl_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_selbal",
#     "_0.5",
#     ".rds"))
slbl_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_selbal",
    "_0.5",
    ".rds"))
slbl_0.5_coefs = getCoefsSelbal(
  X = X_0.5, y = Y, selbal.fit = slbl_0.5, classification = FALSE, check = TRUE)
rownames(slbl_0.5_coefs$llc.coefs)[slbl_0.5_coefs$llc.coefs != 0]



################################################################################
# Strategy 2: GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

slrcor_gbm = slrmatrix(x = X_gbm, y = Y)
pheatmap(slrcor_gbm)

##### slr with approximation step #####
# slr0approx_gbm = slr(
#   x = X_gbm, y = Y, num.clusters = 2, classification = FALSE, approx = TRUE)
# saveRDS(
#   slr0approx_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_slr_approx",
#     "_gbm",
#     ".rds"))
slr0approx_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slr_approx",
    "_gbm",
    ".rds"))
slr0approx_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0approx_gbm$model), sbp = slr0approx_gbm$sbp)
rownames(slr0approx_gbm_coefs$llc.coefs)[slr0approx_gbm_coefs$llc.coefs != 0]

##### slr (no approximation step) #####
# slr0_gbm = slr(
#   x = X_gbm, y = Y, num.clusters = 2, classification = FALSE, approx = TRUE)
# saveRDS(
#   slr0_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_slr",
#     "_gbm",
#     ".rds"))
slr0_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slr",
    "_gbm",
    ".rds"))
slr0_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0_gbm$model), sbp = slr0_gbm$sbp)
rownames(slr0_gbm_coefs$llc.coefs)[slr0_gbm_coefs$llc.coefs != 0]

##### selbal #####
# slbl_gbm = selbal.cv(x = X_gbm, y = Y, n.fold = K)
# saveRDS(
#   slbl_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_selbal",
#     "_gbm",
#     ".rds"))
slbl_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_selbal",
    "_gbm",
    ".rds"))
slbl_gbm_coefs = getCoefsSelbal(
  X = X_gbm, y = Y, selbal.fit = slbl_gbm, classification = FALSE, check = TRUE)
rownames(slbl_gbm_coefs$llc.coefs)[slbl_gbm_coefs$llc.coefs != 0]


