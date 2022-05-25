# Purpose: compare slr to selbal on data sets
# Date: 3/30/2022
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
sum(slr0approx_0.5_coefs$llc.coefs != 0)

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
sum(slr0_0.5_coefs$llc.coefs != 0)

##### cv.slr with approximation step #####
# slrcv0approx_0.5 = cv.slr(
#   x = X_0.5, y = Y, max.clusters = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = TRUE)
# saveRDS(
#   slrcv0approx_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_slrcv_approx",
#     "_0.5",
#     ".rds"))
slrcv0approx_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slrcv_approx",
    "_0.5",
    ".rds"))
slrcv0approx_0.5$nclusters_1se_idx
slrcv0approx_0.5$nclusters_min_idx
slrcv0approx_0.5_cvdata = data.frame(
  `T` = slrcv0approx_0.5$nclusters, 
  cvm = slrcv0approx_0.5$cvm, 
  cvse = slrcv0approx_0.5$cvse
)
ggplot(slrcv0approx_0.5_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(slrcv0approx_0.5_cvdata$`T`))
slrcv0approx_0.5_selected = slrcv0approx_0.5$models[[
  slrcv0approx_0.5$nclusters_1se_idx]]
slrcv0approx_0.5_coefs = getCoefsBM(
  coefs = coefficients(slrcv0approx_0.5_selected$model), 
  sbp = slrcv0approx_0.5_selected$sbp)
rownames(slrcv0approx_0.5_coefs$llc.coefs)[
  slrcv0approx_0.5_coefs$llc.coefs != 0]

##### cv.slr (no approximation step) #####
# slrcv0_0.5 = cv.slr(
#   x = X_0.5, y = Y, max.clusters = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = FALSE)
# saveRDS(
#   slrcv0_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_slrcv",
#     "_0.5",
#     ".rds"))
slrcv0_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slrcv",
    "_0.5",
    ".rds"))
slrcv0_0.5$nclusters_1se_idx
slrcv0_0.5$nclusters_min_idx
slrcv0_0.5_cvdata = data.frame(
  `T` = slrcv0_0.5$nclusters, 
  cvm = slrcv0_0.5$cvm, 
  cvse = slrcv0_0.5$cvse
)
ggplot(slrcv0_0.5_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(slrcv0_0.5_cvdata$`T`))
slrcv0_0.5_selected = slrcv0_0.5$models[[
  slrcv0_0.5$nclusters_1se_idx]]
slrcv0_0.5_coefs = getCoefsBM(
  coefs = coefficients(slrcv0_0.5_selected$model), 
  sbp = slrcv0_0.5_selected$sbp)
rownames(slrcv0_0.5_coefs$llc.coefs)[
  slrcv0_0.5_coefs$llc.coefs != 0]

##### cv.hslr with approximation step #####
# hslrcv0approx_0.5 = cv.hslr(
#   x = X_0.5, y = Y, max.levels = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = TRUE)
# saveRDS(
#   hslrcv0approx_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_hslrcv_approx",
#     "_0.5",
#     ".rds"))
hslrcv0approx_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_hslrcv_approx",
    "_0.5",
    ".rds"))
hslrcv0approx_0.5$nclusters_1se_idx
hslrcv0approx_0.5$nclusters_min_idx
hslrcv0approx_0.5_cvdata = data.frame(
  `T` = hslrcv0approx_0.5$nclusters, 
  cvm = hslrcv0approx_0.5$cvm, 
  cvse = hslrcv0approx_0.5$cvse
)
ggplot(hslrcv0approx_0.5_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(hslrcv0approx_0.5_cvdata$`T`))
hslrcv0approx_0.5_selected = hslrcv0approx_0.5$models[[
  hslrcv0approx_0.5$nclusters_1se_idx]]
hslrcv0approx_0.5_coefs = getCoefsBM(
  coefs = coefficients(hslrcv0approx_0.5_selected$model), 
  sbp = hslrcv0approx_0.5_selected$sbp)
rownames(hslrcv0approx_0.5_coefs$llc.coefs)[
  hslrcv0approx_0.5_coefs$llc.coefs != 0]

##### cv.hslr (no approximation step) #####
# hslrcv0_0.5 = cv.hslr(
#   x = X_0.5, y = Y, max.levels = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = FALSE)
# saveRDS(
#   hslrcv0_0.5,
#   paste0(
#     output_dir, "/sCD14",
#     "_hslrcv",
#     "_0.5",
#     ".rds"))
hslrcv0_0.5 = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_hslrcv",
    "_0.5",
    ".rds"))
hslrcv0_0.5$nclusters_1se_idx
hslrcv0_0.5$nclusters_min_idx
hslrcv0_0.5_cvdata = data.frame(
  `T` = hslrcv0_0.5$nclusters, 
  cvm = hslrcv0_0.5$cvm, 
  cvse = hslrcv0_0.5$cvse
)
ggplot(hslrcv0_0.5_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(hslrcv0_0.5_cvdata$`T`))
hslrcv0_0.5_selected = hslrcv0_0.5$models[[
  hslrcv0_0.5$nclusters_1se_idx]]
hslrcv0_0.5_coefs = getCoefsBM(
  coefs = coefficients(hslrcv0_0.5_selected$model), 
  sbp = hslrcv0_0.5_selected$sbp)
rownames(hslrcv0_0.5_coefs$llc.coefs)[
  hslrcv0_0.5_coefs$llc.coefs != 0]

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
sum(slbl_0.5_coefs$llc.coefs != 0)


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
sum(slr0approx_gbm_coefs$llc.coefs != 0)

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
sum(slr0_gbm_coefs$llc.coefs != 0)

##### cv.slr with approximation step #####
# slrcv0approx_gbm = cv.slr(
#   x = X_gbm, y = Y, max.clusters = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = TRUE)
# saveRDS(
#   slrcv0approx_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_slrcv_approx",
#     "_gbm",
#     ".rds"))
slrcv0approx_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slrcv_approx",
    "_gbm",
    ".rds"))
slrcv0approx_gbm$nclusters_1se_idx
slrcv0approx_gbm$nclusters_min_idx
slrcv0approx_gbm_cvdata = data.frame(
  `T` = slrcv0approx_gbm$nclusters, 
  cvm = slrcv0approx_gbm$cvm, 
  cvse = slrcv0approx_gbm$cvse
)
ggplot(slrcv0approx_gbm_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(slrcv0approx_gbm_cvdata$`T`))
slrcv0approx_gbm_selected = slrcv0approx_gbm$models[[
  slrcv0approx_gbm$nclusters_1se_idx]]
slrcv0approx_gbm_coefs = getCoefsBM(
  coefs = coefficients(slrcv0approx_gbm_selected$model), 
  sbp = slrcv0approx_gbm_selected$sbp)
rownames(slrcv0approx_gbm_coefs$llc.coefs)[
  slrcv0approx_gbm_coefs$llc.coefs != 0]

##### cv.slr (no approximation step) #####
# slrcv0_gbm = cv.slr(
#   x = X_gbm, y = Y, max.clusters = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = FALSE)
# saveRDS(
#   slrcv0_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_slrcv",
#     "_gbm",
#     ".rds"))
slrcv0_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slrcv",
    "_gbm",
    ".rds"))
slrcv0_gbm$nclusters_1se_idx
slrcv0_gbm$nclusters_min_idx
slrcv0_gbm_cvdata = data.frame(
  `T` = slrcv0_gbm$nclusters, 
  cvm = slrcv0_gbm$cvm, 
  cvse = slrcv0_gbm$cvse
)
ggplot(slrcv0_gbm_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(slrcv0_gbm_cvdata$`T`))
slrcv0_gbm_selected = slrcv0_gbm$models[[
  slrcv0_gbm$nclusters_1se_idx]]
slrcv0_gbm_coefs = getCoefsBM(
  coefs = coefficients(slrcv0_gbm_selected$model), 
  sbp = slrcv0_gbm_selected$sbp)
rownames(slrcv0_gbm_coefs$llc.coefs)[
  slrcv0_gbm_coefs$llc.coefs != 0]

##### cv.hslr with approximation step #####
# hslrcv0approx_gbm = cv.hslr(
#   x = X_gbm, y = Y, max.levels = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = TRUE)
# saveRDS(
#   hslrcv0approx_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_hslrcv_approx",
#     "_gbm",
#     ".rds"))
hslrcv0approx_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_hslrcv_approx",
    "_gbm",
    ".rds"))
hslrcv0approx_gbm$nclusters_1se_idx
hslrcv0approx_gbm$nclusters_min_idx
hslrcv0approx_gbm_cvdata = data.frame(
  `T` = hslrcv0approx_gbm$nclusters, 
  cvm = hslrcv0approx_gbm$cvm, 
  cvse = hslrcv0approx_gbm$cvse
)
ggplot(hslrcv0approx_gbm_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(hslrcv0approx_gbm_cvdata$`T`))
hslrcv0approx_gbm_selected = hslrcv0approx_gbm$models[[
  hslrcv0approx_gbm$nclusters_1se_idx]]
hslrcv0approx_gbm_coefs = getCoefsBM(
  coefs = coefficients(hslrcv0approx_gbm_selected$model), 
  sbp = hslrcv0approx_gbm_selected$sbp)
rownames(hslrcv0approx_gbm_coefs$llc.coefs)[
  hslrcv0approx_gbm_coefs$llc.coefs != 0]

##### cv.hslr (no approximation step) #####
# hslrcv0_gbm = cv.hslr(
#   x = X_gbm, y = Y, max.levels = ceiling(sqrt(ncol(X))), nfolds = K,
#   classification = FALSE, approx = FALSE)
# saveRDS(
#   hslrcv0_gbm,
#   paste0(
#     output_dir, "/sCD14",
#     "_hslrcv",
#     "_gbm",
#     ".rds"))
hslrcv0_gbm = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_hslrcv",
    "_gbm",
    ".rds"))
hslrcv0_gbm$nclusters_1se_idx
hslrcv0_gbm$nclusters_min_idx
hslrcv0_gbm_cvdata = data.frame(
  `T` = hslrcv0_gbm$nclusters, 
  cvm = hslrcv0_gbm$cvm, 
  cvse = hslrcv0_gbm$cvse
)
ggplot(hslrcv0_gbm_cvdata, aes(x = `T`, y = cvm)) + 
  geom_path() + 
  geom_point() + 
  geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.2) + 
  scale_x_continuous(breaks = c(hslrcv0_gbm_cvdata$`T`))
hslrcv0_gbm_selected = hslrcv0_gbm$models[[
  hslrcv0_gbm$nclusters_1se_idx]]
hslrcv0_gbm_coefs = getCoefsBM(
  coefs = coefficients(hslrcv0_gbm_selected$model), 
  sbp = hslrcv0_gbm_selected$sbp)
rownames(hslrcv0_gbm_coefs$llc.coefs)[
  hslrcv0_gbm_coefs$llc.coefs != 0]

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
sum(slbl_gbm_coefs$llc.coefs != 0)

