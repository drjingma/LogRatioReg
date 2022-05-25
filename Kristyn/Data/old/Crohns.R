# Purpose: compare slr to selbal on data sets
# Date: 4/13/2022
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

library(ggplot2)
library(ggpubr)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/codalasso.R")

# helper functions
source("Kristyn/Functions/util.R")

getTPlots = function(cvslr_fit){
  cv_data = data.frame(
    `T` = cvslr_fit$nclusters, 
    cvm = cvslr_fit$cvm, 
    cvse = cvslr_fit$cvse,
    Active.Size = sapply(
      cvslr_fit$models, function(elt) sum(elt$sbp != 0))
  )
  cvm.plt = ggplot(cv_data, aes(x = `T`, y = cvm)) + 
    geom_path() + 
    geom_point() + 
    geom_errorbar(aes(ymin = cvm - cvse, ymax = cvm + cvse), width = 0.5) + 
    scale_x_continuous(breaks = c(cv_data$`T`))
  size.plt = ggplot(cv_data, aes(x = `T`, y = Active.Size)) + 
    geom_path() + 
    geom_point() + 
    scale_x_continuous(breaks = c(cv_data$`T`))
  return(ggarrange(cvm.plt, size.plt, nrow = 1))
}

# for hslr tree plot
library(ggraph) # make tree/dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidyverse)
library(tidygraph)

# tuning parameter settings
K = 10
slrmax = 10


################################################################################
# Crohn: a data set in selbal package
#   n = 975 samples, 
#   p = 48 taxa (counts for microbial taxa at genus level), 
#   1 response (y - binary)
W = selbal::Crohn[, 1:48]
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::Crohn[, 49]
Y2 = ifelse(Y == "CD", 1, 0)
################################################################################
# need to handle 0's




################################################################################
# Strategy 1: Replace 0's with 0.5 (used in Lin et al. 2014 [classo])
W_0.5 = W
W_0.5[W_0.5 == 0] = 0.5
X_0.5 = sweep(W_0.5, 1, rowSums(W_0.5), FUN='/')

slrcor_0.5 = slrmatrix(x = X_0.5, y = Y2)
pheatmap(slrcor_0.5)

# ##### slr with approximation step #####
# # slr0approx_0.5 = slr(
# #   x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, approx = TRUE)
# # saveRDS(
# #   slr0approx_0.5,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slr_approx",
# #     "_0.5",
# #     ".rds"))
# slr0approx_0.5 = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slr_approx",
#     "_0.5",
#     ".rds"))
# slr0approx_0.5_coefs = getCoefsBM(
#   coefs = coefficients(slr0approx_0.5$model), sbp = slr0approx_0.5$sbp)
# rownames(slr0approx_0.5_coefs$llc.coefs)[slr0approx_0.5_coefs$llc.coefs != 0]
# sum(slr0approx_0.5_coefs$llc.coefs != 0)
# slr0approx_0.5$sbp[slr0approx_0.5$sbp != 0, , drop = FALSE]
# 
# ##### slr (no approximation step) #####
# # slr0_0.5 = slr(
# #   x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, approx = TRUE)
# # saveRDS(
# #   slr0_0.5,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slr",
# #     "_0.5",
# #     ".rds"))
# slr0_0.5 = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slr",
#     "_0.5",
#     ".rds"))
# slr0_0.5_coefs = getCoefsBM(
#   coefs = coefficients(slr0_0.5$model), sbp = slr0_0.5$sbp)
# rownames(slr0_0.5_coefs$llc.coefs)[slr0_0.5_coefs$llc.coefs != 0]
# sum(slr0_0.5_coefs$llc.coefs != 0)
# slr0_0.5$sbp[slr0_0.5$sbp != 0, , drop = FALSE]
# 
# ##### cv.slr with approximation step #####
# # slrcv0approx_0.5 = cv.slr(
# #   x = X_0.5, y = Y2, max.clusters = slrmax, nfolds = K,
# #   classification = TRUE, approx = TRUE)
# # saveRDS(
# #   slrcv0approx_0.5,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slrcv_approx",
# #     "_0.5",
# #     ".rds"))
# slrcv0approx_0.5 = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slrcv_approx",
#     "_0.5",
#     ".rds"))
# slrcv0approx_0.5$nclusters_1se_idx
# slrcv0approx_0.5$nclusters_min_idx
# getTPlots(slrcv0approx_0.5)
# slrcv0approx_0.5_selected = slrcv0approx_0.5$models[[
#   slrcv0approx_0.5$nclusters_1se_idx]]
# slrcv0approx_0.5_coefs = getCoefsBM(
#   coefs = coefficients(slrcv0approx_0.5_selected$model),
#   sbp = slrcv0approx_0.5_selected$sbp)
# rownames(slrcv0approx_0.5_coefs$llc.coefs)[
#   slrcv0approx_0.5_coefs$llc.coefs != 0]
# sum(slrcv0approx_0.5_coefs$llc.coefs != 0)
# slrcv0approx_0.5$models[[slrcv0approx_0.5$nclusters_1se_idx]]$sbp[
#   slrcv0approx_0.5$models[[slrcv0approx_0.5$nclusters_1se_idx]]$sbp != 0, ,
#   drop = FALSE]
# 
# ##### cv.slr (no approximation step) #####
# # slrcv0_0.5 = cv.slr(
# #   x = X_0.5, y = Y2, max.clusters = slrmax, nfolds = K,
# #   classification = TRUE, approx = FALSE)
# # saveRDS(
# #   slrcv0_0.5,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slrcv",
# #     "_0.5",
# #     ".rds"))
# slrcv0_0.5 = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slrcv",
#     "_0.5",
#     ".rds"))
# slrcv0_0.5$nclusters_1se_idx
# slrcv0_0.5$nclusters_min_idx
# getTPlots(slrcv0_0.5)
# slrcv0_0.5_selected = slrcv0_0.5$models[[
#   slrcv0_0.5$nclusters_1se_idx]]
# slrcv0_0.5_coefs = getCoefsBM(
#   coefs = coefficients(slrcv0_0.5_selected$model),
#   sbp = slrcv0_0.5_selected$sbp)
# rownames(slrcv0_0.5_coefs$llc.coefs)[
#   slrcv0_0.5_coefs$llc.coefs != 0]
# sum(slrcv0_0.5_coefs$llc.coefs != 0)
# slrcv0_0.5$models[[slrcv0_0.5$nclusters_1se_idx]]$sbp[
#   slrcv0_0.5$models[[slrcv0_0.5$nclusters_1se_idx]]$sbp != 0, ,
#   drop = FALSE]
# 
# ##### cv.hslr with approximation step #####
# # hslrcv0approx_0.5 = cv.hslr(
# #   x = X_0.5, y = Y2, max.levels = slrmax, nfolds = K,
# #   classification = TRUE, approx = TRUE)
# # saveRDS(
# #   hslrcv0approx_0.5,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_hslrcv_approx",
# #     "_0.5",
# #     ".rds"))
# hslrcv0approx_0.5 = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_hslrcv_approx",
#     "_0.5",
#     ".rds"))
# hslrcv0approx_0.5$nclusters_1se_idx
# hslrcv0approx_0.5$nclusters_min_idx
# getTPlots(hslrcv0approx_0.5)
# hslrcv0approx_0.5_selected = hslrcv0approx_0.5$models[[6]]
# hslrcv0approx_0.5_coefs = getCoefsBM(
#   coefs = coefficients(hslrcv0approx_0.5_selected$model),
#   sbp = hslrcv0approx_0.5_selected$sbp)
# rownames(hslrcv0approx_0.5_coefs$llc.coefs)[
#   hslrcv0approx_0.5_coefs$llc.coefs != 0]
# sum(hslrcv0approx_0.5_coefs$llc.coefs != 0)
# 
# hslrcv0approx_0.5_sbps = sapply(hslrcv0approx_0.5$models, function(model) model$sbp)
# rownames(hslrcv0approx_0.5_sbps) = rownames(hslrcv0approx_0.5$models[[1]]$sbp)
# hslrcv0approx_0.5_sbps
# 
# hslrcv0approx_0.5_cors = sapply(hslrcv0approx_0.5$models, function(model) model$cors)
# 
# ##### cv.hslr (no approximation step) #####
# # hslrcv0_0.5 = cv.hslr(
# #   x = X_0.5, y = Y2, max.levels = slrmax, nfolds = K,
# #   classification = TRUE, approx = FALSE)
# # saveRDS(
# #   hslrcv0_0.5,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_hslrcv",
# #     "_0.5",
# #     ".rds"))
# hslrcv0_0.5 = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_hslrcv",
#     "_0.5",
#     ".rds"))
# hslrcv0_0.5$nclusters_1se_idx
# hslrcv0_0.5$nclusters_min_idx
# getTPlots(hslrcv0_0.5)
# hslrcv0_0.5_selected = hslrcv0_0.5$models[[
#   hslrcv0_0.5$nclusters_1se_idx]]
# hslrcv0_0.5_coefs = getCoefsBM(
#   coefs = coefficients(hslrcv0_0.5_selected$model),
#   sbp = hslrcv0_0.5_selected$sbp)
# rownames(hslrcv0_0.5_coefs$llc.coefs)[
#   hslrcv0_0.5_coefs$llc.coefs != 0]
# sum(hslrcv0_0.5_coefs$llc.coefs != 0)

##### selbal #####
# slbl_0.5 = selbal.cv(x = X_0.5, y = Y, n.fold = K)
# saveRDS(
#   slbl_0.5,
#   paste0(
#     output_dir, "/Crohns",
#     "_selbal",
#     "_0.5",
#     ".rds"))
slbl_0.5 = readRDS(
  paste0(
    output_dir, "/Crohns",
    "_selbal",
    "_0.5",
    ".rds"))
slbl_0.5_coefs = getCoefsSelbal(
  X = X_0.5, y = Y, selbal.fit = slbl_0.5, classification = TRUE, check = TRUE)
rownames(slbl_0.5_coefs$llc.coefs)[slbl_0.5_coefs$llc.coefs != 0]
sum(slbl_0.5_coefs$llc.coefs != 0) # 12

##### slr #####
slr0_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = FALSE, amini.regularization = FALSE, highdegree.regularization = FALSE,
  include.leading.eigenvector = FALSE)
slr0_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0_0.5$model), sbp = slr0_0.5$sbp)
rownames(slr0_0.5_coefs$llc.coefs)[slr0_0.5_coefs$llc.coefs != 0]
sum(slr0_0.5_coefs$llc.coefs != 0) # 27

##### slr - amini #####
slr0am_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = FALSE, amini.regularization = TRUE, highdegree.regularization = FALSE,
  include.leading.eigenvector = FALSE)
slr0am_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0am_0.5$model), sbp = slr0am_0.5$sbp)
rownames(slr0am_0.5_coefs$llc.coefs)[slr0am_0.5_coefs$llc.coefs != 0]
sum(slr0am_0.5_coefs$llc.coefs != 0) # 27

##### slr - amini + rank 1 approx #####
slr0amap_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = TRUE, highdegree.regularization = FALSE,
  include.leading.eigenvector = FALSE)
slr0amap_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0amap_0.5$model), sbp = slr0amap_0.5$sbp)
rownames(slr0amap_0.5_coefs$llc.coefs)[slr0amap_0.5_coefs$llc.coefs != 0]
sum(slr0amap_0.5_coefs$llc.coefs != 0) # 22

##### slr - high-degree reg #####
slr0hdr_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = FALSE, amini.regularization = FALSE, highdegree.regularization = TRUE,
  include.leading.eigenvector = FALSE)
slr0hdr_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0hdr_0.5$model), sbp = slr0hdr_0.5$sbp)
rownames(slr0hdr_0.5_coefs$llc.coefs)[slr0hdr_0.5_coefs$llc.coefs != 0]
sum(slr0hdr_0.5_coefs$llc.coefs != 0) # 27

##### slr - high-degree reg + rank 1 approx #####
slr0hdrap_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = FALSE, highdegree.regularization = TRUE,
  include.leading.eigenvector = FALSE)
slr0hdrap_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0hdrap_0.5$model), sbp = slr0hdrap_0.5$sbp)
rownames(slr0hdrap_0.5_coefs$llc.coefs)[slr0hdrap_0.5_coefs$llc.coefs != 0]
sum(slr0hdrap_0.5_coefs$llc.coefs != 0) # 47

##### slr - high-degree reg + amini + rank 1 approx #####
slr0hdramap_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = TRUE, highdegree.regularization = TRUE,
  include.leading.eigenvector = FALSE)
slr0hdramap_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0hdramap_0.5$model), sbp = slr0hdramap_0.5$sbp)
rownames(slr0hdramap_0.5_coefs$llc.coefs)[slr0hdramap_0.5_coefs$llc.coefs != 0]
sum(slr0hdramap_0.5_coefs$llc.coefs != 0) # 22

##### slr - eig12 + high-degree reg + amini + approx #####
slr0all_0.5 = slr(
  x = X_0.5, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = TRUE, highdegree.regularization = TRUE,
  include.leading.eigenvector = TRUE)
slr0all_0.5_coefs = getCoefsBM(
  coefs = coefficients(slr0all_0.5$model), sbp = slr0all_0.5$sbp)
rownames(slr0all_0.5_coefs$llc.coefs)[slr0all_0.5_coefs$llc.coefs != 0]
sum(slr0all_0.5_coefs$llc.coefs != 0) # 22

################################################################################
# Strategy 2: GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

slrcor_gbm = slrmatrix(x = X_gbm, y = Y2)
pheatmap(slrcor_gbm)
# 
# ##### slr with approximation step #####
# slr0approx_gbm = slr(
#   x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, approx = TRUE, 
#   amini.regularization = TRUE, highdegree.regularization = TRUE, 
#   include.leading.eigenvector = TRUE)
# # saveRDS(
# #   slr0approx_gbm,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slr_approx",
# #     "_gbm",
# #     ".rds"))
# # slr0approx_gbm = readRDS(
# #   paste0(
# #     output_dir, "/Crohns", 
# #     "_slr_approx", 
# #     "_gbm", 
# #     ".rds"))
# slr0approx_gbm_coefs = getCoefsBM(
#   coefs = coefficients(slr0approx_gbm$model), sbp = slr0approx_gbm$sbp)
# rownames(slr0approx_gbm_coefs$llc.coefs)[slr0approx_gbm_coefs$llc.coefs != 0]
# sum(slr0approx_gbm_coefs$llc.coefs != 0)
# 
# ##### slr (no approximation step) #####
# # slr0_gbm = slr(
# #   x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, approx = TRUE)
# # saveRDS(
# #   slr0_gbm, 
# #   paste0(
# #     output_dir, "/Crohns", 
# #     "_slr", 
# #     "_gbm", 
# #     ".rds"))
# slr0_gbm = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slr",
#     "_gbm",
#     ".rds"))
# slr0_gbm_coefs = getCoefsBM(
#   coefs = coefficients(slr0_gbm$model), sbp = slr0_gbm$sbp)
# rownames(slr0_gbm_coefs$llc.coefs)[slr0_gbm_coefs$llc.coefs != 0]
# sum(slr0_gbm_coefs$llc.coefs != 0)
# 
# ##### cv.slr with approximation step #####
# # slrcv0approx_gbm = cv.slr(
# #   x = X_gbm, y = Y2, max.clusters = ceiling(sqrt(ncol(X))), nfolds = K,
# #   classification = TRUE, approx = TRUE)
# # saveRDS(
# #   slrcv0approx_gbm,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slrcv_approx",
# #     "_gbm",
# #     ".rds"))
# slrcv0approx_gbm = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slrcv_approx",
#     "_gbm",
#     ".rds"))
# slrcv0approx_gbm$nclusters_1se_idx
# slrcv0approx_gbm$nclusters_min_idx
# getTPlots(slrcv0approx_gbm)
# slrcv0approx_gbm_selected = slrcv0approx_gbm$models[[
#   slrcv0approx_gbm$nclusters_1se_idx]]
# slrcv0approx_gbm_coefs = getCoefsBM(
#   coefs = coefficients(slrcv0approx_gbm_selected$model), 
#   sbp = slrcv0approx_gbm_selected$sbp)
# rownames(slrcv0approx_gbm_coefs$llc.coefs)[
#   slrcv0approx_gbm_coefs$llc.coefs != 0]
# 
# ##### cv.slr (no approximation step) #####
# # slrcv0_gbm = cv.slr(
# #   x = X_gbm, y = Y2, max.clusters = ceiling(sqrt(ncol(X))), nfolds = K,
# #   classification = TRUE, approx = FALSE)
# # saveRDS(
# #   slrcv0_gbm,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_slrcv",
# #     "_gbm",
# #     ".rds"))
# slrcv0_gbm = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_slrcv",
#     "_gbm",
#     ".rds"))
# slrcv0_gbm$nclusters_1se_idx
# slrcv0_gbm$nclusters_min_idx
# getTPlots(slrcv0_gbm)
# slrcv0_gbm_selected = slrcv0_gbm$models[[
#   slrcv0_gbm$nclusters_1se_idx]]
# slrcv0_gbm_coefs = getCoefsBM(
#   coefs = coefficients(slrcv0_gbm_selected$model), 
#   sbp = slrcv0_gbm_selected$sbp)
# rownames(slrcv0_gbm_coefs$llc.coefs)[
#   slrcv0_gbm_coefs$llc.coefs != 0]
# 
# ##### cv.hslr with approximation step #####
# # hslrcv0approx_gbm = cv.hslr(
# #   x = X_gbm, y = Y2, max.levels = ceiling(sqrt(ncol(X))), nfolds = K,
# #   classification = TRUE, approx = TRUE)
# # saveRDS(
# #   hslrcv0approx_gbm,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_hslrcv_approx",
# #     "_gbm",
# #     ".rds"))
# hslrcv0approx_gbm = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_hslrcv_approx",
#     "_gbm",
#     ".rds"))
# hslrcv0approx_gbm$nclusters_1se_idx
# hslrcv0approx_gbm$nclusters_min_idx
# getTPlots(hslrcv0approx_gbm)
# hslrcv0approx_gbm_selected = hslrcv0approx_gbm$models[[
#   hslrcv0approx_gbm$nclusters_1se_idx]]
# hslrcv0approx_gbm_coefs = getCoefsBM(
#   coefs = coefficients(hslrcv0approx_gbm_selected$model), 
#   sbp = hslrcv0approx_gbm_selected$sbp)
# rownames(hslrcv0approx_gbm_coefs$llc.coefs)[
#   hslrcv0approx_gbm_coefs$llc.coefs != 0]
# 
# ##### cv.hslr (no approximation step) #####
# # hslrcv0_gbm = cv.hslr(
# #   x = X_gbm, y = Y2, max.levels = ceiling(sqrt(ncol(X))), nfolds = K,
# #   classification = TRUE, approx = FALSE)
# # saveRDS(
# #   hslrcv0_gbm,
# #   paste0(
# #     output_dir, "/Crohns",
# #     "_hslrcv",
# #     "_gbm",
# #     ".rds"))
# hslrcv0_gbm = readRDS(
#   paste0(
#     output_dir, "/Crohns",
#     "_hslrcv",
#     "_gbm",
#     ".rds"))
# hslrcv0_gbm$nclusters_1se_idx
# hslrcv0_gbm$nclusters_min_idx
# getTPlots(hslrcv0_gbm)
# hslrcv0_gbm_selected = hslrcv0_gbm$models[[
#   hslrcv0_gbm$nclusters_1se_idx]]
# hslrcv0_gbm_coefs = getCoefsBM(
#   coefs = coefficients(hslrcv0_gbm_selected$model), 
#   sbp = hslrcv0_gbm_selected$sbp)
# rownames(hslrcv0_gbm_coefs$llc.coefs)[
#   hslrcv0_gbm_coefs$llc.coefs != 0]
# 
##### selbal #####
# slbl_gbm = selbal.cv(x = X_gbm, y = Y, n.fold = K)
# saveRDS(
#   slbl_gbm,
#   paste0(
#     output_dir, "/Crohns",
#     "_selbal",
#     "_gbm",
#     ".rds"))
slbl_gbm = readRDS(
  paste0(
    output_dir, "/Crohns",
    "_selbal",
    "_gbm",
    ".rds"))
slbl_gbm_coefs = getCoefsSelbal(
  X = X_gbm, y = Y, selbal.fit = slbl_gbm, classification = TRUE, check = TRUE)
rownames(slbl_gbm_coefs$llc.coefs)[slbl_gbm_coefs$llc.coefs != 0]
sum(slbl_gbm_coefs$llc.coefs != 0)

##### slr #####
slr0_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = FALSE, amini.regularization = FALSE, highdegree.regularization = FALSE,
  include.leading.eigenvector = FALSE)
slr0_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0_gbm$model), sbp = slr0_gbm$sbp)
rownames(slr0_gbm_coefs$llc.coefs)[slr0_gbm_coefs$llc.coefs != 0]
sum(slr0_gbm_coefs$llc.coefs != 0) # 25

##### slr - amini #####
slr0am_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = FALSE, amini.regularization = TRUE, highdegree.regularization = FALSE,
  include.leading.eigenvector = FALSE)
slr0am_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0am_gbm$model), sbp = slr0am_gbm$sbp)
rownames(slr0am_gbm_coefs$llc.coefs)[slr0am_gbm_coefs$llc.coefs != 0]
sum(slr0am_gbm_coefs$llc.coefs != 0) # 25

##### slr - amini + rank 1 approx #####
slr0amap_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = TRUE, highdegree.regularization = FALSE,
  include.leading.eigenvector = FALSE)
slr0amap_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0amap_gbm$model), sbp = slr0amap_gbm$sbp)
rownames(slr0amap_gbm_coefs$llc.coefs)[slr0amap_gbm_coefs$llc.coefs != 0]
sum(slr0amap_gbm_coefs$llc.coefs != 0) # 15

##### slr - high-degree reg #####
slr0hdr_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = FALSE, amini.regularization = FALSE, highdegree.regularization = TRUE,
  include.leading.eigenvector = FALSE)
slr0hdr_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0hdr_gbm$model), sbp = slr0hdr_gbm$sbp)
rownames(slr0hdr_gbm_coefs$llc.coefs)[slr0hdr_gbm_coefs$llc.coefs != 0]
sum(slr0hdr_gbm_coefs$llc.coefs != 0) # 25

##### slr - high-degree reg + rank 1 approx #####
slr0hdrap_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = FALSE, highdegree.regularization = TRUE,
  include.leading.eigenvector = FALSE)
slr0hdrap_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0hdrap_gbm$model), sbp = slr0hdrap_gbm$sbp)
rownames(slr0hdrap_gbm_coefs$llc.coefs)[slr0hdrap_gbm_coefs$llc.coefs != 0]
sum(slr0hdrap_gbm_coefs$llc.coefs != 0) # 29

##### slr - high-degree reg + amini + rank 1 approx #####
slr0hdramap_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = TRUE, highdegree.regularization = TRUE,
  include.leading.eigenvector = FALSE)
slr0hdramap_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0hdramap_gbm$model), sbp = slr0hdramap_gbm$sbp)
rownames(slr0hdramap_gbm_coefs$llc.coefs)[slr0hdramap_gbm_coefs$llc.coefs != 0]
sum(slr0hdramap_gbm_coefs$llc.coefs != 0) # 15

##### slr - eig12 + high-degree reg + amini + approx #####
slr0all_gbm = slr(
  x = X_gbm, y = Y2, num.clusters = 2, classification = TRUE, 
  approx = TRUE, amini.regularization = TRUE, highdegree.regularization = TRUE,
  include.leading.eigenvector = TRUE)
slr0all_gbm_coefs = getCoefsBM(
  coefs = coefficients(slr0all_gbm$model), sbp = slr0all_gbm$sbp)
rownames(slr0all_gbm_coefs$llc.coefs)[slr0all_gbm_coefs$llc.coefs != 0]
sum(slr0all_gbm_coefs$llc.coefs != 0) # 15
