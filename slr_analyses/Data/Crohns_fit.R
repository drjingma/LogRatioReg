# Purpose: compare slr to other methods on data sets
# Date: 6/16/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "Kristyn/Data/outputs"

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(selbal)

library(pROC)

# library(zCompositions)
library(pheatmap)

library(tidyverse)

source("RCode/func_libs.R")
source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/slrscreen.R")
source("Kristyn/Functions/codalasso.R")
source("Kristyn/Functions/util.R")

# tuning parameter settings
K = 10
scaling = TRUE

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
#     output_dir, "/Crohns",
#     "_classo",
#     "_gbm",
#     ".rds"))
cl_gbm = readRDS(
  paste0(
    output_dir, "/Crohns",
    "_classo",
    "_gbm",
    ".rds"))

# slr ##########################################################################
# slrscreen0cv = cv.slr.screen(
#   x = X_gbm, y = Y2, method = "wald", 
#   response.type = "binary", s0.perc = 0, zeta = 0, 
#   nfolds = K, type.measure = "mse", 
#   parallel = FALSE, scale = scaling, trace.it = FALSE)
# slrscreen0 = slr.screen(
#   x = X_gbm, y = Y2, method = "wald", 
#   response.type = "binary", s0.perc = 0, zeta = 0, 
#   threshold = slrscreen0cv$threshold[slrscreen0cv$index["1se",]])
# saveRDS(
#   slrscreen0,
#   paste0(
#     output_dir, "/Crohns",
#     "_slrscreen",
#     "_gbm",
#     ".rds"))
slr0 = readRDS(
  paste0(
    output_dir, "/Crohns",
    "_slrscreen",
    "_gbm",
    ".rds"))

# selbal #######################################################################
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

# codacore #####################################################################
# library(codacore)
# codacore0 = codacore(
#   x = X_gbm, y = Y2, logRatioType = "ILR", 
#   objective = "binary classification", cvParams = list(numFolds = K))
# saveRDS(
#   codacore0,
#   paste0(
#     output_dir, "/Crohns",
#     "_codacore",
#     "_gbm",
#     ".rds"))
codacore0 = readRDS(
  paste0(
    output_dir, "/Crohns",
    "_codacore",
    "_gbm",
    ".rds"))









