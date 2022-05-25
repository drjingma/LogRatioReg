# Purpose: compare slr to selbal on data sets
# Date: 5/25/2022
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
# Strategy 2: 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

################################################################################
# fit methods
################################################################################

# classo #######################################################################
cl_gbm = codalasso(X_gbm, Y2, numFolds = K)
saveRDS(
  cl_gbm,
  paste0(
    output_dir, "/Crohns",
    "_classo",
    "_gbm",
    ".rds"))

# slr -- alpha = 0.01 ##########################################################
slr0.01 = slr(x = X_gbm, y = Y2, alpha = 0.01, classification = TRUE)
saveRDS(
  slr0.01,
  paste0(
    output_dir, "/Crohns",
    "_slr_alpha0.01",
    "_gbm",
    ".rds"))

# slr -- alpha = 0.05 ##########################################################
slr0.05 = slr(x = X_gbm, y = Y2, alpha = 0.05, classification = TRUE)
saveRDS(
  slr0.05,
  paste0(
    output_dir, "/Crohns",
    "_slr_alpha0.05",
    "_gbm",
    ".rds"))

# slr -- screen ################################################################
slrscreen0cv = cv.slr.screen(
  x = X_gbm, y = Y2, method = "wald", 
  response.type = "binary", s0.perc = 0, zeta = 0, 
  nfolds = K, type.measure = "mse", 
  parallel = FALSE, scale = scaling, trace.it = FALSE)
slrscreen0 = slr.screen(
  x = X_gbm, y = Y2, method = "wald", 
  response.type = "binary", s0.perc = 0, zeta = 0, 
  threshold = slrscreen0cv$threshold[slrscreen0cv$index["1se",]])
saveRDS(
  slrscreen0,
  paste0(
    output_dir, "/Crohns",
    "_slrscreen",
    "_gbm",
    ".rds"))

# selbal #######################################################################
slbl_gbm = selbal.cv(x = X_gbm, y = Y, n.fold = K)
saveRDS(
  slbl_gbm,
  paste0(
    output_dir, "/Crohns",
    "_selbal",
    "_gbm",
    ".rds"))

# codacore #####################################################################
library(codacore)
codacore0 = codacore(
  x = X_gbm, y = Y2, logRatioType = "ILR", 
  objective = "binary classification", cvParams = list(numFolds = K))
saveRDS(
  codacore0,
  paste0(
    output_dir, "/Crohns",
    "_codacore",
    "_gbm",
    ".rds"))










