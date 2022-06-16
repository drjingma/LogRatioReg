rm(list=ls())
# Purpose: compare methods on prediction accuracy & timing using Crohns data
# Date: 6/13/2022

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs_mse"

source("slr_analyses/Functions/util.R")

library(tidyverse)
library(reshape2)

numSplits = 20

# tuning parameter settings
K = 10
scaling = TRUE

file.end0 = paste0(
  # "_sim", b,
  "_Crohns", 
  "_gbm")

################################################################################
# find the data split that results in low AUC for slr

bad.split.idx = NA
for(i in 1:numSplits){
  
  # slr
  slrscreen_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(slrscreen_tmp) = NULL
  if(data.table::data.table(slrscreen_tmp)$auc < 0.65){
    print(paste0("sim index that results in low AUC for slr is: ", i, " ###"))
    bad.split.idx = i
  }
  
}

################################################################################
# import all data/methods for that data split to compare them

# data
data.tmp = readRDS(paste0(
  output_dir, "/data", "_sim", bad.split.idx, file.end0, ".rds"
))

# compositional lasso
cl.tmp = readRDS(paste0(
  output_dir, "/classo_metrics",
  "_sim", bad.split.idx, file.end0, ".rds"
))
# slr
slr.tmp = readRDS(paste0(
  output_dir, "/slr_metrics",
  "_sim", bad.split.idx, file.end0, ".rds"
))
# selbal
slbl.tmp = readRDS(paste0(
  output_dir, "/selbal_metrics",
  "_sim", bad.split.idx, file.end0, ".rds"
))
# codacore
cdcr.tmp = readRDS(paste0(
  output_dir, "/codacore_metrics",
  "_sim", bad.split.idx, file.end0, ".rds"
))
metrics.tmp = t(data.frame(
  classo = cl.tmp, 
  slr = slr.tmp, 
  selbal = slbl.tmp, 
  codacore = cdcr.tmp
))

################################################################################
# see how many (and which) variables were selected by each method

XTr = data.tmp$XTr
YTr = data.tmp$YTr
Y2Tr = data.tmp$Y2Tr
XTe = data.tmp$XTe
YTe = data.tmp$YTe
Y2Te = data.tmp$Y2Te

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(selbal)

library(pROC)

source("RCode/func_libs.R")
source("RCode/SLR.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

# classo #####################################################################
start.time = Sys.time()
classo = codalasso(XTr, Y2Tr, numFolds = K)
end.time = Sys.time()
cl.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

# get prediction error on test set
classo.Yhat.test = predict(classo, XTe)

cl.metrics = c(
  acc = mean((classo.Yhat.test > 0) == Y2Te),
  auc = roc(
    Y2Te, classo.Yhat.test, levels = c(0, 1), direction = "<")$auc,
  time = cl.timing
)

# slr ########################################################################
start.time = Sys.time()
slrscreen0cv = cv.slr(
  x = XTr, y = Y2Tr, method = "wald", 
  response.type = "binary", s0.perc = 0, zeta = 0, 
  nfolds = K, type.measure = "mse", 
  parallel = FALSE, scale = scaling, trace.it = FALSE)
slrscreen0 = slr(
  x = XTr, y = Y2Tr, method = "wald", 
  response.type = "binary", s0.perc = 0, zeta = 0, 
  threshold = slrscreen0cv$threshold[slrscreen0cv$index["1se",]])
end.time = Sys.time()
slrscreen0.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

slrscreen0.fullSBP = matrix(0, nrow = ncol(XTr), ncol = 1)
rownames(slrscreen0.fullSBP) = colnames(XTr)
slrscreen0.fullSBP[match(
  names(slrscreen0$sbp), rownames(slrscreen0.fullSBP))] = slrscreen0$sbp

# get prediction error on test set
slrscreen0.Yhat.test = predict(
  slrscreen0$fit, 
  data.frame(balance = balance::balance.fromSBP(
    x = XTe, y = slrscreen0.fullSBP)), 
  type = "response")

slrscreen0.metrics = c(
  acc = mean((slrscreen0.Yhat.test > 0.5) == Y2Te),
  auc = roc(
    Y2Te, slrscreen0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
  time = slrscreen0.timing
)

# selbal #####################################################################
start.time = Sys.time()
slbl = selbal.cv(x = XTr, y = YTr, n.fold = K)
end.time = Sys.time()
slbl.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

slbl.coefs = getCoefsSelbal(
  X = XTr, y = YTr, selbal.fit = slbl, classification = TRUE, 
  check = TRUE)

# get prediction error on test set
slbl.Yhat.test = predict.glm(
  slbl$glm, 
  newdata = data.frame(V1 = balance::balance.fromSBP(
    x = XTe, y = slbl.coefs$sbp)), 
  type = "response")

slbl.metrics = c(
  acc = mean((slbl.Yhat.test < 0.5) == Y2Te),
  auc = roc(
    YTe, slbl.Yhat.test, levels = levels(YTe), direction = "<")$auc,
  time = slbl.timing
)

# codacore ###################################################################
library(codacore)
start.time = Sys.time()
codacore0 = codacore(
  x = XTr, y = Y2Tr, logRatioType = "ILR", 
  objective = "binary classification", cvParams = list(numFolds = K))
end.time = Sys.time()
codacore0.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

if(length(codacore0$ensemble) > 0){ # at least 1 log-ratio found
  # get prediction error on test set
  codacore0.Yhat.test = predict(codacore0, XTe)
  
} else{
  print(paste0("sim ", i, " -- codacore has no log-ratios"))
  codacore0model = stats::glm(Y ~ 1, family = "binomial")
  
  # get prediction error on test set
  codacore0.Yhat.test = predict(codacore0model, data.frame(XTe))
}

codacore0.metrics = c(
  acc = mean((codacore0.Yhat.test > 0.5) == Y2Te),
  auc = roc(
    Y2Te, codacore0.Yhat.test, levels = c(0, 1), direction = "<")$auc,
  time = codacore0.timing
)


