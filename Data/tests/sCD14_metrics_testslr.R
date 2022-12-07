# Purpose: test slr on sCD14 data set
# Date: 8/30/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "Data/outputs_metrics"

# Other simulation settings
numSplits = 20

################################################################################
# 70/30 train/test splits with sparsity filtering at 80%
################################################################################

b = 4
library(mvtnorm)

library(Matrix)
library(glmnet)

library(pROC)

source("Functions/slrs.R")
source("Functions/util.R")

getF1 = function(y, yhat){
  TPR = sum((yhat > 0) & (y == 1)) / sum(y == 1) # TPR (recall)
  prec = sum((yhat > 0) & (y == 1)) / sum(yhat > 0) # precision
  F1 = 2 / (1/TPR + 1/prec) # f1 is harmonic mean of precision & recall
  return(F1)
}

# tuning parameter settings
hparam = "1se"
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4

filter.perc = 0.8 # 0.8, 1
split.perc = 0.7 # 0.7, 0.8

file.end = paste0(
  "_sCD14", 
  "_split", split.perc, 
  "_filter", filter.perc, 
  "_hparam", hparam,
  "_gbm",
  "_sim", b,
  ".rds")

##############################################################################
# sCD14: another HIV data set in selbal package
#   n = 151 samples (a subset from sCD14 data set), 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
W = selbal::sCD14[, 1:60]
W.origin <- W
W <- W[,apply(W==0,2,mean)<filter.perc]
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::sCD14[, 61]
p = ncol(W)

##############################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = selbal::cmultRepl2(W, zero.rep = "bayes")

##############################################################################
# Train/Test Split
#   Following Gordon-Rodriguez et al. 2022, fit each method on random
#     train/test splits, 
#     -- since response is continuous, no need to stratify by case-control.
numObs = nrow(X_gbm)
inputDim = ncol(X_gbm)
data.tmp = readRDS(paste0(output_dir, "/data", file.end))
XTr = data.tmp$XTr
XTe = data.tmp$XTe
YTr = data.tmp$YTr
YTe = data.tmp$YTe

set.seed(1)
# slr - spectral #############################################################
start.time = Sys.time()
slrspeccv = cv.slr(
  x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  nfolds = K, type.measure = "mse",
  scale = scaling, trace.it = FALSE)
if(hparam == "min"){
  slrspec = slr(
    x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0,
    threshold = slrspeccv$threshold[slrspeccv$index["min",]],
    positive.slope = TRUE)
} else if(hparam == "1se"){
  slrspec = slr(
    x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0,
    threshold = slrspeccv$threshold[slrspeccv$index["1se",]],
    positive.slope = TRUE)
} else{
  stop("invalid hparam setting (method for selecting hyperparameter(s)).")
}
end.time = Sys.time()
slrspec.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

# get SBP
slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrspec.fullSBP) = colnames(X)
slrspec.fullSBP[match(
  names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp

slrspec$sbp
# why is g_Collinsella +?
pheatmap::pheatmap(slrspec$cluster.mat) # there aren't 2 clusters (more like 3)

# lines of slr()
x = XTr
y = YTr
screen.method = "wald"
cluster.method = "spectral"
response.type = "continuous"
s0.perc = 0
zeta = 0
threshold = slrspeccv$threshold[slrspeccv$index["1se",]]
positive.slope = TRUE

if(!("data.frame" %in% class(x))) x = data.frame(x)

n <- length(y)

feature.scores = getFeatureScores(x, y, screen.method, response.type, s0.perc)
which.features <- (abs(feature.scores) >= threshold)

x.reduced <- x[,which.features] # reduced data matrix
Aitchison.var = getAitchisonVar(x.reduced)
rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(x.reduced)

Aitchison.sim <- max(Aitchison.var) - Aitchison.var

## Perform spectral clustering
sbp.est <- spectral.clustering(Aitchison.sim, zeta = zeta)
sbp.est
graph <- learn_k_component_graph(Aitchison.sim, k = 2, beta = 1, verbose = FALSE, abstol = 1e-8)
graph
