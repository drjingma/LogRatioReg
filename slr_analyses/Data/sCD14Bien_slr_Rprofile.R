# Purpose: use sCD14 data set in Bien et al., 2020 (trac)
# Date: 9/13/2022
rm(list=ls())

################################################################################
# libraries and settings

library(mvtnorm)

library(Matrix)
library(glmnet)

library(pROC)

source("RCode/func_libs.R")
source("slr_analyses/Functions/slrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

getF1 = function(y, yhat){
  TPR = sum((yhat > 0) & (y == 1)) / sum(y == 1) # TPR (recall)
  prec = sum((yhat > 0) & (y == 1)) / sum(yhat > 0) # precision
  F1 = 2 / (1/TPR + 1/prec) # f1 is harmonic mean of precision & recall
  return(F1)
}

# tuning parameter settings
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4

################################################################################
# sCD14: another HIV data set, this time given by Bien 2020 (trac paper)
#   n = 152 samples (a subset from sCD14 data set), 
#   p = 539 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
sCD14 = readRDS("Data/sCD14.RDS")
W.orig = sCD14$x
Y = sCD14$y

# filtering
W.has = W.orig != 0
W.prop.samples = apply(W.has, 2, mean)

W = W.orig[, W.prop.samples >= 0.25]
X = sweep(W, 1, rowSums(W), FUN='/')
p = ncol(W)

##############################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = selbal::cmultRepl2(W, zero.rep = "bayes")

##############################################################################
# Train/Test Split
#   Following Gordon-Rodriguez et al. 2022, fit each method on 20 random 80/20
#     train/test splits, 
#     -- since response is continuous, no need to stratify by case-control.

numObs = nrow(X_gbm)
inputDim = ncol(X_gbm)
trainIdx = sample(cut(1:numObs, breaks=5, labels=F))
XTr = X_gbm[trainIdx != 1,]
YTr = Y[trainIdx != 1]
XTe = X_gbm[trainIdx == 1,]
YTe = Y[trainIdx == 1]

##############################################################################
# fit methods
##############################################################################

# slr - spectral #############################################################
# slrspeccv = cv.slr(
#   x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "mse",
#   scale = scaling, trace.it = FALSE)
# saveRDS(
#   slrspeccv, 
#   "slr_analyses/Data/slrspeccv1.rds"
# )
slrspeccv = readRDS("slr_analyses/Data/slrspeccv1.rds")

Rprof()
slrspec = slr(
  x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  threshold = slrspeccv$threshold[slrspeccv$index["1se",]],
  positive.slope = TRUE)
Rprof(NULL)
summaryRprof()

# slr - hierarchical #########################################################
# slrhiercv = cv.slr(
#   x = XTr, y = YTr, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "mse",
#   scale = scaling, trace.it = FALSE)
# saveRDS(
#   slrspeccv, 
#   "slr_analyses/Data/slrhiercv1.rds"
# )
slrhiercv = readRDS("slr_analyses/Data/slrhiercv1.rds")

Rprof()
slrhier = slr(
  x = XTr, y = YTr, screen.method = "wald", cluster.method = "hierarchical",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  threshold = slrhiercv$threshold[slrhiercv$index["1se",]],
  positive.slope = TRUE)
Rprof(NULL)
summaryRprof()

# benchmarking some functions ##################################################

# spectral clustering methods #

# matrix to be spectrally-clustered
x = XTr
y = YTr
screen.method = "wald"
response.type = "continuous"
s0.perc = 0
zeta = 0
threshold = slrspeccv$threshold[slrspeccv$index["1se",]]
n <- length(y)
feature.scores = getFeatureScores(x, y, screen.method, response.type, s0.perc)
which.features <- (abs(feature.scores) >= threshold)
x.reduced <- x[,which.features] # reduced data matrix
Aitchison.var = getAitchisonVar(x.reduced)
rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(x.reduced)
Aitchison.sim <- max(Aitchison.var) - Aitchison.var 

# our spectral clustering
sbp.est <- spectral.clustering(Aitchison.sim, zeta = zeta)
sbp.est

# sClust
library(sClust)
fastClustering(Aitchison.sim, smplPoint = 3)

# kernlab
library(kernlab)
kernlab.sc <- specc(x = Aitchison.sim, centers = 2, kernel = "laplacedot")
plot(my.data, col=sc, pch=4)            # estimated classes (x)
points(my.data, col=obj$classes, pch=5) # true classes (<>)




