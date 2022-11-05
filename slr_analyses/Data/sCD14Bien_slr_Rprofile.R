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

# Rprof()
slrspec = slr(
  x = XTr, y = YTr, screen.method = "wald", cluster.method = "spectral",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  threshold = slrspeccv$threshold[slrspeccv$index["1se",]],
  positive.slope = TRUE)
# Rprof(NULL)
# summaryRprof()

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

# Rprof()
slrhier = slr(
  x = XTr, y = YTr, screen.method = "wald", cluster.method = "hierarchical",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  threshold = slrhiercv$threshold[slrhiercv$index["1se",]],
  positive.slope = TRUE)
# Rprof(NULL)
# summaryRprof()

# benchmark! ###################################################################
library(microbenchmark)

x = XTr
y = YTr
screen.method = "wald"
response.type = "continuous"
s0.perc = 0
zeta = 0
threshold = slrspeccv$threshold[slrspeccv$index["1se",]]

# benchmarking getAitchisonVar #################################################

n <- length(y)

feature.scores = getFeatureScores(x, y, screen.method, response.type, s0.perc)
which.features <- (abs(feature.scores) >= threshold)
x.reduced <- x[,which.features, drop = FALSE] # reduced data matrix

getAitchisonVar = function(x){
  p = ncol(x)
  A <- matrix(0, p, p)
  for (j in 1:p){
    for (k in 1:p){
      if (k == j){
        next
      }
      else{
        A[j,k] <- stats::var(log(x[,j])-log(x[,k])) # Aitchison variation
      }
    }
  }
  return(A)
}

Aitchison.var = getAitchisonVar(x.reduced)

# try using outer() instead
AitchVar = function(x, y) stats::var(log(x) - log(y))
AitchVarVec = Vectorize(AitchVar)
Aitchison.var2 = outer(X = x.reduced, Y = x.reduced, FUN = AitchVarVec)
rownames(Aitchison.var2) <- colnames(Aitchison.var2) <- NULL

all.equal(Aitchison.var, Aitchison.var2, check.attributes = FALSE)

Outer <- function(x,y,fun) {
  mat <- matrix(mapply(fun, rep(x, length(y)), 
                       rep(y, each=length(x))),
                length(x), length(y))
}

# using outer() is faster than using a double for loop!
microbenchmark(
  getAitchisonVar(x.reduced), 
  outer(X = x.reduced, Y = x.reduced, FUN = AitchVarVec),
  Outer(x = x.reduced, y = x.reduced, fun = AitchVar),
  times = 1000
)

# benchmarking slr.fromContrast ################################################

rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(x.reduced)
Aitchison.sim <- max(Aitchison.var) - Aitchison.var
sbp.est <- spectral.clustering(Aitchison.sim, zeta = zeta)
cluster.mat = Aitchison.sim

slr.fromContrast <- function(x, contrast){
  
  if(length(contrast) != ncol(x)) stop("Contrast must have length ncol(x) = D.")
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")
  
  # lpos <- sum(contrast == 1)
  # lneg <- sum(contrast == -1)
  # const <- sqrt((lpos*lneg)/(lpos+lneg))
  logX <- log(x)
  ipos <- rowMeans(logX[, contrast == 1, drop = FALSE])
  ineg <- rowMeans(logX[, contrast == -1, drop = FALSE])
  
  ipos - ineg
}

balance <- slr.fromContrast(x.reduced, sbp.est) 

# try using apply instead

slr.fromContrast2 <- function(x, contrast){
  if(length(contrast) != ncol(x)) stop("Contrast must have length ncol(x) = D.")
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")
  logX <- log(x)
  ipos <- apply(logX[, contrast == 1, drop = FALSE], 1, mean)
  ineg <- apply(logX[, contrast == -1, drop = FALSE], 1, mean)
  ipos - ineg
}

balance2 <- slr.fromContrast2(x.reduced, sbp.est) 
all.equal(balance, balance2)

# rowMeans is actually faster in this case -- keep it.
microbenchmark(
  slr.fromContrast(x.reduced, sbp.est),
  slr.fromContrast2(x.reduced, sbp.est),
  times = 100
)

# check spectral clustering ####################################################
#   is ordering necessary?

spectral.clustering <- function(W, n_eig = 2, zeta = 0) {
  L = graph.laplacian(W,zeta = zeta) # compute graph Laplacian
  ei = eigen(L, symmetric = TRUE)    # Compute the eigenvectors and values of L
  # we will use k-means to cluster the eigenvectors corresponding to
  # the leading smallest eigenvalues
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  obj <- kmeans(ei$vectors[, 1:n_eig], centers = n_eig, nstart = 100, algorithm = "Lloyd")
  if (n_eig==2){
    cl <- 2*(obj$cluster - 1) - 1 
  } else {
    cl <- obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl) 
}

sbp.est <- spectral.clustering(Aitchison.sim, zeta = zeta)

spectral.clustering2 <- function(W, n_eig = 2, zeta = 0) {
  L = graph.laplacian(W,zeta = zeta) # compute graph Laplacian
  ei = eigen(L, symmetric = TRUE)    # Compute the eigenvectors and values of L
  # we will use k-means to cluster the eigenvectors corresponding to
  # the leading smallest eigenvalues
  # ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  obj <- kmeans(ei$vectors[, 1:n_eig], centers = n_eig, nstart = 100, algorithm = "Lloyd")
  if (n_eig==2){
    cl <- 2*(obj$cluster - 1) - 1 
  } else {
    cl <- obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl) 
}

sbp.est2 <- spectral.clustering2(Aitchison.sim, zeta = zeta)

# they are the same
sbp.est
-sbp.est2

# somehow, max is much larger for spectral.clustering2 -- might as well keep it.
microbenchmark(
  spectral.clustering(Aitchison.sim, zeta = zeta),
  spectral.clustering2(Aitchison.sim, zeta = zeta),
  times = 1000
)

# benchmarking spectral clustering functions ###################################

# spectral clustering methods #

# matrix to be spectrally-clustered

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




