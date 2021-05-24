getwd()

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)
library(balance) # for sbp.fromHclust()
library(rare)

# plotting libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)
library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Dr. Ma sources
library(Matrix)
library(glmnet)
library(compositions)
library(stats)
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "supervisedlogratiosalpha.R"))

# Method Settings
tol = 1e-4
nlam = 10# 100
intercept = TRUE
K = 3 #10
rho.type = "square"
linkage = "average"

# Simulation settings
# numSims = 100
n = 6 # 100
p = 8 # 200
rho = 0.5 # 0.2, 0.5
# which beta?
beta.settings = "new"
if(beta.settings == "old" | beta.settings == "linetal2014"){
  beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
} else{
  beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
}
non0.beta = (beta != 0)
sigma_eps = 0.5
seed = 1
muW = c(
  rep(log(p), 5), 
  rep(0, p - 5)
)
SigmaW = matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    SigmaW[i, j] = rho^abs(i - j)
  }
}

################################################################################
# Simulated Data #
################################################################################
set.seed(123)
# set.seed(1)
# set.seed(16) # slr0 is more sparse

# simulate training data #
# generate W
W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
# let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
V = exp(W)
rowsumsV = apply(V, 1, sum)
X = V / rowsumsV
epsilon = rnorm(n, 0, sigma_eps)
Z = log(X)
# generate Y
Y = Z %*% beta + epsilon

# simulate test data #
# simulate independent test set of size n
# generate W
W.test = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
# let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
V.test = exp(W.test)
rowsumsV.test = apply(V.test, 1, sum)
X.test = V.test / rowsumsV.test
epsilon.test = rnorm(n, 0, sigma_eps)
Z.test = log(X.test)
# generate Y
Y.test = Z.test %*% beta + epsilon.test

# get foldid like cv.glmnet does it
#   https://github.com/cran/glmnet/blob/master/R/cv.glmnet.R
nfolds = K
foldid = sample(rep(seq(nfolds), length = n))

################################################################################
# Applying SLR #
################################################################################

# apply the old slr #
slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
            rho.type = rho.type, linkage = linkage, foldid = foldid, 
            standardize = TRUE)

# choose lambda
lam.min.idx0 = which.min(slr$cvm)
lam.min0 = slr$lambda[lam.min.idx0]

# plot cvm vs lambda
# plot(slr$lambda, slr$cvm, #xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
#      col = rgb(0,0,0, alpha = 0.25))
ggplot(data.frame(lambda = slr$lambda, cvm = slr$cvm), 
       aes(x = lambda, y = cvm)) + 
  geom_path() + 
  geom_point(alpha = 0.25) + 
  theme_classic()
# what is the minimizing lambda?
lam.min0

# betahat?
thetahat0 = slr$bet[, lam.min.idx0]
betahat0 = getBeta(thetahat0, btree = slr$btree)
sum(betahat0 != 0)

# old slr, with some changes
slr.v2 = cvSLR0(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
                rho.type = rho.type, linkage = linkage, foldid = foldid, 
                scaling = TRUE)
all.equal(slr$lambda, slr.v2$lambda)
all.equal(slr$bet, slr.v2$bet)
all.equal(slr$a0, slr.v2$a0)
all.equal(slr$cvm, as.numeric(slr.v2$cvm)) # not the same!
cbind(slr$cvm, as.numeric(slr.v2$cvm))
ggplot(data.frame(lambda = c(slr$lambda, slr.v2$lambda), 
                  cvm = c(slr$cvm, as.numeric(slr.v2$cvm)), 
                  type = c(rep("SLR cv.glmnet()", length(slr$lambda)),
                           rep("SLR manual CV", length(slr.v2$lambda))
                           )
                  ),
       aes(x = lambda, y = cvm, color = type, linetype = type)) + 
  geom_path(size = 1) + 
  geom_point(alpha = 0.25, size = 2) + 
  theme_classic()



################################################################################
# line-by-line
y = Y; nfolds = K; lambda = NULL

################################################################################
# calculat MSE by hand...

# the large p is, the higher the discrepancy between the ROC curves...
all.equal(slr$int, slr.v2$int) # TRUE
all.equal(slr$bet, slr.v2$bet) # TRUE
all.equal(slr$lambda, slr.v2$lambda) # TRUE
all.equal(as.vector(slr$cvm), as.vector(slr.v2$cvm)) # "Mean relative difference: 0.007774666"
as.vector(slr$cvm) # uses cv.glmnet
as.vector(slr.v2$cvm) # calculates by hand the same way as slralpha

# calculating by hand a different way:
folds = list()
for(i in 1:nfolds) folds[[i]] = which(foldid == i)

lambdas = slr$lambda
fit.preval = matrix(NA, nrow = n, ncol = nlam)

# fold 1 #
i = 1
# fit the model, ommitting the fold
fold1.fit = glmnet(x = X[-folds[[i]], ], y = Y[-folds[[i]]], 
                   lambda = lambdas, nlambda = nlam, 
                   intercept = intercept)
# predicted values on the fold
fold1.pred = X[folds[[i]], ] %*% fold1.fit$beta + 
  rep(fold1.fit$a0, each = length(folds[[i]]))
# error = squared deviation from the predicted value and true y
fold1.err = fold1.pred
for(j in 1:nlam) fold1.err[, j] = (fold1.err[, j] - Y[folds[[i]]])^2
# fold1.err
# Y[folds[[i]]]
# fold1.err - Y[folds[[i]]]
# (fold1.err - Y[folds[[i]]])^2
fold1.errmeans = colMeans(fold1.err)
slr.v2$errs[, , i] # matches slr.v2
# slr$cv.glmnet$fit.preval
# what is fit.preval?
# a prevalidated array is returned containing fitted values for each observation 
# and each value of ‘lambda’. This means these fits are computed with this 
# observation and the rest of its fold omitted
for(j in 1:length(folds[[i]])){
  fit.preval[folds[[i]][j], ] = fold1.pred[j, ]
}
# predicted values on all data
fold1.pred.all = X %*% fold1.fit$beta + rep(fold1.fit$a0, each = n)
# error = squared deviation from the predicted value and true y
fold1.err.all = fold1.pred.all
for(j in 1:nlam) fold1.pred.all[, j] = (fold1.pred.all[, j] - Y)^2
fold1.err.all

# fold 2 #
i = 2
# fit the model, ommitting the fold
fold2.fit = glmnet(x = X[-folds[[i]], ], y = Y[-folds[[i]]], 
                   lambda = lambdas, nlambda = nlam, 
                   intercept = intercept)
# predicted values on the fold
fold2.pred = X[folds[[i]], ] %*% fold2.fit$beta + 
  rep(fold2.fit$a0, each = length(folds[[i]]))
# error = squared deviation from the predicted value and true y
fold2.err = fold2.pred
for(j in 1:nlam) fold2.err[, j] = (fold2.err[, j] - Y[folds[[i]]])^2
fold2.errmeans = colMeans(fold2.err)
slr.v2$errs[, , i] # matches slr.v2
for(j in 1:length(folds[[i]])){
  fit.preval[folds[[i]][j], ] = fold2.pred[j, ]
}
# predicted values on all data
fold2.pred.all = X %*% fold2.fit$beta + rep(fold2.fit$a0, each = n)
# error = squared deviation from the predicted value and true y
fold2.err.all = fold2.pred.all
for(j in 1:nlam) fold2.pred.all[, j] = (fold2.pred.all[, j] - Y)^2

# fold 3 #
i = 3
# fit the model, ommitting the fold
fold3.fit = glmnet(x = X[-folds[[i]], ], y = Y[-folds[[i]]], 
                   lambda = lambdas, nlambda = nlam, 
                   intercept = intercept)
# predicted values on the fold
fold3.pred = X[folds[[i]], ] %*% fold3.fit$beta + 
  rep(fold3.fit$a0, each = length(folds[[i]]))
# error = squared deviation from the predicted value and true y
fold3.err = fold3.pred
for(j in 1:nlam) fold3.err[, j] = (fold3.err[, j] - Y[folds[[i]]])^2

fold3.errmeans = colMeans(fold3.err)
slr.v2$errs[, , i] # matches slr.v2
for(j in 1:length(folds[[i]])){
  fit.preval[folds[[i]][j], ] = fold3.pred[j, ]
}
# predicted values on all data
fold3.pred.all = X %*% fold3.fit$beta + rep(fold3.fit$a0, each = n)
# error = squared deviation from the predicted value and true y
fold3.err.all = fold3.pred.all
for(j in 1:nlam) fold3.pred.all[, j] = (fold3.pred.all[, j] - Y)^2
pred.all.foldmean = (fold1.pred.all + fold2.pred.all + fold3.pred.all) / nfolds
all.equal(as.matrix(slr$cv.glmnet$fit.preval), fit.preval)

mean((slr$cv.glmnet$fit.preval[, nlam] - Y)^2)
slr$cv.glmnet$cvm[nlam]

# turns out the problem was that I was fitting X instead of Xb at each fold!

# scaling issues... posssibly centering issues too?

btree = slr$btree
U = getU(btree = btree)
logX = log(X)
logX.cent = logX - matrix(rep(1, times = n), ncol = 1) %*% 
  matrix(colMeans(logX), nrow = 1)
logX.cent.U = logX.cent %*% U
Xb = computeBalances(X, btree)
ilrX.cent = Xb - matrix(rep(1, times = n), ncol = 1) %*% 
  matrix(colMeans(Xb), nrow = 1)
all.equal(logX.cent.U, ilrX.cent) # TRUE

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# apply the new slr, with alpha = 1 -- should be equivalent to old slr #
slr.a1 = cvSLRalpha(
  y = Y, X = X, linkage = linkage, rho.type = rho.type, 
  intercept = intercept, lambda = slr$lambda, 
  alpha = 1, nlam = nlam, nfolds = K, foldid = foldid, scaling = TRUE
)

# choose lambda
lam.min.idx = which.min(slr.a1$cvm)
lam.min = slr.a1$lambda[lam.min.idx]

# plot cvm vs lambda
# plot(slr.a1$lambda, slr.a1$cvm, col = rgb(0,0,0, alpha = 0.25))
ggplot(data.frame(lambda = slr.a1$lambda, cvm = slr.a1$cvm), 
       aes(x = lambda, y = cvm)) + 
  geom_path() + 
  geom_point(alpha = 0.25) + 
  theme_classic()
# what is the minimizing lambda?
lam.min


################################################


# plot both together
ggplot(data.frame(lambda = c(slr$lambda, slr.v2$lambda, slr.a1$lambda), 
                  cvm = c(slr$cvm, slr.v2$cvm, slr.a1$cvm), 
                  type = c(rep("SLR, cv.glmnet", length(slr$lambda)),
                           rep("SLR, manual cvm", length(slr.v2$lambda)), 
                           rep("SLRalpha1", length(slr.a1$lambda)))),
       aes(x = lambda, y = cvm, color = type, linetype = type)) + 
  geom_path(size = 1) + 
  geom_point(alpha = 0.25, size = 2) + 
  theme_classic()

# thetahat?
slr$bet
slr.v2$bet
slr.a1$theta
all.equal(unname(as.matrix(slr$bet)), unname(as.matrix(slr.v2$bet)))
all.equal(unname(as.matrix(slr$bet)), unname(as.matrix(slr.a1$theta)))

#betahat?
U %*% slr$bet
U %*% slr.v2$bet
slr.a1$beta

betahat = slr$bet[[1]][, lam.min.idx]
thetahat = slr$thet[[1]][, lam.min.idx]
sum(betahat != 0)

# cbind(betahat0, betahat)
# sum(betahat0 != 0)
# sum(betahat != 0)

info = data.frame(
  "lambda min" = c(lam.min0, lam.min), 
  "min cvm" = c(min(slr0$cvm), min(slr$cvm)), 
  "S.hat" = c(sum(betahat0 != 0), sum(betahat != 0)), 
  "S0" = rep(sum(beta != 0), 2))
rownames(info) = c("SLR", "SLR2")
info
# library(knitr)
# kable(info, "latex")




