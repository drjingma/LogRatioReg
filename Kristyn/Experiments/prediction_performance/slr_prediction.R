getwd()

# libraries
library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
library(limSolve) # for constrained lm, lsei()
library(stats) # for hclust()
library(balance) # for sbp.fromHclust()

# set up parallelization
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(doRNG)
rng.seed = 123
registerDoRNG(rng.seed)

# Dr. Ma sources
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))

# settings
tol = 1e-4
get_lambda = "glmnet"

# Cross-validation
cv.n_lambda = 100
cv.K = 5

# Repetitions
rep.n = 100
n.train = 70
n.test = 28

# data
# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
load(paste0("Data/", "BMI.rda"))
# dim(raw_data) # 98 x 89
# dim(X) # 98 x 87
# dim(X.prop) # 98 x 87
log.X.prop = log(X.prop)
n = dim(X)[1]
num.genera = dim(X)[2]

################################################################################
# Experiments
################################################################################

# They generate 100 bootstrap samples and use the same CV procedure to select 
#   the genera (for stable selection results)
pred.err = foreach(
  b = 1:rep.n, 
  .combine = cbind#, .noexport = c("ConstrLassoC0", "ConstrLaso", "cv.func")
) %dopar% {
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
  library(limSolve) # for constrained lm, lsei()
  library(stats) # for hclust()
  library(balance) # for sbp.fromHclust()
  
  # split into train and test sets
  train.idx = sample(1:n, n.train)
  Xtrain = X.prop[train.idx, ]
  Ytrain = y[train.idx]
  Xtest = X.prop[-train.idx, ]
  Ytest = y[-train.idx]
  
  # refitted CV
  # Split the data into 10 folds
  
  # Fit Lasso on training set
  cv.fits = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K, 
                  get_lambda = get_lambda)
  lambdamin.idx = which.min(cv.fits$cvm)
  betahat = as.matrix(cv.fits$bet[, lambdamin.idx])
  names(betahat) = colnames(betahat)
  # get fitted model
  predictCLM = function(X){
    cv.fits$int[lambdamin.idx] + computeBalances(X, cv.fits$btree) %*% betahat
  }
  # calculate squared error (prediction error?)
  Ypred = predictCLM(Xtest)
  (Ytest - Ypred)^2
}
dim(pred.err)
mse = colMeans(pred.err)
print(paste0("mean prediction error: ", mean(mse)))
print(paste0("standard deviation: ", (sd(mse))))
print(paste0("standard error: ", (sd(mse)) / sqrt(rep.n)))
print(paste0(
  "95% CI: (", 
  mean(mse) - 2 * (sd(mse)) / sqrt(rep.n), 
  ", ",
  mean(mse) + 2 * (sd(mse)) / sqrt(rep.n), ")"
))

saveRDS(pred.err,
        file = paste0("Kristyn/Experiments/output",
                      "/slr_prediction", 
                      "_K", cv.K, 
                      "_seed", rng.seed,
                      ".rds"))
