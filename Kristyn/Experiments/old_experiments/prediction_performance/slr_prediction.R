getwd()
output_dir = "Kristyn/Experiments/prediction_performance/output"

# libraries
library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
library(limSolve) # for constrained lm, lsei()
library(stats) # for hclust()
library(balance) # for sbp.fromHclust()

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

# settings
tol = 1e-4
cv.n_lambda = 200
cv.K = 5
intercept = TRUE

# Repetitions
numReps = 100
n.train = 70
n.test = 28

# data
load(paste0("Data/", "BMI.rda"))
log.X.prop = log(X.prop)
n = dim(X)[1]
num.genera = dim(X)[2]

################################################################################
# Experiments
################################################################################

# They generate 100 bootstrap samples and use the same CV procedure to select 
#   the genera (for stable selection results)
pred.err = foreach(
  b = 1:numReps, 
  .combine = cbind
) %dorng% {
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  library(mvtnorm) # for rmvnorm if allow.noise in fitSLR()
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
                  intercept = intercept)
  lambdamin.idx = which.min(cv.fits$cvm)
  betahat0 = as.numeric(cv.fits$int[lambdamin.idx])
  betahat = as.matrix(cv.fits$bet[, lambdamin.idx])
  # get fitted model
  predictSLR = function(X){
    betahat0 + computeBalances(X, cv.fits$btree) %*% betahat
  }
  # calculate prediction error
  Ypred = predictSLR(Xtest)
  (Ytest - Ypred)^2
}
dim(pred.err)
mse = colMeans(pred.err)
mse.mean = mean(mse)
mse.sd = sd(mse)
mse.se = mse.sd / sqrt(numReps)
data.frame(mean = mse.mean, sd = mse.sd, se = mse.se, 
           lower = mse.mean - 2 * mse.se, upper = mse.mean + 2 * mse.se)
### K = 5
# intercept = TRUE
# mean       sd       se    lower    upper
# 1 39.80968 32.34372 3.234372 33.34093 46.27842
# intercept = FALSE
# 1 31.41447 11.82889 1.182889 29.04869 33.78024
saveRDS(
  pred.err, 
  file = paste0(output_dir,
                "/slr_prediction", 
                "_int", intercept,
                "_K", cv.K, 
                "_seed", rng.seed,
                ".rds"))
