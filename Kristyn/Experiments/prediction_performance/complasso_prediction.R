# ``To compare the prediction performance of the two methods, we randomly divided the data into a
# training set of 70 subjects and a test set of 28 subjects, and used the fitted model chosen by cross-
#   validation based on the training set to evaluate the prediction error on the test set. The prediction
# error averaged over 100 replicates was 30·30 for the proposed method and 30·55 for lasso (ii),
# with standard errors of 0·97 and 1·04, respectively, suggesting that the prediction performance
# of the proposed method is similar to or better than that of lasso (ii).''
#  -- Lin et al 2014

getwd()

# libraries
library(limSolve) # for constrained lm

# set up parallelization
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Dr. Ma sources
source("RCode/func_libs.R")

# settings
tol = 1e-4

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
  library(limSolve)
  
  # split into train and test sets
  train.idx = sample(1:n, n.train)
  Xtrain = log.X.prop[train.idx, ]
  Ytrain = y[train.idx]
  Xtest = log.X.prop[-train.idx, ]
  Ytest = y[-train.idx]
  
  # refitted CV
  # Split the data into 10 folds
  
  # Fit Lasso on training set
  cv.fits = cv.func(
    method="ConstrLasso", y = Ytrain, x = Xtrain, 
    Cmat = matrix(1, dim(Xtrain)[2], 1), nlam = cv.n_lambda, 
    nfolds = cv.K, tol = tol)
  lambdamin.idx = which.min(cv.fits$cvm)
  betabar = cv.fits$bet[, lambdamin.idx]
  names(betabar) = colnames(betabar)
  # get fitted model (no refitting)
  intercept = cv.fits$int[lambdamin.idx]
  predictCLM = function(X){
    intercept + X %*% betabar
  }
  print(paste("sum of betas = ", sum(betabar)))
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
                      "/complasso_prediction", 
                      "_K", cv.K, 
                      "_seed", rng.seed,
                      ".rds"))
