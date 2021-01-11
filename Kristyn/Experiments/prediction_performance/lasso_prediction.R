getwd()

# libraries
library(limSolve) # for constrained lm
library(glmnet)

# set up parallelization
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# settings
tol = 1e-4

# Cross-validation
cv.n_lambda = 100
cv.K = 10

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
  cv.fits = cv.glmnet(x = log.X.prop, y = y, nlambda = cv.n_lambda, nfolds = cv.K)
  cv.exact = glmnet(x = log.X.prop, y = y, lambda = cv.fits$lambda.min)
  betahat = as.matrix(cv.exact$beta)
  names(betahat) = colnames(betahat)
  # get fitted model (no refitting)
  intercept = cv.exact$a0
  predictCLM = function(X){
    intercept + X %*% betahat
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
                      "/lasso_prediction", 
                      "_seed", rng.seed,
                      ".rds"))
