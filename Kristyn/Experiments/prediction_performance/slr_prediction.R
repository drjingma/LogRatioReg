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
refit = FALSE
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
  Xtrain = X.prop[train.idx, ]
  Ytrain = y[train.idx]
  Xtest = X.prop[-train.idx, ]
  Ytest = y[-train.idx]
  
  # refitted CV
  # Split the data into 10 folds
  
  # Fit Lasso on training set
  cv.fits = cvSLR(y = Ytrain, X = Xtrain, nlam = cv.n_lambda, nfolds = cv.K)
  lambdamin.idx = which.min(cv.fits$cvm)
  betabar = cv.fits$bet[, lambdamin.idx]
  names(betabar) = colnames(betabar)
  # get fitted model (no refitting)
  if(refit){
    selected_variables = (betabar != 0)
    # get refitted coefficients, after model selection and w/o penalization
    if(all(!selected_variables)){ # if none selected
      lmfit = lm(Ytrain ~ 1)
      predictCLM = function(X) coefficients(lmfit)
    } else{ # otherwise, fit on selected variables
      # fit to the subsetted data
      Xtrain.sub = Xtrain[, selected_variables, drop = FALSE]
      Xtrain.sub2 = cbind(1, Xtrain.sub)
      Q = as.matrix(rep(1, sum(selected_variables)))
      Q2 = rbind(0, Q)
      colnames(Xtrain.sub2)[1] = "Intercept"
      rownames(Q2) = colnames(Xtrain.sub2)
      lsei.fit = lsei(A = Xtrain.sub2, B = Ytrain, E = t(Q2), F = 0)
      predictCLM = function(X){
        lsei.fit$X[1] +
          X[, selected_variables, drop = FALSE] %*% lsei.fit$X[-1]
      }
      print(paste("sum of betas = ", sum(lsei.fit$X[-1])))
    }
  } else{
    intercept = cv.fits$int[lambdamin.idx]
    predictCLM = function(X){
      intercept + X %*% betabar
    }
  }
  print(paste("sum of betas = ", sum(betabar)))
  # calculate squared error (prediction error?)
  Ypred = predictCLM(Xtest)
  (Ytest - Ypred)^2
}
dim(pred.err)

mean.pred.err = colMeans(pred.err)
print(paste0("mean prediction error: ", mean(mean.pred.err)))
print(paste0("standard deviation: ", (sd(mean.pred.err))))
print(paste0("standard error: ", (sd(mean.pred.err)) / rep.n))

saveRDS(pred.err,
        file = paste0("Kristyn/Experiments/output",
                      "/slr_prediction", 
                      "_refit", refit,
                      "_seed", rng.seed,
                      ".rds"))
