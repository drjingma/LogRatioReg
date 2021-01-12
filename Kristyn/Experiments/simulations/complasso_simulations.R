# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit compositional Lasso to the data
# Date: 1/2021
# Notes: 

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

# Method Settings
tol = 1e-4

# Cross-validation settings
cv.n_lambda = 100
cv.K = 10

# Simulation settings
numSims = 500
n = 100
p = 200
rho = 0.2
beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
sigma_epsilon = 0.5
seed = 1
muW = c(
  rep(log(p), 5), 
  rep(0, p - 5)
)
SigmaW = matrix(0, p, p)
for(j in 1:p){
  for(k in j:p){
    SigmaW[j, k] = rho^(k - j)
  }
}
SigmaW = SigmaW + t(SigmaW) - diag(diag(SigmaW))

################################################################################
# Simulations #
################################################################################

bs.selected_variables = foreach(
  b = 1:numSims, 
  .combine = cbind, 
  .noexport = c("ConstrLassoC0")
) %dopar% {
  source("RCode/func_libs.R")
  library(limSolve)
  
  # simulate data
  # simulate covariates (compositional data)
  V = rmvnorm(n = 1, mean = muW, sigma = SigmaW) # after all that, n = numSims
  W = exp(V) # counts
  rowsumsW = apply(W, 1, sum)
  X = W / rowsumsW # compositions
  epsilon = rnorm(n, 0, sigma_epsilon) # noise
  Z = log(X)
  Y = Z %*% beta + epsilon # response
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # refitted CV
  # Split the data into 10 folds
  shuffle = sample(1:n)
  idfold = (shuffle %% cv.K) + 1
  n_fold = as.vector(table(idfold))
  
  # Do cross-validation
  # calculate squared error (prediction error?) for each fold, 
  #   needed for CV(lambda) calculation
  cvm = rep(NA, cv.n_lambda) # want to have CV(lambda)
  cvm_sqerror = matrix(NA, cv.K, cv.n_lambda)
  # Fit Lasso for each fold removed
  for (j in 1:cv.K){
    # Training data
    Xtrain = log.X.prop.bs[idfold != j, ]
    Ytrain = y.bs[idfold != j]
    # Test data
    Xtest = log.X.prop.bs[idfold == j, ]
    Ytest = y.bs[idfold == j]
    
    # Fit LASSO on that fold using fitLASSOcompositional
    Lasso_j = ConstrLasso(
      Ytrain, Xtrain, Cmat = matrix(1, dim(Xtrain)[2], 1), 
      nlam = cv.n_lambda, tol = tol)
    non0.betas = Lasso_j$bet != 0 # diff lambda = diff col
    for(m in 1:cv.n_lambda){
      selected_variables = non0.betas[, m]
      # get refitted coefficients, after model selection and w/o penalization
      if(all(!selected_variables)){ # if none selected
        refit = lm(Ytrain ~ 1)
        predictCLM = function(X) coefficients(refit)
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
      }
      # calculate squared error (prediction error?)
      Ypred = predictCLM(Xtest)
      cvm_sqerror[j, m] = sum(crossprod(Ytest - Ypred))
    }
  }
  
  # Calculate CV(lambda) for each value of lambda
  cvm = colMeans(cvm_sqerror)
  
  # Find lambda_min = argmin{CV(lambda)}
  lambda_min_index = which.min(cvm)
  lambda_min = Lasso_j$lambda[lambda_min_index]
  
  # final fit
  Lasso_select.bs = ConstrLasso(
    y.bs, log.X.prop.bs, Cmat = matrix(1, dim(log.X.prop)[2], 1), 
    lambda = lambda_min, nlam = 1, tol=tol)
  XYdata = data.frame(log.X.prop.bs, y = y.bs)
  non0.betas = Lasso_select.bs$bet != 0 # diff lambda = diff col
  selected_variables = non0.betas
  # record which variables were selected
  selected_variables
  
}
rownames(bs.selected_variables) = colnames(X)

bs.selected_variables_numeric = apply(bs.selected_variables, 2, as.numeric)
bs.selection_percentages = apply(bs.selected_variables_numeric, 1, FUN = 
                                   function(x) sum(x, na.rm = TRUE) / bs.n)
names(bs.selection_percentages) = rownames(bs.selected_variables)
bs.results = list(
  seed = rng.seed,  
  selected_variables = bs.selected_variables, 
  selection_percentages = bs.selection_percentages
)

saveRDS(bs.results,
        file = paste0("Kristyn/Experiments/output",
                      "/complasso_selection", 
                      "_refit",
                      "_B", bs.n, 
                      "_seed", rng.seed,
                      ".rds"))

sort(bs.selection_percentages)


