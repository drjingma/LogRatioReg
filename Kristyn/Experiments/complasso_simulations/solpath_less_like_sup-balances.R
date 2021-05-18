# last updated: 05/11/2021
# simulate a data set using settings in Lin et al 2014, 
# fit Complasso and SLR, and plot the solution path.
# Unlike solpath_like_sup-balances, I'm fitting something more like my slr code.

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output"

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)
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

# ggplot
library(ggplot2)

# Method Settings
tol = 1e-4
nlam = 3 # for testing
intercept = TRUE
K = 10#5
rho.type = "abs" #"square"
linkage = "average"

# Simulation settings
numSims = 50#100
n = 50#100
p = 20#200
rho = 0.5 # 0.2, 0.5
# which beta?
beta.settings = "new"
if(beta.settings == "old" | beta.settings == "linetal2014"){
  beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
} else{
  beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
}
names(beta) = paste0('s',1:p)
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
# Sigma0 = rgExpDecay(p,rho)$Sigma
# all.equal(SigmaW, Sigma0)

################################################################################
# Simulations -- different lambda sequence #
################################################################################

registerDoRNG(rng.seed)
sims = foreach(
  b = 1:numSims
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
  # simulate training data #
  # generate W
  W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
  colnames(W) = paste0('s',1:p)
  # let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
  V = exp(W)
  rowsumsV = apply(V, 1, sum)
  X = V / rowsumsV
  epsilon = rnorm(n, 0, sigma_eps)
  Z = log(X)
  # generate Y
  Y = Z %*% beta + epsilon
  
  # apply compositional lasso
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), 
    nlam = nlam, nfolds = K, tol = tol, intercept = intercept)
  
  # apply supervised log-ratios
  # slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
  #             rho.type = rho.type)
  # S <- matrix(0,p,p)
  # rownames(S) <- colnames(S) <- colnames(X)
  # for (j in 1:(p-1)){
  #   for (k in (j+1):p){
  #     newx <- log(X[1:n,j]) - log(X[1:n,k])
  #     newx <- newx - mean(newx)
  #     newy <- Y - mean(Y)
  #     S[j,k] <- S[k,j] <- abs(cor(newx,Y))
  #     # S[j,k] <- S[k,j] <- abs(crossprod(newx,newy)/(sqrt(crossprod(newx)) * sqrt(crossprod(newy))))
  #   }
  # }
  # h_supervised <- hclust(as.dist(1-S),method = linkage)
  h_supervised <- getSupervisedTree(Y, X, linkage = linkage, rho.type = rho.type)
  # plot(h_supervised)
  # ggtree(h_supervised) + geom_point(aes(shape=isTip, color=isTip), size=3)
  sbp_supervised <- sbp.fromHclust(h_supervised)
  # ba_supervised <- balance.fromSBP(X,sbp_supervised)
  Xb = computeBalances(X, sbp = sbp_supervised)
  cv_exact = cv.glmnet(x = Xb, y = Y, nlambda = nlam, nfolds = K)
  # cv_exact <- cv.glmnet(x=ba_supervised[1:n,],y=Y,lambda=NULL)
  lambda <- log(cv_exact$lambda)
  lambda_new <- exp(seq(max(lambda),min(lambda),length.out = nlam))
  cv_exact <- cv.glmnet(x = Xb, y = Y, lambda = lambda_new, 
                        type.measure = c("mse"))
  
  list(X = X, Y = Y,
       fit.cl = complasso, fit.slr = cv_exact
  )
}



# organize the simulation results to have more easily-accessible matrices #
# cl
lambda.cl.mat = matrix(NA, nlam, numSims)
# slr
lambda.slr.mat = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  sim.tmp = sims[[i]]
  # cl
  lambda.cl.mat[, i] = sim.tmp$fit.cl$lambda
  # slr
  lambda.slr.mat[, i] = sim.tmp$fit.slr$lambda
}



################################################################################
# using the same lambda sequence in each method
nlam = 100#200
lambda.df = data.frame(
  complasso = exp(seq(max(log(lambda.cl.mat)), 
                      min(log(lambda.cl.mat)), length.out = nlam)),
  slr = exp(seq(max(log(lambda.slr.mat)), 
                min(log(lambda.slr.mat)),length.out = nlam))
)


registerDoRNG(rng.seed)
sims3 = foreach(
  b = 1:numSims
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
  print(paste0("starting sim ", b))
  
  # get simulated data #
  X = sims[[b]]$X
  Y = sims[[b]]$Y
  
  # apply compositional lasso
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), 
    lambda = lambda.df$complasso, nfolds = K, tol = tol, intercept = intercept)
  
  cl.res = apply(complasso$bet, 2, function(a) tpr.for.coef(beta, a))
  S.hat.cl = cl.res[1, ]
  TPR.cl = cl.res[2, ]
  
  # apply supervised log-ratios
  # slr = cvSLR(
  #   y = Y, X = X, lambda = lambda.df$slr, nfolds = K, intercept = intercept, 
  #   rho.type = rho.type)
  # S <- matrix(0,p,p)
  # rownames(S) <- colnames(S) <- colnames(X)
  # for (j in 1:(p-1)){
  #   for (k in (j+1):p){
  #     newx <- log(X[1:n,j]) - log(X[1:n,k])
  #     newx <- newx - mean(newx)
  #     newy <- Y - mean(Y)
  #     S[j,k] <- S[k,j] <- abs(cor(newx,Y))
  #     # S[j,k] <- S[k,j] <- abs(crossprod(newx,newy)/(sqrt(crossprod(newx)) * sqrt(crossprod(newy))))
  #   }
  # }
  # h_supervised <- hclust(as.dist(1-S),method = "average") # linkage = average
  h_supervised <- getSupervisedTree(Y, X, linkage = linkage, rho.type = rho.type)
  # plot(h_supervised)
  # ggtree(h_supervised) + geom_point(aes(shape=isTip, color=isTip), size=3)
  
  sbp_supervised <- sbp.fromHclust(h_supervised)
  # ba_supervised <- balance.fromSBP(X, sbp_supervised)
  Xb <- computeBalances(X, sbp = sbp_supervised)
  slr <- cv.glmnet(x = Xb, y = Y, lambda = lambda.df$slr, nfolds = K, 
                   type.measure = c("mse"))
  
  slr.res = apply(slr$glmnet.fit$beta, 2, 
                  function(a) tpr.for.coef.ilr(beta, a, sbp_supervised))
  S.hat.slr = slr.res[1, ]
  TPR.slr = slr.res[2, ]
  
  list(X = X, Y = Y, 
       fit.cl = complasso, TPR.cl = TPR.cl, S.hat.cl = S.hat.cl, 
       fit.slr = slr, TPR.slr = TPR.slr, S.hat.slr = S.hat.slr)
}

