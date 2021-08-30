# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit supervised log-ratios method to the data
# Date: 08/29/2021
# Notes: cv.glmnet for predictions, but beta from refit.glmnet, like in aim31

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output_aim31_metrics_roc"

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

# Other simulation settings
numSims = 100

################################################################################
# Simulations #
################################################################################

registerDoRNG(rng.seed)
res = foreach(
  b = 1:numSims
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  
  library(balance) # for sbp.fromHclust()
  source("RCode/func_libs.R")
  source("Kristyn/Functions/supervisedlogratios.R")
  
  # helper functions
  run.glmnet <- function(x,y,xt,yt,lambda=NULL){
    cv_exact <- cv.glmnet(x=x,y=y,lambda=lambda)
    # lambda <- log(cv_exact$lambda)
    # lambda_new <- exp(seq(max(lambda),min(lambda)+2,length.out = 100))
    # cv_exact <- cv.glmnet(x=x,y=y,lambda=lambda_new)
    refit_exact <- glmnet(x=x,y=y,family='gaussian',lambda=cv_exact$lambda.min)
    pred_exact <- predict(cv_exact, newx = xt, type = "response", s = 'lambda.min')
    return(list(beta.min=refit_exact$beta,
                beta=cv_exact$glmnet.fit$beta,
                lambda=cv_exact$lambda,
                mse.pred = mean((pred_exact-yt)^2), 
                cv_exact = cv_exact, 
                refit_exact = refit_exact))
  }
  
  roc.for.coef <- function(beta_hat, beta, eps = 1e-08){
    TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
    FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
    tpr <- TP/(TP + FN)
    S_hat <- sum((abs(beta_hat) > eps))
    out <- c(S_hat,tpr)
    names(out) <- c('S_hat','tpr')
    return(out)
  }
  
  roc.for.coef.LR <- function(beta_hat,beta,sbp,eps=1e-08){
    
    if (is.null(sbp)){
      stop('A sequential binary partition tree is needed for roc evaluation!')
    }
    
    # first identify the variable at the LR scale
    index <- which(abs(beta_hat) > eps)
    
    # map to original variable
    if (length(index)==0){
      S_hat <- NULL
    } else  if (length(index)==1){
      S_hat <- names(which(abs(sbp[,index])>0))
    } else {
      S_hat <- names(which(rowSums(abs(sbp[,index]))>0))
    }
    S0 <- names(which((abs(beta) > eps)))
    TP <- intersect(S_hat, S0)
    tpr <- length(TP)/length(S0)
    out <- c(length(S_hat),tpr)
    names(out) <- c('S_hat','tpr')
    return(out)
  }
  
  # Settings to toggle with
  rho.type = "square" # 1 = "absolute value", 2 = "square"
  beta.settings = "new"
  linkage = "average"
  tol = 1e-4
  nlam = 100
  intercept = TRUE
  K = 10
  n = 100
  p = 200
  rho = 0.5 # 0.2, 0.5
  scaling = TRUE
  
  # which beta?
  if(beta.settings == "old" | beta.settings == "linetal2014"){
    beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
  } else{
    beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
  }
  
  # Population parameters
  non0.beta = (beta != 0)
  sigma_eps = 0.5
  seed = 1
  muW = c(rep(log(p), 5), rep(0, p - 5))
  SigmaW <- rgExpDecay(p,rho)$Sigma
  
  # for beta selection accuracy metrics
  slr.non0.beta = abs(beta) > 10e-8
  is0.beta = abs(beta) <= 10e-8
  
  # beta sparsity
  bspars = sum(non0.beta)
  
  file.end = paste0(
    "_dim", n, "x", p, 
    "_", beta.settings, 
    "_rho", rho, 
    "_int", intercept,
    "_scale", scaling,
    "_K", K,
    "_seed", rng.seed,
    ".rds")
  
  ##############################################################################
  
  # simulate data
  logW <- mvrnorm(n = n * 2, mu = muW, Sigma = SigmaW) 
  W <- exp(logW)
  colnames(W) <- paste0('s', 1:p)
  XAll <- sweep(W, 1, rowSums(W), FUN='/')
  ZAll = log(XAll)
  names(beta) <- colnames(W)
  yAll <-  ZAll %*% beta + rnorm(n) * sigma_eps
  
  # subset out training and test sets
  X = XAll[1:n, ]
  X.test = XAll[-(1:n), ]
  Z = ZAll[1:n, ]
  Z.test = ZAll[-(1:n)] 
  Y <- yAll[1:n,]
  Y.test <- yAll[-(1:n),]
  
  ##############################################################################
  
  # tuning parameter
  nlam <- 100
  lambda <- exp(seq(from=2, to=log(1e-4), length.out=100))
  
  # apply supervised log-ratios, using CV to select lambda
  S <- matrix(0,p,p)
  rownames(S) <- colnames(S) <- colnames(W)
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      newx <- log(X[1:n,j]) - log(X[1:n,k])
      newx <- newx - mean(newx)
      newy <- Y - mean(Y)
      S[j,k] <- S[k,j] <- abs(cor(newx,Y))^2
    }
  }
  btree <- hclust(as.dist(1-S),method = linkage)
  
  sbp_supervised = sbp.fromHclust(btree)
  ba_supervised <- balance.fromSBP(XAll, sbp_supervised)
  fit_supervised <- run.glmnet(
    x = ba_supervised[1:n,], y = Y, xt = ba_supervised[-(1:n),], yt = Y.test, 
    lambda = lambda)
  
  # evaluate model #
  
  # roc
  slr.roc <- apply(fit_supervised$beta, 2, function(a) 
    roc.for.coef.LR(a, beta, sbp_supervised)) # used cv.glmnet beta matrix
  
  saveRDS(slr.roc, file = paste0(output_dir, "/slr_roc_samelam", b, file.end))
  
  ##############################################################################
  
  # apply compositional lasso, using CV to select lambda
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = Z, Cmat = matrix(1, p, 1), nlam = nlam, 
    lambda = lambda, 
    nfolds = K, tol = tol, intercept = intercept, scaling = scaling)
  
  # choose lambda
  cl.lam.min.idx = which.min(complasso$cvm)
  cl.lam.min = complasso$lambda[cl.lam.min.idx]
  cl.a0 = complasso$int[cl.lam.min.idx]
  cl.betahat = complasso$bet[, cl.lam.min.idx]
  
  # evaluate model #
  
  # roc
  cl.roc <- apply(complasso$bet, 2, function(a) 
    roc.for.coef(a, beta)) # used cv.glmnet beta matrix
  
  saveRDS(cl.roc, file = paste0(output_dir, "/classo_roc_samelam", b, file.end))
  
}


