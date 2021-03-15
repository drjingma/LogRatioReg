################################################################################
# Not the same simulated data for classic & compositional lasso tests
# from LogRatioReg > RCode > sup-balances.R
################################################################################

# Supervised Log Ratios Regression
# do hierarchical clustering with specified linkage for compositional data X and y
getSupervisedTree = function(
  y, X, linkage = "complete", allow.noise = FALSE, noise = 1e-12){
  n = dim(X)[1]
  p = dim(X)[2]
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(1, p, p) # diagonal == 1
  y_demeaned = y - mean(y)
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      Zjk_demeaned = Zjk - mean(Zjk)
      # if(all(Zjk == 0)){
      #   # do not do anything, except maybe warn
      # } else{
      #   val = abs(cor(Zjk_demeaned, y_demeaned))
      # }
      val = abs(cor(Zjk_demeaned, y_demeaned))
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # give the rows and columns the names of taxa in X, for sbp.fromHclust()
  rownames(cormat) = colnames(X)
  rownames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  Gammamat = 1 - cormat
  
  # get tree from hierarchical clustering
  btree_slr = hclust(as.dist(Gammamat), method = linkage)
  
  return(btree_slr)
}

# do hierarchical clustering with specified linkage for compositional data X and y
getSupervisedTree.old = function(
  y, X, linkage = "complete", allow.noise = FALSE, noise = 1e-12){
  # browser()
  n = dim(X)[1]
  p = dim(X)[2]
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(1, p, p) # diagonal == 1
  y_demeaned = y - mean(y)
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      Zjk_demeaned = Zjk - mean(Zjk)
      if(all(Zjk == 0)){ # add noise (hopefully never necessary! hopefully never in this situation)
        if(allow.noise) Zjk = Zjk + rmvnorm(n, rep(0, n), noise * diag(n))
        # do not do anything, except maybe warn
      } else{
        val = abs(cor(Zjk_demeaned, y_demeaned))
      }
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # find out which columns give na
  # for(i in 1:p) if(!all(!is.na(cormat[, i]))) print(paste0(i, ":", which(is.na(cormat[, i]))))
  # 
  # give the rows and columns the names of taxa in X, 
  #   so that sbp.fromHclust works later
  rownames(cormat) = colnames(X)
  rownames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  Gammamat = 1 - cormat
  
  # get tree from hierarchical clustering
  btree_slr = hclust(as.dist(Gammamat), method = linkage)
  
  return(btree_slr)
}

# using a hierarchical tree, compute the balances for X
computeBalances = function(X, btree){
  # checks
  if(class(X) != "matrix"){
    if(class(X) == "numeric"){
      X = matrix(X, ncol = length(X))
    } else{
      warning("computeBalances: X is neither of class matrix nor numeric.")
    }
  }
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:ncol(X), sep = "")
  # compute balances from hclust object using balance pkg:
  # 1. build SBP (serial binary partition) matrix from hclust object
  U = sbp.fromHclust(btree) # U = basis vectors
  # 2. calculate balances from SBP matrix
  # balances = balance.fromSBP(X, U)
  balances = log(X) %*% U
  return(balances)
}

fitSLR = function(
  y, X, linkage = "complete", lambda = NULL, nlam = 20, intercept = TRUE,
  allow.noise = FALSE, noise = 1e-12
){
  # checks
  if(length(y) != nrow(X)) stop("fitSLR : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:ncol(X), sep = "")
  # check if lambda is given, assign get_lambda accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
  # compute balances
  btree = getSupervisedTree(y, X, linkage, allow.noise, noise)
  Xb = computeBalances(X, btree)
  
  # run glmnet
  glmnet.fit = glmnet(x = Xb, y = y, lambda = lambda, nlambda = nlam, 
                      intercept = intercept)
  # check lambda length
  if(nlam != length(glmnet.fit$lambda)){
    lambda <- log(glmnet.fit$lambda)
    lambda_new <- exp(seq(max(lambda), min(lambda),length.out = nlam))
    glmnet.fit = glmnet(x = Xb, y = y, lambda = lambda_new)
  }
  return(list(
    int = glmnet.fit$a0, 
    bet = glmnet.fit$beta, 
    lambda = glmnet.fit$lambda,
    glmnet = glmnet.fit, 
    btree = btree
  ))
}

cvSLR = function(
  y, X, linkage = "complete", lambda = NULL, nlam = 20, nfolds = 10, 
  foldid = NULL, intercept = TRUE, allow.noise = FALSE, noise = 1e-12
){
  
  # compute balances
  btree = getSupervisedTree(y, X, linkage, allow.noise, noise)
  Xb = computeBalances(X, btree)
  
  # call cv.glmnet
  cv_exact = cv.glmnet(
    x = Xb, y = y, lambda = lambda, nlambda = nlam, nfolds = nfolds, 
    foldid = foldid, intercept = intercept, type.measure = c("mse"))
  return(list(
    int = cv_exact$glmnet.fit$a0,
    bet = cv_exact$glmnet.fit$beta,
    lambda = cv_exact$lambda,
    cv.glmnet = cv_exact,
    glmnet = cv_exact$glmnet.fit,
    btree = btree,
    cvm = cv_exact$cvm
  ))
}

cvSLR.old = function(
  y, X, linkage = "complete", lambda = NULL, nlam = 20, nfolds = 10, 
  intercept = TRUE, allow.noise = FALSE, noise = 1e-12
){
  n <- nrow(X)
  p <- ncol(X)
  
  # Fit to the original data
  slr = fitSLR(y = y, X = X, linkage = linkage, lambda = lambda, nlam = nlam, 
               intercept = intercept, allow.noise = allow.noise, noise = noise)
  bet = slr$bet
  int = slr$int
  lambda = slr$lambda
  nlam = length(lambda)
  btree = slr$btree
  
  # Split the data into K folds
  shuffle = sample(1:n)
  idfold = (shuffle %% nfolds) + 1 # IDs for which fold each obs goes into
  n_fold = as.vector(table(idfold)) # number of things in each fold
  
  cvm_sqerror = matrix(rep(NA, nfolds * nlam), nfolds, nlam) # calculate squared error for each fold, needed for CV(lambda) calculation
  cvse_errorsd =  matrix(rep(NA, nfolds * nlam), nfolds, nlam) # calculate sd of error for each fold, needed for SE_CV(lambda) calculation
  # Calculate Lasso for each fold removed
  for (j in 1:nfolds){
    # Training data
    Xtrain = X[idfold != j, ]
    Ytrain = y[idfold != j]
    # Test data
    Xtest = X[idfold == j, ]
    Ytest = y[idfold == j]
    
    # Calculate LASSO on that fold using fitLASSO
    slr_j = fitSLR(y = Ytrain, X = Xtrain, linkage = linkage, lambda = lambda, 
                   nlam = nlam, intercept = intercept, 
                   allow.noise = allow.noise, noise = noise)
    # Any additional calculations that are needed for calculating CV and SE_CV(lambda)
    for(m in 1:nlam){
      Ypred = slr_j$int[m] + 
        computeBalances(Xtest, btree) %*% as.matrix(slr_j$bet[ , m])
      cvm_sqerror[j, m] = sum(crossprod(Ytest - Ypred)) / length(Ytest)
      cvse_errorsd[j, m] = cvm_sqerror[j, m] / n_fold[j]
    }
  }
  
  # Calculate CV(lambda) and SE_CV(lambda) for each value of lambda
  cvse = sqrt(apply(cvse_errorsd, 2, var)) / sqrt(nfolds)
  cvm = colMeans(cvm_sqerror)
  
  # Find lambda_min = argmin{CV(lambda)}
  lambda_min_idx = which.min(cvm)
  lambda_min = lambda[lambda_min_idx]
  
  # Find lambda_1se = maximal lambda s.t. CV(lambda) <= CV(lambda_min) + CV_SE(lambda_min)
  bound = cvm[lambda_min_idx] + cvse[lambda_min_idx]
  lambda_1se_idx = which(cvm <= bound)
  lambda_1se = lambda[min(lambda_1se_idx)]
  
  return(list(
    int = int, 
    bet = bet, 
    lambda = lambda, 
    btree = btree,
    lambda_min = lambda_min, 
    lambda_1se = lambda_1se, 
    cvm = cvm, 
    cvm_idx = lambda_min_idx,
    cvse = cvse,
    cvse_idx = lambda_1se_idx, 
    foldid = idfold
  ))
}

# 
LRtoLC = function(
  LRcoefficients, btree
){
  U = sbp.fromHclust(btree) # matrix of basis vectors, p x (p - 1)
  LCcoefficients = U %*% LRcoefficients
  return(LCcoefficients)
}
