################################################################################
# Not the same simulated data for classic & compositional lasso tests
# from LogRatioReg > RCode > sup-balances.R
################################################################################

# Supervised Log Ratios Regression

# do hierarchical clustering with specified linkage for data X and y
getSupervisedTree = function(X, y, linkage, allow.noise = FALSE){
  # browser()
  n = dim(X)[1]
  p = dim(X)[2]
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(NA, p, p)
  y_demeaned = y - mean(y)
  for(j in 1:p){
    for(k in 1:p){
      Zjk = log(X[, j]) - log(X[, k])
      # diagonal = 1
      if(j == k){
        cormat[j, k] = 1
      } else{
        if(all(Zjk == 0)){
          if(allow.noise){
            Zjk = rnorm(n, 0, 1e-10) # add noise if Zjk = 0vec
          } else{
            print(j)
            print(k)
            print(Zjk)
            stop(paste0("Zjk == 0vec for j = ", j, " and k = ", k ," leads to correlation NA"))
          }
        }
        # off-diagonal = cor(Zjk, y)
        Zjk_demeaned = Zjk - mean(Zjk)
        cormat[j, k] = cor(Zjk_demeaned, y_demeaned)
      }
    }
  }
  # give the rows and columns the names of taxa in X, 
  #   so that sbp.fromHclust works later
  rownames(cormat) = colnames(X)
  rownames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  Gammamat = 1 - abs(cormat)
  
  # get tree from hierarchical clustering
  btree_slr = hclust(as.dist(Gammamat), method = linkage)
  
  return(btree_slr)
}

# using a hierarchical tree, compute the balances for X
computeBalances = function(btree_slr, X){
  # compute balances from hclust object using balance pkg:
  # 1. build SBP (serial binary partition) matrix from hclust object
  sbp_slr = sbp.fromHclust(btree_slr)
  # 2. calculate balances from SBP matrix
  balances_slr = balance.fromSBP(X, sbp_slr)
  return(balances_slr)
}

fitGlmnet = function(X, y, lambda){
  cv_exact = cv.glmnet(x = X, y = y, lambda = lambda)
  lambda = log(cv_exact$lambda)
  lambda_new = exp(seq(max(lambda), min(lambda) + 2, length.out = 100))
  cv_exact2 = cv.glmnet(x = X, y = y, lambda = lambda_new)
  refit_exact = glmnet(x = X,y = y, family = 'gaussian', lambda = cv_exact2$lambda.min)
  return(list(beta = refit_exact$beta, lambda = refit_exact$lambda))
}

# fit SLR given data (X, y), linkage, and (optional) lambda
# this version is the clean version
fitSLRLasso = function(X, y, linkage, lambda = NULL, allow.noise = FALSE){
  btree = getSupervisedTree(X, y, linkage, allow.noise)
  Xb = computeBalances(btree, X)
  
  # lasso fit
  cv_exact = cv.glmnet(x = Xb, y = y, lambda = lambda)
  lambda = log(cv_exact$lambda)
  lambda_new = exp(seq(max(lambda), min(lambda) + 2, length.out = 100))
  cv_exact2 = cv.glmnet(x = Xb, y = y, lambda = lambda_new)
  refit_exact = glmnet(x = Xb,y = y, family = 'gaussian', lambda = cv_exact2$lambda.min)
  
  return(list(
    btree = btree,
    betahat = as.vector(refit_exact$beta), 
    lambda = refit_exact$lamdba
  ))
}
# tried using the fitting function in Dr. Ma's code, but changed my mind
fitSLRLasso0 = function(X, y, linkage, lambda = NULL, allow.noise = TRUE){
  # divide dataset into training and validation set
  #   (dataset might be a training set being further divided for validation)
  # n = dim(X)[1]
  # shuffle = sample(1:n)
  # # id_tv = (shuffle %% 2) + 1
  # id_tv = (shuffle %% 3) + 1
  # n_tv = as.vector(table(id_tv))
  # Xtrain.tmp = X[id_tv %in% c(1, 2), ]
  # Xvalid.tmp = X[id_tv == 3, ]
  # ytrain.tmp = y[id_tv %in% c(1, 2)]
  # yvalid.tmp = y[id_tv == 3]
  # # fit
  btree = getSupervisedTree(X, y, linkage, allow.noise)
  # Xtrainb = computeBalances(btree, Xtrain.tmp)
  # Xvalidb = computeBalances(btree, Xvalid.tmp)
  Xb = computeBalances(btree, X)
  
  # lasso fit
  glm.temp = fitGlmnet(Xb, y, lambda)
  
  return(list(
    btree = btree,
    betahat = as.vector(glm.temp$beta), 
    lambda = glm.temp$lamdba
    ))
}
