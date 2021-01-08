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

fitGlmnet = function(X, y, lambda, nlambda = 100, nfolds = 10){
  # apply to user-defined lambda sequence, if given. if not, let glmnet provide.
  cv_exact = cv.glmnet(x = X, y = y, lambda = lambda, nfolds = nfolds)
  # get a new lambda sequence (why?)
  lambda = log(cv_exact$lambda)
  lambda_new = exp(seq(max(lambda), min(lambda) + 2, length.out = nlambda))
  # fit cv.glmnet to choose lambda based on cross-validation
  cv_exact2 = cv.glmnet(x = X, y = y, lambda = lambda_new, nfolds = nfolds)
  # fit glmnet to the chosen lambda
  refit_exact = glmnet(x = X,y = y, family = 'gaussian', 
                       lambda = cv_exact2$lambda.min)
  return(list(cv.glmnet = cv_exact2, 
              glmnet = refit_exact, 
              beta = refit_exact$beta, 
              lambda = refit_exact$lambda))
}

# fit SLR given data (X, y), linkage, and (optional) lambda
# this version is the clean version
fitSLRLasso = function(X, y, linkage, lambda = NULL, nlambda, nfolds){
  btree = getSupervisedTree(X, y, linkage)
  Xb = computeBalances(btree, X)
  
  # lasso fit
  glm.temp = fitGlmnet(Xb, y, lambda, nfolds)
  
  return(list(
    glmnet = glm.temp$glmnet, 
    btree = btree,
    betahat = as.vector(glm.temp$beta), 
    lambda = glm.temp$lamdba
  ))
}
