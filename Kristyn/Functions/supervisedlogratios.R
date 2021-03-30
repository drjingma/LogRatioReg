# Supervised Log Ratios Regression
# do hierarchical clustering with specified linkage for compositional data X and y
getSupervisedTree = function(y, X, linkage = "complete", rho.type = "square"){
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
      if(rho.type == "square" | rho.type == "squared" | rho.type == "s" | 
         rho.type == 1){
        val = (cor(Zjk_demeaned, y_demeaned))^2
      } else{
        val = abs(cor(Zjk_demeaned, y_demeaned))
      }
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

# # old version, with noise arguments
# getSupervisedTree.old = function(
#   y, X, linkage = "complete",
#   allow.noise = FALSE, noise = 1e-12
# ){
#   # browser()
#   n = dim(X)[1]
#   p = dim(X)[2]
#   if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
#   
#   # calculate correlation of each pair of log-ratios with response y
#   cormat = matrix(1, p, p) # diagonal == 1
#   y_demeaned = y - mean(y)
#   for (j in 1:(p - 1)){
#     for (k in (j + 1):p){
#       Zjk = log(X[, j]) - log(X[, k])
#       Zjk_demeaned = Zjk - mean(Zjk)
#       if(all(Zjk == 0)){ # add noise (hopefully never necessary! hopefully never in this situation)
#         if(allow.noise) Zjk = Zjk + rmvnorm(n, rep(0, n), noise * diag(n))
#         # do not do anything, except maybe warn
#       } else{
#         val = abs(cor(Zjk_demeaned, y_demeaned))
#       }
#       cormat[j, k] = val
#       cormat[k, j] = val
#     }
#   }
#   # find out which columns give na
#   # for(i in 1:p) if(!all(!is.na(cormat[, i]))) print(paste0(i, ":", which(is.na(cormat[, i]))))
#   #
#   # give the rows and columns the names of taxa in X,
#   #   so that sbp.fromHclust works later
#   rownames(cormat) = colnames(X)
#   rownames(cormat) = colnames(X)
#   
#   # get dissimilarity matrix
#   Gammamat = 1 - cormat
#   
#   # get tree from hierarchical clustering
#   btree_slr = hclust(as.dist(Gammamat), method = linkage)
#   
#   return(btree_slr)
# }

# using a hierarchical tree, compute the balances for X
computeBalances = function(X, btree = NULL, sbp = NULL, U = NULL){
  # checks
  if(class(X) != "matrix"){
    if(class(X) == "numeric"){
      X = matrix(X, ncol = length(X))
    } else{
      warning("computeBalances: X is neither of class matrix nor numeric.")
    }
  }
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:ncol(X), sep = "")
  
  # get U
  if(is.null(U)){
    U = getU(btree, sbp)
  }
  
  # calculate balances from U
  balances = log(X) %*% U
  return(balances)
}

getU = function(btree = NULL, sbp = NULL){
  # get SBP matrix of contrasts (if not given)
  if(is.null(btree) & is.null(sbp)) stop("getU: need one of btree or sbp arg.")
  if(is.null(sbp) & !is.null(btree)) sbp = sbp.fromHclust(btree)
  
  # normalize SBP matrix to get U
  r.vec = colSums(sbp == 1)
  s.vec = colSums(sbp == -1)
  U = sbp
  for(i in 1:ncol(U)){
    U[(U[, i] == 1), i] = 1 / r.vec[i]
    # divide -1s by s
    U[(U[, i] == -1), i] = -1 / s.vec[i]
  }
  # multiply by sqrt(rs / (r + s))
  norm.const = sqrt((r.vec * s.vec) / (r.vec + s.vec))
  U = sweep(U, MARGIN = 2, STATS = norm.const, FUN = "*")
  # U2 = t(t(U) * norm.const) # equal to U
  return(U)
}


fitSLR = function(
  y, X, linkage = "complete", lambda = NULL, nlam = 20, intercept = TRUE,
  rho.type = "squared"
){
  # checks
  if(length(y) != nrow(X)) stop("fitSLR : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:ncol(X), sep = "")
  # check if lambda is given, assign get_lambda accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
  # compute balances
  btree = getSupervisedTree(y, X, linkage, rho.type)
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
  foldid = NULL, intercept = TRUE, rho.type = "squared"
){
  
  # compute balances
  btree = getSupervisedTree(y, X, linkage, rho.type)
  Xb = computeBalances(X, btree)
  
  # run cv.glmnet
  cv_exact = cv.glmnet(
    x = Xb, y = y, lambda = lambda, nlambda = nlam, nfolds = nfolds, 
    foldid = foldid, intercept = intercept, type.measure = c("mse"))
  
  # check lambda length
  if(nlam != length(cv_exact$lambda)){
    lambda <- log(cv_exact$lambda)
    lambda_new <- exp(seq(max(lambda), min(lambda),length.out = nlam))
    cv_exact = cv.glmnet(
      x = Xb, y = y, lambda = lambda_new, nfolds = nfolds, 
      foldid = foldid, intercept = intercept, type.measure = c("mse"))
  }
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

fitILR = function(
  y = NULL, X, btree = NULL, sbp = NULL, U = NULL, lambda = NULL, nlam = 20, 
  intercept = TRUE
){
  # checks
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:ncol(X), sep = "")
  # check if lambda is given, assign get_lambda accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
  # compute balances
  Xb = computeBalances(X, btree, sbp, U)
  
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

cvILR = function(
  y = NULL, X, btree = NULL, sbp = NULL, U = NULL, lambda = NULL, nlam = 20, 
  nfolds = 10, foldid = NULL, intercept = TRUE
){
  
  # compute balances
  Xb = computeBalances(X, btree, sbp, U)
  
  # run cv.glmnet
  cv_exact = cv.glmnet(
    x = Xb, y = y, lambda = lambda, nlambda = nlam, nfolds = nfolds, 
    foldid = foldid, intercept = intercept, type.measure = c("mse"))
  
  # check lambda length
  if(nlam != length(cv_exact$lambda)){
    lambda <- log(cv_exact$lambda)
    lambda_new <- exp(seq(max(lambda), min(lambda),length.out = nlam))
    cv_exact = cv.glmnet(
      x = Xb, y = y, lambda = lambda_new, nfolds = nfolds, 
      foldid = foldid, intercept = intercept, type.measure = c("mse"))
  }
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



getBeta = function(
  theta, btree = NULL, sbp = NULL, U = NULL
){
  # get U
  if(is.null(U)){
    U = getU(btree, sbp)
  }
  
  # calculate beta
  beta = U %*% theta
  return(beta)
}
