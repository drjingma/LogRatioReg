# Supervised Log Ratios Regression



# do hierarchical clustering with specified linkage 
#   for compositional data X and y
getSupervisedTree = function(y, X, linkage = "complete", rho.type = "square"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(0, p, p) # diagonal == 1
  y_demeaned = y - mean(y)
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      Zjk_demeaned = Zjk - mean(Zjk)
      if(rho.type == "square" | rho.type == "squared" | rho.type == "s" | 
         rho.type == 2){
        # val = (cor(Zjk_demeaned, y_demeaned))^2
        val = (cor(Zjk, y))^2
      } else{
        # val = abs(cor(Zjk_demeaned, y_demeaned))
        val = abs(cor(Zjk, y))
      }
      # if(is.na(val)) stop("getSupervisedTree() : correlation = 0")
      # hopefully we never have to use this line below.
      # if(is.na(val)) val = 1 #######################################################################################
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



# Compute the balances for X
# using a hierarchical tree or SBP matrix or transformation matrix U
computeBalances = function(X, btree = NULL, sbp = NULL, U = NULL){
  p = dim(X)[2]
  
  # checks
  if(!("matrix" %in% class(X))){
    if("numeric" %in% class(X)){
      X = matrix(X, ncol = length(X))
    } else{
      warning("computeBalances: X is neither of class matrix nor numeric.")
    }
  }
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  
  # get U
  if(is.null(U)){
    U = getU(btree, sbp)
  }
  
  # calculate balances from U
  balances = log(X) %*% U
  return(balances)
}



# Compute the transformation matrix U
#   using a hierarchical tree or SBP matrix
getU = function(btree = NULL, sbp = NULL){
  
  # get SBP matrix of contrasts (if not given)
  if(is.null(btree) & is.null(sbp)) stop("getU: need one of btree or sbp arg.")
  if(is.null(sbp) & !is.null(btree)) sbp = sbp.fromHclust(btree)
  
  # normalize SBP matrix to get U
  r.vec = colSums(sbp == 1)
  s.vec = colSums(sbp == -1)
  U = sbp
  for(i in 1:ncol(U)){
    # divide 1s by r
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



# Fit supervised log-ratios model to compositional data X and response y
fitSLR = function(
  y, X, linkage = "complete", lambda = NULL, nlam = 20, intercept = TRUE,
  rho.type = "squared"
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("fitSLR : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  # check if lambda is given, assign nlam accordingly
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
    # if it's not nlam, refit to nlam lambdas
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
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("fitSLR : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  # check if lambda is given, assign nlam accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
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
    # note: one difference (that causes results to be slightly different)
    #   between my implementation here and the one in sup-balances.R is that
    #   in sup-balances.R the below line is instead
    #   lambda_new <- exp(seq(max(lambda), min(lambda)+2,length.out = nlam))
    #   -- it doesn't seem to be the biggest cause for difference, however.
    lambda_new <- exp(seq(max(lambda), min(lambda),length.out = nlam))
    cv_exact = cv.glmnet(
      x = Xb, y = y, lambda = lambda_new, nfolds = nfolds, 
      # note: yet another difference is that here, I specify type.measure to be 
      #   "mse", rather than the default, "deviance",
      #   but I think this should not be an issue since we are using a gaussian
      #   model so they end up the same.
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


# Fit any balance regression model to compositional data X and response y
#   using the balances' associated binary tree, SBP matrix, or 
#   transformation matrix U
fitILR = function(
  y = NULL, X, btree = NULL, sbp = NULL, U = NULL, lambda = NULL, nlam = 20, 
  intercept = TRUE
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("fitILR : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  # check if lambda is given, assign nlam accordingly
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
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("fitILR : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  # check if lambda is given, assign nlam accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
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



# Obtain beta (or betahat)
#   by transforming from balance regression model's theta (or thetahat)
#   to the linear log-contrasts model with parameters beta
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

getTPR = function(
  type = "ilr", beta, thetahat = NULL, sbp = NULL, betahat = NULL
){
  if(is.null(beta)) stop("getTPR: can't compute TPR without beta!")
  if(type == "ilr" | type == 1 | type == "balance"){
    if(is.null(sbp)) stop("getTPR: for ilr model, need sbp matrix!")
    if(is.null(thetahat)) stop("getTPR: for ilr model, need thethat!")
    # betahat = getBeta(thetahat, btree.slr)
    non0.beta = (beta != 0)
    non0.thetahat = (thetahat != 0)
    sel.cols.SBP = sbp[, non0.thetahat, drop = FALSE]
    non0.betahat = apply(sel.cols.SBP, 1, function(row) any(row != 0))
    # TPR & S.hat
    TPR = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat = sum(non0.betahat)
  } else if(type == "llc" | type == 2 | type == "lc" | type == "logcontrast" | 
            type == "linearlogcontrast"){
    if(is.null(betahat)) stop("getTPR: for linear log contrast model, need betahat!")
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat = sum(non0.betahat)
  } else{
    stop("getTPR: invalid type -- can only be ilr model or linear log contrast model")
  }
  return(c("S.hat" = S.hat, "TPR" = TPR))
}
