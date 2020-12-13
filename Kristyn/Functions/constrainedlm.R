# helper functions
standardize <- function(X, Y, center = FALSE, scale = FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # Center Y
  Ymean = mean(Y)
  if(center){
    Ytilde = Y - mean(Y)
  } else{
    Ytilde = Y
  }
  
  # Center and scale X
  Xmeans = colMeans(X)
  Xcen = X - matrix(Xmeans, n, p, byrow=T)
  normsX2 = colSums(Xcen^2) / n
  weights = 1 / sqrt(normsX2)
  if(center & scale){
    Xtilde = Xcen %*% diag(weights)
  } else if(center & !scale){
    Xtilde = Xcen
  } else if(!center & scale){
    Xtilde = X %*% diag(weights)
  } else{
    Xtilde = X
  }
  colnames(Xtilde) = colnames(X)
  
  # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

backStandardize <- function(stdXY, betahat, scale = FALSE){
  which.covariates = rownames(betahat)
  if(scale){
    betahat.tilde = betahat
    wts = stdXY$weights[which.covariates]
    betahat = diag(wts) %*% betahat.tilde
  }
  betahat0 = as.numeric(stdXY$Ymean - crossprod(stdXY$Xmeans[which.covariates], 
                                                betahat))
  return(list(
    betahat0 = betahat0, 
    betahat = betahat))
}

clm <- function(X, Y, Q){
  betahat = coefficients(lm(Y ~ -1 + X))
  names(betahat) = colnames(X)
  betahat = na.omit(betahat)
  if(!is.null(na.action(betahat))){
    warning("clm : some parameters are ommitted because NA in betahat")
    X = X[, names(betahat)]
    Q = as.matrix(rep(1, length(betahat)))
  }
  # compute beta.bar, constrained fit
  XtXinvQ = solve(crossprod(X), Q)
  betabar = betahat - XtXinvQ %*% 
    solve(crossprod(Q, XtXinvQ), crossprod(Q, betahat))
  return(betabar)
}
