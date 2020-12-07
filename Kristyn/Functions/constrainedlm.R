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
  
  # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

backStandardize <- function(stdXY, betahat, scale = FALSE){
  if(scale){
    betahat.tilde = betahat
    betahat = diag(stdXY$weights) %*% betahat.tilde
  }
  betahat0 = as.numeric(stdXY$Ymean - crossprod(stdXY$Xmeans, betahat))
  return(list(
    betahat0 = betahat0, 
    betahat = betahat))
}

clm <- function(X, Y, Q){
  betahat = solve(crossprod(X), crossprod(X, Y))
  # compute beta.bar, constrained fit
  Q = as.matrix(rep(1, dim(X)[2]))
  betabar = betahat - solve(crossprod(X), Q) %*% 
    solve(crossprod(Q, solve(crossprod(X), Q)), crossprod(Q, betahat))
  return(betabar)
}
