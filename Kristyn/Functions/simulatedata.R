generateCovariates = function(mu, Sigma, sbp, n){
  logW <- mvrnorm(n = n, mu = mu, Sigma = Sigma) 
  W <- exp(logW)
  X <- sweep(W, 1, rowSums(W), FUN='/')
  ilrX = getIlrX(X = X, sbp = sbp)
  return(list("X" = X, "ilrX" = ilrX))
}

simulateBalanceReg = function(
  mu, Sigma, sbp, n, theta, sigma.noise
){
  # generate X
  covariates = generateCovariates(
    mu = mu, Sigma = Sigma, sbp = sbp, n = n)
  
  # generate y
  y = covariates$ilrX %*% theta + rnorm(n) * sigma.noise
  
  return(list("X" = covariates$X, "ilrX" = covariates$ilrX, "y" = y))
}

