generateCovariates = function(mu, Sigma, U, n){
  logW <- mvrnorm(n = n, mu = mu, Sigma = Sigma) 
  W <- exp(logW)
  X <- sweep(W, 1, rowSums(W), FUN='/')
  ilrX = computeBalances(X, U = U)
  return(list("X" = X, "ilrX" = ilrX))
}

simulateBalanceReg = function(mu, Sigma, U, n, theta, sigma.noise){
  # generate X
  covariates = generateCovariates(mu, Sigma, U, n)
  
  # generate y
  y = covariates$ilrX %*% theta + rnorm(n) * sigma.noise
  
  return(list("X" = covariates$X, "ilrX" = covariates$ilrX, "y" = y))
}

