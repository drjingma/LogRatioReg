generateCovariates = function(mu, Sigma, btree, sbp, U, n){
  logW <- mvrnorm(n = n, mu = mu, Sigma = Sigma) 
  W <- exp(logW)
  X <- sweep(W, 1, rowSums(W), FUN='/')
  ilrX = computeBalances(X = X, btree = btree, sbp = sbp, U = U)
  return(list("X" = X, "ilrX" = ilrX))
}

simulateBalanceReg = function(
  mu, Sigma, btree = NULL, sbp = NULL, U = NULL, n, theta, sigma.noise
){
  # generate X
  covariates = generateCovariates(
    mu = mu, Sigma = Sigma, btree = btree, sbp = sbp, U = U, n = n)
  
  # generate y
  y = covariates$ilrX %*% theta + rnorm(n) * sigma.noise
  
  return(list("X" = covariates$X, "ilrX" = covariates$ilrX, "y" = y))
}

