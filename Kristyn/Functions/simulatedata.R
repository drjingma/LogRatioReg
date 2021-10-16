simulateBalanceReg = function(mu, Sigma, U, n, theta, sigma.noise){
  # make X
  logW <- mvrnorm(n = n, mu = mu, Sigma = Sigma) 
  W <- exp(logW)
  X <- sweep(W, 1, rowSums(W), FUN='/')
  ilrX = computeBalances(X, U = U)
  
  # generate y
  y = ilrX %*% theta + rnorm(n) * sigma.noise
  
  return(list("X" = X, "ilrX" = ilrX, "y" = y))
}
