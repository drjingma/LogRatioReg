simulateBalanceReg = function(mu, Sigma, U, n){
  # make X
  logW <- mvrnorm(n = n, mu = mu, Sigma = Sigma) 
  W <- exp(logW)
  X <- sweep(W, 1, rowSums(W), FUN='/')
  ilrX = computeBalances(X, U = U)
  
  # generate y
  y = ilrX %*% theta + rnorm(n) * sigma_eps
  
  return(list("X" = X, "ilrX" = ilrX, "y" = y))
}
