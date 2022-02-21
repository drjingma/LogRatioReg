get_sigma_eps = function(sbp1, ilr.trans.constant, theta, Rsq, rho){
  # calculate Var[ilr(X)]
  pos.idx = which(sbp1 == 1)
  num.pos = length(pos.idx)
  neg.idx = which(sbp1 == -1)
  num.neg = length(neg.idx)
  pairs.pos <- pairs.neg <- c()
  if(num.pos > 1){
    pairs.pos = t(combn(pos.idx, 2))
  }
  if(num.neg > 1){
    pairs.neg = t(combn(neg.idx, 2))
  }
  Sigma.ij.pos = rho^abs(pairs.pos[, 1] - pairs.pos[, 2])
  Sigma.ij.neg = rho^abs(pairs.neg[, 1] - pairs.neg[, 2])
  varIlrX = ilr.trans.constant^2 * (
    1 / num.pos + 1 / num.neg + 
      2 * sum(Sigma.ij.pos) / num.pos^2 + 2 * sum(Sigma.ij.neg) / num.neg^2
  )
  # get sigma.eps.sq
  sigma.eps.sq = theta^2 * varIlrX * (1 / Rsq - 1)
  # return sigma.eps
  return(sqrt(sigma.eps.sq))
}

#####

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

#####

alrinv <- function(y) {
  x <- cbind(exp(y), 1)
  x / rowSums(x)
}

simulateLatentVarReg = function(
  theta1, coeffX1, sbp, n, sigma.noise
){
  # get info from sbp
  p = ncol(sbp)
  num.pos = length(which(SBP.true == 1))
  num.neg = length(which(SBP.true == -1))
  # get latent variable
  latent.var = runif(n)
  # simulate y from latent variable
  y = theta1 * latent.var + rnorm(n) * sigma.noise
  # simulate X
  alrX.noises = as.matrix(rnorm(2 * n * (p - 1)), nrow = 2 * n) * sigma.noise
  alrX.covariates = as.matrix(
    c(rep(1, num.pos) / num.pos, 
      -rep(1, num.neg) / num.neg, 
      rep(0, p - 1 - num.pos - num.neg)))
  alrX = as.matrix(latent.var) %*% coeffX1 * t(alrX.covariates) + alrX.noises
}
