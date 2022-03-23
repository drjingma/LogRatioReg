getLogLik = function(residuals, n = NULL, weights = NULL){
  if(is.null(n)) n = length(residuals)
  if(is.null(weights)) weights = rep(1, n)
  ll = 0.5 * (
    sum(log(weights)) - 
      n * (log(2 * pi) + 1 - log(n) + log(sum(weights * residuals^2)))
  )
  return(ll)
}
getBIC = function(num.parameters, residuals, n = NULL, weights = NULL){
  # https://stackoverflow.com/questions/35131450/calculating-bic-manually-for-lm-object
  if(is.null(n)) n = length(residuals)
  if(is.null(weights)) weights = rep(1, n)
  df.ll = num.parameters + 1
  ll = getLogLik(residuals = residuals, n = n, weights = weights)
  bic = -2 * ll + log(n) * df.ll
  return(bic)
}

getBICseq = function(
  y, predMat, betahat0.vec, betahat.mat
){
  bics = rep(NA, length(betahat0.vec))
  for(i in 1:length(betahat0.vec)){
    a0.i = betahat0.vec[i]
    bet.i = betahat.mat[, i]
    yhat.i = a0.i + predMat %*% bet.i
    resid.i = y - yhat.i
    num.params.i = sum(abs(bet.i) > 10e-8)
    bics[i] = getBIC(num.params.i, resid.i)
  }
  return(bics)
}
