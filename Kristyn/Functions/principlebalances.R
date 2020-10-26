fitPBLasso = function(X, y, lambda = NULL){
  # get principle balances
  sbp_pba = pba(X)
  # lasso_fit
  glm.temp = run.glmnet(sbp_pba@pba, y, lambda)
  return(list(
    betahat = as.vector(glm.temp$beta), 
    lambda = glm.temp$lamdba
  ))
}
