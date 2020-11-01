fitPBLasso = function(X, y, lambda = NULL){
  # get principle balances
  sbp_pba = pba(X)
  # lasso_fit
  Xb = sbp_pba@pba
  glm.temp = fitGlmnet(Xb, y, lambda)
  
  return(list(
    glmnet = glm.temp$glmnet,
    betahat = as.vector(glm.temp$beta), 
    lambda = glm.temp$lamdba
  ))
}
