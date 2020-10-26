fitCOATLasso = function(X, y, linkage, lambda = NULL){
  
  # get COAT tree
  d_coat = 1 - coat(X)$corr
  rownames(d_coat) = colnames(X)
  colnames(d_coat) = colnames(X)
  coat_tree = hclust(as.dist(d_coat),method = linkage)
  # compute balances
  sbp_coat = sbp.fromHclust(coat_tree)
  Xb = balance.fromSBP(X, sbp_coat)
  
  # lasso fit
  glm.temp = fitGlmnet(Xb, y, lambda)
  
  return(list(
    btree = coat_tree,
    betahat = as.vector(glm.temp$beta), 
    lambda = glm.temp$lamdba
  ))
}
