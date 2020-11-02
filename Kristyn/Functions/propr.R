fitproprLasso = function(X, y, linkage, lambda = NULL){
  
  # get propr tree
  pr <- propr::propr(X, metric = "phs")
  propr_tree = hclust(as.dist(pr@matrix),method = linkage)
  # compute balances
  sbp_propr = sbp.fromHclust(propr_tree)
  Xb = balance.fromSBP(X, sbp_propr)
  
  # lasso fit
  glm.temp = fitGlmnet(Xb, y, lambda)
  
  return(list(
    glmnet = glm.temp$glmnet,
    btree = propr_tree,
    betahat = as.vector(glm.temp$beta), 
    lambda = glm.temp$lamdba
  ))
}
