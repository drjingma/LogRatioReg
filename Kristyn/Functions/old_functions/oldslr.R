# Supervised Log Ratios Regression

# do hierarchical clustering with specified linkage 
#   for compositional data X and y
getSlrMatrix = function(y, X, type = "similarity"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("getSlrMatrix(): dim(X)[1] != length(y)!!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(0, p, p) # diagonal == 1
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      val = abs(stats::cor(log(X[, j]) - log(X[, k]), y))
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # give the rows and columns the names of taxa in X, for sbp.fromHclust()
  rownames(cormat) = colnames(X)
  colnames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  if(type != "similarity"){
    cormat = 1 - cormat
  } 
  return(cormat)
}

# getSupervisedTree = function(
#   y, X, linkage = "complete", rho.type = "square"
# ){
#   Gammamat = getSupervisedMatrix(
#     y = y, X = X, rho.type = rho.type, type = "distance")
#   # get tree from hierarchical clustering
#   btree_slr = hclust(as.dist(Gammamat), method = linkage)
#   return(btree_slr)
# }
# # old version, with noise arguments
# getSupervisedTree.old = function(
#   y, X, linkage = "complete",
#   allow.noise = FALSE, noise = 1e-12
# ){
#   # browser()
#   n = dim(X)[1]
#   p = dim(X)[2]
#   if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
#   
#   # calculate correlation of each pair of log-ratios with response y
#   cormat = matrix(1, p, p) # diagonal == 1
#   y_demeaned = y - mean(y)
#   for (j in 1:(p - 1)){
#     for (k in (j + 1):p){
#       Zjk = log(X[, j]) - log(X[, k])
#       Zjk_demeaned = Zjk - mean(Zjk)
#       if(all(Zjk == 0)){ # add noise (hopefully never necessary! hopefully never in this situation)
#         if(allow.noise) Zjk = Zjk + rmvnorm(n, rep(0, n), noise * diag(n))
#         # do not do anything, except maybe warn
#       } else{
#         val = abs(cor(Zjk_demeaned, y_demeaned))
#       }
#       cormat[j, k] = val
#       cormat[k, j] = val
#     }
#   }
#   # find out which columns give na
#   # for(i in 1:p) if(!all(!is.na(cormat[, i]))) print(paste0(i, ":", which(is.na(cormat[, i]))))
#   #
#   # give the rows and columns the names of taxa in X,
#   #   so that sbp.fromHclust works later
#   rownames(cormat) = colnames(X)
#   rownames(cormat) = colnames(X)
#   
#   # get dissimilarity matrix
#   Gammamat = 1 - cormat
#   
#   # get tree from hierarchical clustering
#   btree_slr = hclust(as.dist(Gammamat), method = linkage)
#   
#   return(btree_slr)
# }
