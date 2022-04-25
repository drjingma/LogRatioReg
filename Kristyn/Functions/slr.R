

slrmatrix = function(x, y){
  p = ncol(x)
  
  rhoMat = matrix(0, p, p)
  rownames(rhoMat) <- colnames(rhoMat) <- colnames(x)
  for (j in 1:p){
    for (k in 1:p){
      if (k==j){next}
      else {
        rhoMat[j,k] <- abs(stats::cor(log(x[,j])-log(x[,k]),y)) # correlation
      }
    }
  }
  return(rhoMat)
}

# spectral.clustering() function, with automatic splitting of two things
spectral.clustering2 = function(
  W, n_eig = 2, reindex = FALSE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01,
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
  ) {
  # compute graph laplacian
  L = graph.laplacian2(
    W = W, amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter,
    highdegree.regularization.summary = highdegree.regularization.summary,
    highdegree.regularization = highdegree.regularization)          
  ei = eigen(L, symmetric = TRUE)
  # compute the eigenvectors and values of L
  # we will use k-means to cluster the data
  # using the leading eigenvalues in absolute values
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  if(nrow(W) == n_eig){
    obj = list(cluster= 1:n_eig)
  } else{
    if(include.leading.eigenvector){
      eigenvector.indices = 1:n_eig
    } else{ # just the second largest eigenvector
      eigenvector.indices = 2:n_eig
    }
    obj <- kmeans(
      # ei$vectors[, 2:n_eig, drop = FALSE], centers = n_eig, nstart = 100)
      ei$vectors[, eigenvector.indices, drop = FALSE], centers = n_eig, 
      nstart = 100)
  }
  if (reindex){
    cl <- 2*(obj$cluster - 1) - 1
  } else{
    cl = obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl)
}

graph.laplacian2 = function(
  W, # normalized = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01,
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE
){
  stopifnot(nrow(W) == ncol(W)) 
  n = nrow(W)    # number of vertices
  degrees <- colSums(W) # degrees of vertices
  if(highdegree.regularization){
    if(highdegree.regularization.summary == "maximal"){
      maximaldegree = max(n * W) #* highdegree.regularization.parameter
    } else{
      maximaldegree = do.call(highdegree.regularization.summary, list(degrees))
    }
  }
  W.tmp = W
  
  # Amini et al., 2016 regularization method: perturb the network by adding 
  #   some links with low edge weights
  if(amini.regularization){
    W.tmp <- W.tmp + 
      amini.regularization.parameter * mean(degrees) / n * tcrossprod(rep(1,n))
  }
  
  # high-degree regularization: reduce the weights of edges proportionally to 
  #   the excess of degrees
  if(highdegree.regularization){
    lambdas = sapply(degrees, function(x) min(2 * maximaldegree / x, 1))
    for(i in 1:n){
      for(j in 1:n){
          W.tmp[i, j] = sqrt(lambdas[i] * lambdas[j]) * W.tmp[i, j]
      }
    }
  }
  
  # if(normalized){
  D_half = diag(1 / sqrt(degrees))
  # } else {
  #   return(W.tmp)
  # }
  return(D_half %*% W.tmp %*% D_half)
}

slr <- function(
  x, y, num.clusters = 2, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE, 
  subtractFrom1 = FALSE
){
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  rhoMat <- slrmatrix(x = x, y = y)
  if(subtractFrom1){
    rhoMat = 1 - rhoMat
  }
  out <- list()
  out$kernel <- rhoMat
  
  ## Split into active/inactive sets
  if(approx){
    rhoMat.svd <- svd(rhoMat)
    rhoMat_approx <- tcrossprod(
      rhoMat.svd$u[, 1], 
      rhoMat.svd$v[, 1]) *
      rhoMat.svd$d[1]
    rownames(rhoMat_approx) <- colnames(rhoMat_approx) <- rownames(rhoMat)
    clusters1 <- spectral.clustering2(
      rhoMat_approx, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  } else{
    clusters1 <- spectral.clustering2(
      rhoMat, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  }
  cluster.lengths = table(clusters1)
  clusters = names(cluster.lengths)
  out$num.clusters = num.clusters
  
  ## Find which set has the more predictive balance
  cors = rep(NA, num.clusters)
  Rsqs = rep(NA, num.clusters)
  sbp.ests = matrix(0, nrow = p, ncol = num.clusters)
  bal.ests = matrix(NA, nrow = n, ncol = num.clusters)
  rownames(sbp.ests) <- colnames(x)
  for(i in 1:num.clusters){ 
    cluster.label = as.numeric(clusters[i])
    index = which(clusters1 == cluster.label)
    if(length(index) >= 2){ # if a log-ratio can be made from the variables in this cluster
      ## Perform spectral clustering to get the numerator/denominator groups
      subset = spectral.clustering2(
        rhoMat[index, index], n_eig = 2, reindex = TRUE, 
        amini.regularization = amini.regularization,
        amini.regularization.parameter = amini.regularization.parameter,
        highdegree.regularization.summary = highdegree.regularization.summary, 
        highdegree.regularization = highdegree.regularization,
        include.leading.eigenvector = include.leading.eigenvector)
      ## calculate balance and its correlation with y
      sbp.ests[match(names(subset), rownames(sbp.ests)), i] = subset
      bal.ests[, i] = balance::balance.fromSBP(
        x = x, y = sbp.ests[, i, drop = FALSE])
    } else{ # otherwise, if there's just one variable in the cluster, set all other variables as -1
      sbp.ests[, i] = rep(NA, p)
    }
    cors[i] = stats::cor(bal.ests[, i], y)
    Rsqs[i] = cors[i]^2
  }
  out$cors = cors
  out$Rsqs = Rsqs
  # The correct active set should have largest correlation magnitude, i.e. Rsq.
  ## We refit the linear model on the balance from the set with the 
  ##    largest correlation magnitude.
  selected.cluster = which.max(Rsqs)
  out$index = sbp.ests[, selected.cluster]
  refit_data = data.frame(V1 = bal.ests[, selected.cluster], y = y)
  if(!classification){
    refit <- lm(y~V1, data = refit_data)
  } else{
    refit = stats::glm(
      y~V1, data = refit_data, family = binomial(link = "logit"))
  }
  out$model <- refit
  
  # return the full SBP vector (with entries for all p variables)
  full.sbp.est = sbp.ests[, selected.cluster, drop = FALSE]
  rownames(full.sbp.est) = colnames(x)
  out$sbp = full.sbp.est
  
  # return the slr object
  return(out)
}



























# # slr() function with the following changes:
# #   1. rank1approx option added
# #   2. using balances correlations instead of effect sizes to choose active subset
# #   3. spectral.clustering2() instead of spectral.clustering()
# slrT2 <- function(
#   x, y, rank1approx = FALSE, classification = FALSE, alpha = 0.05){
#   p <- ncol(x)
#   
#   ## Compute pairwise correlation 
#   rhoMat <- slrmatrix(x = x, y = y)
#   out <- list()
#   out$kernel <- rhoMat
#   
#   ## Split into active/inactive sets
#   if(rank1approx){
#     rhoMat.svd <- svd(rhoMat)
#     rhoMat_approx_1 <- tcrossprod(
#       rhoMat.svd$u[, 1], rhoMat.svd$v[, 1]) * rhoMat.svd$d[1]
#     index <- which(spectral.clustering2(rhoMat_approx_1, reindex = TRUE) == 1)
#   } else{
#     index <- which(spectral.clustering2(rhoMat, reindex = TRUE) == 1)
#   }
#   
#   ## Find which set has the more predictive balance
#   # if(length(index) == 1 | length(index) == p - 1){ # one valid cluster
#   #   warning("slr(): Only one covariate was subsetted. No log-ratios can be made.")
#   #   # fit the subset with length p - 1
#   #   if(length(index) == p - 1){
#   #     subset <- spectral.clustering2(rhoMat[index, index], reindex = TRUE)
#   #   } else{
#   #     subset <- spectral.clustering2(rhoMat[-index, -index], reindex = TRUE)
#   #   }
#   #   sbp.est <- matrix(0, nrow = p, ncol = 1)
#   #   rownames(sbp.est) <- colnames(x)
#   #   sbp.est[match(names(subset),rownames(sbp.est)), 1] <- subset
#   #   est.balance <- balance::balance.fromSBP(x = x, y = sbp.est)
#   #   if(!classification){
#   #     refit2 <- lm(y~est.balance)
#   #     pval = summary(refit2)$coefficients["est.balance", "Pr(>|t|)"]
#   #   } else{
#   #     refit2 = stats::glm(
#   #       y~est.balance, family = binomial(link = "logit"))
#   #     pval = summary(refit2)$coefficients["est.balance", "Pr(>|z|)"]
#   #   }
#   #   if(pval < alpha){
#   #     out$index = subset
#   #     refit = refit2
#   #   } else{
#   #     # fit only the intercept
#   #     out$index = NA
#   #     if(!classification){
#   #       refit <- lm(y~1)
#   #     } else{
#   #       refit = stats::glm(y~1, family = binomial(link = "logit"))
#   #     }
#   #   }
#   # } else{
#   ## Perform spectral clustering to get the numerator/denominator groups
#   ##    for the active/inactive sets
#   subset1 <- spectral.clustering2(rhoMat[index, index], reindex = TRUE)
#   subset2 <- spectral.clustering2(rhoMat[-index, -index], reindex = TRUE) 
#   
#   # Since we don't know which subset contains the active variables, 
#   ## we first fit a linear model with balances obtained from both subsets. 
#   sbp.est <- matrix(0, ncol = 2, nrow = p)
#   rownames(sbp.est) <- colnames(x)
#   sbp.est[match(names(subset1),rownames(sbp.est)), 1] <- subset1
#   sbp.est[match(names(subset2),rownames(sbp.est)), 2] <- subset2
#   est.balance <- balance::balance.fromSBP(x = x, y = sbp.est)
#   cors = stats::cor(est.balance, y)
#   
#   # The correct subset should have larger coefficient. 
#   ## We refit the linear model on the balance from the correct subset. 
#   out$cors = cors
#   
#   if ( abs(cors[1, ]) > abs(cors[2, ]) ){
#     # pick subset1
#     out$index <- subset1
#     if(!classification){
#       refit <- lm(y~est.balance[, 1])
#     } else{
#       refit = stats::glm(
#         y~est.balance[, 1], family = binomial(link = "logit"))
#     }
#   } else {
#     # pick subset2
#     out$index <- subset2
#     if(!classification){
#       refit <- lm(y~est.balance[, 2])
#     } else{
#       refit = stats::glm(
#         y~est.balance[, 2], family = binomial(link = "logit"))
#     }
#   }
#   # }
#   
#   out$model <- refit
#   out
# }

