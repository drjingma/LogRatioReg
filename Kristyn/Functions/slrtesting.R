# spectral.clustering() function, with automatic splitting of two things
spectral.clustering.kmeans_testing = function(
  W, n_eig = 2, reindex = FALSE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01,
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE
) {
  # compute graph laplacian
  L = graph.laplacian2_testing(
    W = W, amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter,
    highdegree.regularization.summary = highdegree.regularization.summary,
    highdegree.regularization = highdegree.regularization)          
  ei = eigen(L$L, symmetric = TRUE)
  # compute the eigenvectors and values of L
  # we will use k-means to cluster the data
  # using the leading eigenvalues in absolute values
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  if(nrow(W) == n_eig){
    obj = list(cluster= 1:n_eig)
  } else{
    obj <- kmeans(
      ei$vectors[, 1:n_eig, drop = FALSE], centers = n_eig, 
      nstart = 100)
  }
  if (reindex){
    cl <- 2*(obj$cluster - 1) - 1
  } else{
    cl = obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(list(cl = cl, L = L, ei = ei))
}

graph.laplacian2_testing = function(
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
    weights = matrix(NA, nrow = n, ncol = n)
    for(i in 1:n){
      for(j in 1:n){
        weightij = sqrt(lambdas[i] * lambdas[j])
        weights[i, j] = weightij
        W.tmp[i, j] = weightij * W.tmp[i, j]
      }
    }
  }

  D_half = diag(1 / sqrt(degrees)) # Normalize
  L = D_half %*% W.tmp %*% D_half # Laplacian
  return_obj = list(
    L = L, W = W, W.tmp = W.tmp,
    amini.regularization = amini.regularization,
    amini.regularization.parameter = amini.regularization.parameter,
    highdegree.regularization = highdegree.regularization,
    highdegree.regularization.summary = highdegree.regularization.summary
  )
  if(highdegree.regularization){
    return_obj$lambdas = lambdas
    return_obj$weights = weights
  }
  return(return_obj)
}

# spectral clustering using sign check method
spectral.clustering.cut_testing = function(
  W, reindex = FALSE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01,
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE
) {
  stopifnot(nrow(W) == ncol(W)) 
  n = nrow(W)    # number of vertices
  
  # compute graph laplacian
  L = graph.laplacian2_testing(
    W = W, amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter,
    highdegree.regularization.summary = highdegree.regularization.summary,
    highdegree.regularization = highdegree.regularization)
  IminusL = diag(n) - L$L
  
  ei = eigen(IminusL, symmetric = TRUE)
  # compute the eigenvectors and values of L
  # we will use k-means to cluster the data
  # using the leading eigenvalues in absolute values
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing = TRUE)]
  if(nrow(W) == 2){
    obj = list(cluster = 1:2)
  } else{
    eigenvector2 = ei$vectors[, 2, drop = FALSE]
    obj  = list(cluster = as.numeric(eigenvector2 < 0) + 1)
  }
  if (reindex){
    cl <- 2*(obj$cluster - 1) - 1
  } else{
    cl = obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(list(cl = cl, L = L, ei = ei))
}

slr2sc_testing <- function(
  x, y, num.clusters = 2, classification = FALSE, approx = FALSE, 
  amini.regularization = FALSE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  similarity.matrix = TRUE, 
  maxGamma = FALSE,
  spectral.clustering.method = "kmeans" # "kmeans" or "cut"
){
  if(spectral.clustering.method == "cut" & num.clusters != 2){
    stop("spectral.clustering.method == cut requires num.clusters == 2")
  }
  
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  rhoMat <- slrmatrix(x = x, y = y)
  if(similarity.matrix){
    if(maxGamma){
      rhoMat = max(rhoMat) - rhoMat
    } else{
      rhoMat = 1 - rhoMat
    }
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
    affinityMat = rhoMat_approx
  } else{
    affinityMat = rhoMat
  }
  if(spectral.clustering.method == "kmeans"){
    clusters1 <- spectral.clustering.kmeans_testing(
      affinityMat, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization)
  } else if(spectral.clustering.method == "cut"){
    clusters1 <- spectral.clustering.cut_testing(
      affinityMat, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization)
  }
  cluster.lengths = table(clusters1$cl)
  out$spectralclustering1 = clusters1
  clusters = names(cluster.lengths)
  out$num.clusters = num.clusters
  
  ## Find which set has the more predictive balance
  cors = rep(NA, num.clusters)
  Rsqs = rep(NA, num.clusters)
  sbp.ests = matrix(0, nrow = p, ncol = num.clusters)
  bal.ests = matrix(NA, nrow = n, ncol = num.clusters)
  rownames(sbp.ests) <- colnames(x)
  out$spectralclustering2 = list()
  for(i in 1:num.clusters){ 
    cluster.label = as.numeric(clusters[i])
    index = which(clusters1$cl == cluster.label)
    if(length(index) >= 2){ # if a log-ratio can be made from the variables in this cluster
      ## Perform spectral clustering to get the numerator/denominator groups
      if(spectral.clustering.method == "kmeans"){
        subset = spectral.clustering.kmeans_testing(
          rhoMat[index, index], n_eig = 2, reindex = TRUE, 
          amini.regularization = amini.regularization,
          amini.regularization.parameter = amini.regularization.parameter,
          highdegree.regularization.summary = highdegree.regularization.summary, 
          highdegree.regularization = highdegree.regularization)
      } else if(spectral.clustering.method == "cut"){
        subset = spectral.clustering.cut_testing(
          rhoMat[index, index], reindex = TRUE, 
          amini.regularization = amini.regularization,
          amini.regularization.parameter = amini.regularization.parameter,
          highdegree.regularization.summary = highdegree.regularization.summary, 
          highdegree.regularization = highdegree.regularization)
      }
      out$spectralclustering2[[i]] = subset
      ## calculate balance and its correlation with y
      sbp.ests[match(names(subset$cl), rownames(sbp.ests)), i] = subset$cl
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

slr1sc <- function(
    x, y, classification = FALSE, approx = FALSE, 
    amini.regularization = FALSE, 
    amini.regularization.parameter = 0.01, 
    highdegree.regularization.summary = "mean",
    highdegree.regularization = FALSE
){
  
  num.clusters = 3
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  rhoMat0 <- slrmatrix(x = x, y = y)
  rhoMat = max(rhoMat0) - rhoMat0
  
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
    affinityMat = rhoMat_approx
  } else{
    affinityMat = rhoMat
  }
  cluster.labels <- spectral.clustering.kmeans_testing(
    affinityMat, n_eig = num.clusters, 
    amini.regularization = amini.regularization,
    amini.regularization.parameter = amini.regularization.parameter,
    highdegree.regularization.summary = highdegree.regularization.summary,
    highdegree.regularization = highdegree.regularization)
  cluster.lengths = table(cluster.labels$cl)
  clusters = names(cluster.lengths)
  out$num.clusters = num.clusters
  
  ## Find which pair of the 3 clusters has the most predictive balance
  all.pairs = t(combn(clusters, 2))
  cors = rep(NA, nrow(all.pairs))
  Rsqs = rep(NA, nrow(all.pairs))
  sbp.ests = matrix(0, nrow = p, ncol = nrow(all.pairs))
  bal.ests = matrix(NA, nrow = n, ncol = nrow(all.pairs))
  rownames(sbp.ests) <- colnames(x)
  for(i in 1:nrow(all.pairs)){
    cluster1.label = all.pairs[i, 1]
    subset1 = which(cluster.labels$cl == cluster1.label)
    cluster2.label = all.pairs[i, 2]
    subset2 = which(cluster.labels$cl == cluster2.label)
    sbp.ests[match(names(subset1), rownames(sbp.ests)), i] = 1
    sbp.ests[match(names(subset2), rownames(sbp.ests)), i] = -1
    bal.ests[, i] = balance::balance.fromSBP(
      x = x, y = sbp.ests[, i, drop = FALSE])
    cors[i] = stats::cor(bal.ests[, i], y)
    Rsqs[i] = cors[i]^2
  }
  out$cors = cors
  out$Rsqs = Rsqs
  # The correct active set should have largest correlation magnitude, i.e. Rsq.
  ## We refit the linear model on the balance from the set with the 
  ##    largest correlation magnitude.
  selected.bal = which.max(Rsqs)
  out$index = sbp.ests[, selected.bal]
  refit_data = data.frame(V1 = bal.ests[, selected.bal], y = y)
  if(!classification){
    refit <- lm(y~V1, data = refit_data)
  } else{
    refit = stats::glm(
      y~V1, data = refit_data, family = binomial(link = "logit"))
  }
  out$model <- refit
  
  # return the full SBP vector (with entries for all p variables)
  full.sbp.est = sbp.ests[, selected.bal, drop = FALSE]
  rownames(full.sbp.est) = colnames(x)
  out$sbp = full.sbp.est
  
  # return the slr object
  return(out)
}
