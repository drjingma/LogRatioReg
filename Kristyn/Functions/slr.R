

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
spectral.clustering.kmeans = function(
  W, n_eig = 2, reindex = FALSE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01
  ) {
  # compute graph laplacian
  L = graph.laplacian2(
    W = W, amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter)          
  ei = eigen(L, symmetric = TRUE)
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
  return(cl)
}

graph.laplacian2 = function(
  W, # normalized = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01
){
  stopifnot(nrow(W) == ncol(W)) 
  n = nrow(W)    # number of vertices
  degrees <- colSums(W) # degrees of vertices
  W.tmp = W
  
  # Amini et al., 2016 regularization method: perturb the network by adding 
  #   some links with low edge weights
  if(amini.regularization){
    W.tmp <- W.tmp + 
      amini.regularization.parameter * mean(degrees) / n * tcrossprod(rep(1,n))
  }
  
  # if(normalized){
  D_half = diag(1 / sqrt(degrees))
  # } else {
  #   return(W.tmp)
  # }
  return(D_half %*% W.tmp %*% D_half)
}

# spectral clustering using sign check method
spectral.clustering.cut = function(
  W, reindex = FALSE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01
) {
  stopifnot(nrow(W) == ncol(W)) 
  n = nrow(W)    # number of vertices
  
  # compute graph laplacian
  L = graph.laplacian2(
    W = W, amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter)
  IminusL = diag(n) - L
  
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
  return(cl)
}

# slr -- 1 application of spectral clustering
slr <- function(
    x, y, classification = FALSE, approx = FALSE, 
    amini.regularization = FALSE, 
    amini.regularization.parameter = 0.01, 
    selection.crit = "Rsq", # "cor", "Rsq", "selbal"
    ad.hoc = FALSE
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
  cluster.labels <- spectral.clustering.kmeans(
    affinityMat, n_eig = num.clusters, 
    amini.regularization = amini.regularization,
    amini.regularization.parameter = amini.regularization.parameter)
  cluster.lengths = table(cluster.labels)
  clusters = names(cluster.lengths)
  out$num.clusters = num.clusters
  
  ## Find which pair of the 3 clusters has the most predictive balance
  all.pairs = t(combn(clusters, 2))
  cors = rep(NA, nrow(all.pairs))
  Rsqs = rep(NA, nrow(all.pairs))
  # pvals = rep(NA, nrow(all.pairs))
  # adjRsqs = rep(NA, nrow(all.pairs))
  # aics = rep(NA, nrow(all.pairs))
  # bics = rep(NA, nrow(all.pairs))
  sbp.ests = matrix(0, nrow = p, ncol = nrow(all.pairs))
  bal.ests = matrix(NA, nrow = n, ncol = nrow(all.pairs))
  rownames(sbp.ests) <- colnames(x)
  for(i in 1:nrow(all.pairs)){
    cluster1.label = all.pairs[i, 1]
    subset1 = which(cluster.labels == cluster1.label)
    cluster2.label = all.pairs[i, 2]
    subset2 = which(cluster.labels == cluster2.label)
    sbp.ests[match(names(subset1), rownames(sbp.ests)), i] = 1
    sbp.ests[match(names(subset2), rownames(sbp.ests)), i] = -1
    bal.ests[, i] = balance::balance.fromSBP(
      x = x, y = sbp.ests[, i, drop = FALSE])
    # cors and Rsqs
    cors[i] = stats::cor(bal.ests[, i], y)
    Rsqs[i] = cors[i]^2
    # # significance -- didn't fix balance selection problem!
    # refit_data.tmp = data.frame(V1 = bal.ests[, i], y = y)
    # if(!classification){
    #   refit.tmp <- lm(y~V1, data = refit_data.tmp)
    #   pvals[i] = summary(refit.tmp)$coefficients["V1", "Pr(>|t|)"]
    #   adjRsqs[i] = summary(refit.tmp)$adj.r.squared
    #   aics[i] = AIC(refit.tmp)
    #   bics[i] = BIC(refit.tmp)
    # } else{
    #   refit.tmp = stats::glm(
    #     y~V1, data = refit_data.tmp, family = binomial(link = "logit"))
    #   pvals[i] = summary(refit.tmp)$coefficients["V1", "Pr(>|z|)"]
    #   adjRsqs[i] = summary(refit.tmp)$adj.r.squared # is this a thing for glm?
    # }
  }
  out$cors = cors
  out$Rsqs = Rsqs
  # The correct active set should have largest correlation magnitude, i.e. Rsq.
  ## We refit the linear model on the balance from the set with the 
  ##    largest correlation magnitude.
  if(selection.crit %in% c("selbal")){
    # calculate geometric means of I+, I-, inactive set I0
    geometric.means = matrix(NA, nrow = nrow(x), ncol = num.clusters)
    colnames(geometric.means) = clusters
    for(i in 1:num.clusters){
      geometric.means[, i] =  apply(
        x[, which(cluster.labels == clusters[i])], 1,
        function(row) geometric.mean(row))
    }
    # fit selbal
    data_slbl = getSelbalData(
      X = geometric.means, y = y, classification = classification)
    fit_slbl = selbal::selbal.aux(
      x = data_slbl$X, y = data_slbl$y, logt = TRUE, maxV = 2)
    # get balance
    slbl.pair = sort(fit_slbl$Taxa) # assuming each pair in all.pairs is sorted
    selected.pair = which(
      apply(all.pairs, 1, function(row) isTRUE(all.equal(row, slbl.pair))))
    refit_data = data.frame(V1 = bal.ests[, selected.pair], y = y)
  } else{ # selection.crit %in% c("Rsq", "cor")
    selected.pair = which.max(Rsqs)
  }
  if(ad.hoc){ # choose the sparsest balance
    ad.hoc.invoked = FALSE
    # pick the balance with smallest active set
    cardinality.ests = apply(sbp.ests, 2, function(col) sum(col != 0))
    min.card.pair = which.min(cardinality.ests)
    if(cardinality.ests[selected.pair] != cardinality.ests[min.card.pair]){
      # save original sbp
      sbp.original = sbp.ests[, selected.pair, drop = FALSE]
      rownames(sbp.original) = colnames(x)
      out$sbp.original = sbp.original
      # new selected pair, via ad hoc method
      ad.hoc.invoked = TRUE
      selected.pair = min.card.pair
    }
    out$ad.hoc.invoked = ad.hoc.invoked
  }
  
  # fit the balance regression model
  out$index = sbp.ests[, selected.pair] # redundant, may remove
  refit_data = data.frame(V1 = bal.ests[, selected.pair], y = y)
  if(!classification){
    refit <- lm(y~V1, data = refit_data)
  } else{
    refit = stats::glm(
      y~V1, data = refit_data, family = binomial(link = "logit"))
  }
  out$model <- refit
  
  # return the SBP vector (with entries for all p variables)
  full.sbp.est = sbp.ests[, selected.pair, drop = FALSE]
  rownames(full.sbp.est) = colnames(x)
  out$sbp = full.sbp.est
  
  # return the slr object
  return(out)
}





# # used to be slr() -- two applications of spectral clustering
# slr2sc <- function(
#   x, y, num.clusters = 2, classification = FALSE, approx = FALSE, 
#   amini.regularization = FALSE, 
#   amini.regularization.parameter = 0.01, 
#   highdegree.regularization.summary = "mean",
#   highdegree.regularization = FALSE,
#   spectral.clustering.method = "kmeans" # "kmeans" or "cut"
# ){
#   if(spectral.clustering.method == "cut" & num.clusters != 2){
#     stop("spectral.clustering.method == cut requires num.clusters == 2")
#   }
#   
#   n = nrow(x)
#   p = ncol(x)
#   
#   ## Compute pairwise correlation
#   rhoMat0 <- slrmatrix(x = x, y = y)
#   rhoMat = max(rhoMat0) - rhoMat0
#   
#   out <- list()
#   out$kernel <- rhoMat
#   
#   ## Split into active/inactive sets
#   if(approx){
#     rhoMat.svd <- svd(rhoMat)
#     rhoMat_approx <- tcrossprod(
#       rhoMat.svd$u[, 1], 
#       rhoMat.svd$v[, 1]) *
#       rhoMat.svd$d[1]
#     rownames(rhoMat_approx) <- colnames(rhoMat_approx) <- rownames(rhoMat)
#     affinityMat = rhoMat_approx
#   } else{
#     affinityMat = rhoMat
#   }
#   if(spectral.clustering.method == "kmeans"){
#     clusters1 <- spectral.clustering.kmeans(
#       affinityMat, n_eig = num.clusters, 
#       amini.regularization = amini.regularization,
#       amini.regularization.parameter = amini.regularization.parameter,
#       highdegree.regularization.summary = highdegree.regularization.summary,
#       highdegree.regularization = highdegree.regularization)
#   } else if(spectral.clustering.method == "cut"){
#     clusters1 <- spectral.clustering.cut(
#       affinityMat, 
#       amini.regularization = amini.regularization,
#       amini.regularization.parameter = amini.regularization.parameter,
#       highdegree.regularization.summary = highdegree.regularization.summary,
#       highdegree.regularization = highdegree.regularization)
#   }
#   cluster.lengths = table(clusters1)
#   clusters = names(cluster.lengths)
#   out$num.clusters = num.clusters
#   
#   ## Find which set has the more predictive balance
#   cors = rep(NA, num.clusters)
#   Rsqs = rep(NA, num.clusters)
#   sbp.ests = matrix(0, nrow = p, ncol = num.clusters)
#   bal.ests = matrix(NA, nrow = n, ncol = num.clusters)
#   rownames(sbp.ests) <- colnames(x)
#   for(i in 1:num.clusters){ 
#     cluster.label = as.numeric(clusters[i])
#     index = which(clusters1 == cluster.label)
#     if(length(index) >= 2){ # if a log-ratio can be made from the variables in this cluster
#       ## Perform spectral clustering to get the numerator/denominator groups
#       if(spectral.clustering.method == "kmeans"){
#         subset = spectral.clustering.kmeans(
#           rhoMat[index, index], n_eig = 2, reindex = TRUE, 
#           amini.regularization = amini.regularization,
#           amini.regularization.parameter = amini.regularization.parameter,
#           highdegree.regularization.summary = highdegree.regularization.summary, 
#           highdegree.regularization = highdegree.regularization)
#       } else if(spectral.clustering.method == "cut"){
#         subset = spectral.clustering.cut(
#           rhoMat[index, index], reindex = TRUE, 
#           amini.regularization = amini.regularization,
#           amini.regularization.parameter = amini.regularization.parameter,
#           highdegree.regularization.summary = highdegree.regularization.summary, 
#           highdegree.regularization = highdegree.regularization)
#       }
#       ## calculate balance and its correlation with y
#       sbp.ests[match(names(subset), rownames(sbp.ests)), i] = subset
#       bal.ests[, i] = balance::balance.fromSBP(
#         x = x, y = sbp.ests[, i, drop = FALSE])
#     } else{ # otherwise, if there's just one variable in the cluster, set all other variables as -1
#       sbp.ests[, i] = rep(NA, p)
#     }
#     cors[i] = stats::cor(bal.ests[, i], y)
#     Rsqs[i] = cors[i]^2
#   }
#   out$cors = cors
#   out$Rsqs = Rsqs
#   # The correct active set should have largest correlation magnitude, i.e. Rsq.
#   ## We refit the linear model on the balance from the set with the 
#   ##    largest correlation magnitude.
#   selected.cluster = which.max(Rsqs)
#   out$index = sbp.ests[, selected.cluster]
#   refit_data = data.frame(V1 = bal.ests[, selected.cluster], y = y)
#   if(!classification){
#     refit <- lm(y~V1, data = refit_data)
#   } else{
#     refit = stats::glm(
#       y~V1, data = refit_data, family = binomial(link = "logit"))
#   }
#   out$model <- refit
#   
#   # return the full SBP vector (with entries for all p variables)
#   full.sbp.est = sbp.ests[, selected.cluster, drop = FALSE]
#   rownames(full.sbp.est) = colnames(x)
#   out$sbp = full.sbp.est
#   
#   # return the slr object
#   return(out)
# }


