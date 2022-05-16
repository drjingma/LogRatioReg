
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
spectral.clustering.kmeans = function(W, n_eig = 2, reindex = FALSE) {
  # compute graph laplacian
  L = graph.laplacian2(W = W)          
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

graph.laplacian2 = function(W){
  stopifnot(nrow(W) == ncol(W)) 
  n = nrow(W)    # number of vertices
  degrees <- colSums(W) # degrees of vertices
  D_half = diag(1 / sqrt(degrees)) # normalization
  return(D_half %*% W %*% D_half) # return normalized Laplacian
}

# spectral clustering using sign check method
spectral.clustering.cut = function(W, reindex = FALSE) {
  stopifnot(nrow(W) == ncol(W)) 
  n = nrow(W)    # number of vertices
  
  # compute graph laplacian
  L = graph.laplacian2(W = W)
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
    x, y, classification = FALSE, 
    alpha = 0.05, ordered.tests = TRUE
){
  
  num.clusters = 3
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  rhoMat <- slrmatrix(x = x, y = y)
  Ahat = max(rhoMat) - rhoMat
  
  out <- list()
  out$kernel <- Ahat
  
  ## cluster
  cluster.labels <- spectral.clustering.kmeans(Ahat, n_eig = num.clusters)
  cluster.lengths = table(cluster.labels)
  clusters = names(cluster.lengths)
  out$num.clusters = num.clusters
  
  ## Find which pair of the 3 clusters produces the most predictive balance
  all.pairs = t(combn(clusters, 2))
  cors = rep(NA, nrow(all.pairs))
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
    cors[i] = stats::cor(bal.ests[, i], y)
  }
  out$cors = cors
  # balances' active set size
  cardinality.ests = apply(sbp.ests, 2, function(col) sum(col != 0))
  
  # The correct active set should have largest correlation magnitude, i.e. Rsq.
  ## We refit the linear model on the balance from the set with the 
  ##    largest correlation magnitude.
  max.pair = which.max(abs(cors))
  selected.pair = max.pair
  z.scores = 0.5 * log((1 + abs(cors)) / (1 - abs(cors)))
  sigma.z1.minus.z2 = sqrt(2 / (n - 3))
  if(ordered.tests){
    # rank the absolute correlations and comptue their z.scores
    ranks = rank(abs(cors))
    # ordered tests, where |r1| >= |r2| >= |r3|
    test.stats = c(
      z.scores[ranks[1]] - z.scores[ranks[2]], # test(|r1|, |r2|)
      z.scores[ranks[2]] - z.scores[ranks[3]], # test(|r2|, |r3|)
      z.scores[ranks[1]] - z.scores[ranks[3]] # test(|r1|, |r3|)
    ) / sigma.z1.minus.z2
    p.values = rep(NA, length(test.stats))
    for(i in 1:length(test.stats)){ # since lower.tail doesn't take a vector
      p.values[i] = pnorm(q = test.stats[i], lower.tail = test.stats[i] < 0)
    }
    signif.diffs = p.values <= alpha
    # if |r1| & |r2| are different, choose B1 (with max |r|, i.e. |r1|).
    if(!signif.diffs[1]){ 
      # if |r1| & |r2| are similar, check if |r2| & |r3| are similar
      if(signif.diffs[2]){
        # if |r1| & |r2| are similar, but |r2| & |r3| are different,
        #   choose the sparsest among B1 & B2
        if(cardinality.ests[ranks[2]] < cardinality.ests[ranks[1]]){
          selected.pair = ranks[2]
        }
      } else {
        # if |r1| & |r2| are similar, and |r2| & |r3| are similar, 
        #   check if |r1| & |r3| are similar --
        #   -- if |r1| & |r3| are different, choose the sparsest among B1 & B2
        if(cardinality.ests[ranks[2]] < cardinality.ests[ranks[1]]){
          selected.pair = ranks[2]
        }
        if(!signif.diffs[3]){
          # if |r1| & |r2| are similar, and |r2| & |r3| are similar, 
          #   and also |r1| & |r3| are different, 
          #   choose the sparsest among B1, B2, & B3
          if(cardinality.ests[ranks[3]] < cardinality.ests[ranks[2]]){
            selected.pair = ranks[3]
          }
        }
      }
    }
  } else{
    other.pairs = (1:3)[-max.pair]
    # tests (I), (II), (III)
    test.stats = c(
      z.scores[other.pairs] - z.scores[max.pair], 
      diff(z.scores[other.pairs])
    ) / sigma.z1.minus.z2
    p.values = rep(NA, length(test.stats))
    for(i in 1:length(test.stats)){ # since lower.tail doesn't take a vector
      p.values[i] = pnorm(q = test.stats[i], lower.tail = test.stats[i] < 0)
    }
    signif.diffs = p.values <= alpha
    if(all(signif.diffs[c(1, 2)])){ # if (I) & (II) == 1, 
      # choose Bmax
    } else if(all(!(signif.diffs[c(1, 2)]))){ # if(I) & (II) == 0, 
      # choose sparsest of the 3 balances
      selected.pair = which.min(cardinality.ests)
    } else if(sum(signif.diffs[c(1, 2)]) == 1){ # if (I) xor (II) == 1 (other 0),
      # check which is bigger, Bmax or the other balance
      # (either B1 or B2)
      sim.other.pair = other.pairs[!signif.diffs[c(1, 2)]]
      if(cardinality.ests[sim.other.pair] < cardinality.ests[max.pair]){
        selected.pair = sim.other.pair
        if(!signif.diffs[3]){ # if (III) == 0, 
          notsim.other.pair = other.pairs[signif.diffs[c(1, 2)]]
          # check which is bigger between B1 and B2
          if(cardinality.ests[notsim.other.pair] < 
             cardinality.ests[sim.other.pair]){
            selected.pair = notsim.other.pair
          }
        }
      }
    }
  }
  out$adhoc.invoked = (selected.pair == max.pair)
  
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



