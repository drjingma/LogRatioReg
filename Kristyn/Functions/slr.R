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
  include.leading.eigenvector = FALSE
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
        # if(i != j){
          W.tmp[i, j] = sqrt(lambdas[i] * lambdas[j]) * W.tmp[i, j]
        # }
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
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  rhoMat <- slrmatrix(x = x, y = y)
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

# slr with 1 application of spectral clustering (with 3 clusters instead of 2)
slr1sc <- function(
  x, y, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  num.clusters = 3
  n = nrow(x)
  p = ncol(x)
  
  ## Compute pairwise correlation
  rhoMat <- slrmatrix(x = x, y = y)
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
    cluster.labels <- spectral.clustering2(
      rhoMat_approx, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary, 
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  } else{
    cluster.labels <- spectral.clustering2(
      rhoMat, n_eig = num.clusters, 
      amini.regularization = amini.regularization,
      amini.regularization.parameter = amini.regularization.parameter,
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  }
  cluster.lengths = table(cluster.labels)
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
    subset1 = which(cluster.labels == cluster1.label)
    cluster2.label = all.pairs[i, 2]
    subset2 = which(cluster.labels == cluster2.label)
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

cv.slr = function(
  x, y, max.clusters = 5, nfolds = 5, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  numclusters_candidates = 1:max.clusters
  slrmodels = list()
  for(i in numclusters_candidates){
    slrmodels[[i]] = slr(
      x = x, y = y, num.clusters = i, classification = classification, 
      approx = approx, 
      amini.regularization = amini.regularization, 
      amini.regularization.parameter = amini.regularization.parameter, 
      highdegree.regularization.summary = highdegree.regularization.summary,
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
  }
  
  # split the data into nfolds folds
  shuffle = sample(1:n)
  idfold = (shuffle %% nfolds) + 1
  n_fold = table(idfold)
  
  # calculate lasso for each fold removed
  cvm_sqerr = matrix(NA, nfolds, max.clusters) # squared error for each fold
  for(j in 1:nfolds){
    # Training data
    xtr = x[idfold != j, ]
    ytr = y[idfold != j]
    # Test data
    xte = x[idfold == j, ]
    yte = y[idfold == j]
    
    # for num.clusters = 1, ..., max.clusters, calculate squared error
    for(m in numclusters_candidates){
      fit_jm = slr(
        x = xtr, y = ytr, num.clusters = m, classification = classification, 
        approx = approx, 
        amini.regularization = amini.regularization, 
        amini.regularization.parameter = amini.regularization.parameter, 
        highdegree.regularization.summary = highdegree.regularization.summary,
        highdegree.regularization = highdegree.regularization,
        include.leading.eigenvector = include.leading.eigenvector)
      if(!classification){
        fit_jm.coefs = coefficients(fit_jm$model)
        ypred = fit_jm.coefs[1] + fit_jm.coefs[-1] * 
          balance::balance.fromSBP(x = xte, y = fit_jm$sbp)
        cvm_sqerr[j, m] = mean(crossprod(yte - ypred))
      } else{
        # how to do for classification? use auc? accuracy? see selbal.cv #################################################
        ypred = predict.glm(
          fit_jm$model,
          newdata = data.frame(
            V1 = balance::balance.fromSBP(x = xte, y = fit_jm$sbp)),
          type = "response")
        cvm_sqerr[j, m] = sum(crossprod(yte - ypred))
      }
    }
  }
  # Calculate CV(numclusters) and SE_CV(numclusters) for each possible 
  #   number of clusters, up to max.clusters
  # scores = -sqrt(cvm_sqerr) # codacore:::findBestCutoff.CoDaBaseLearner
  cvse = apply(cvm_sqerr, 2, stats::sd)/sqrt(nfolds)
  cvm = colMeans(cvm_sqerr)
  
  # Find numclust_min = argmin{CV(numclusters)}
  numclust_min_index = which.min(cvm)
  numclust_min = numclusters_candidates[numclust_min_index]
  
  # Find numclust_1se = maximal numclusters s.t. CV(numclusters) <= CV(nclust_min) + CV_SE(nclust_min)
  oneserule = cvm[numclust_min_index] + cvse[numclust_min_index]
  which_numclust_1se = which(cvm <= oneserule)
  numclust_1se_index = which_numclust_1se[length(which_numclust_1se)]
  numclust_1se = numclusters_candidates[numclust_1se_index]
  return(
    list(
      nclusters = numclusters_candidates,
      max.clusters = max.clusters, 
      models = slrmodels,
      nclusters_min = numclust_min, 
      nclusters_1se = numclust_1se, 
      cvm = cvm, 
      nclusters_min_idx = numclust_min_index,
      cvse = cvse,
      nclusters_1se_idx = numclust_1se_index
    )
  )
}

hslr <- function(
  x, y, num.levels = 1, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  hslrmodels = list()
  Rsqs = matrix(NA, nrow = num.levels, ncol = 2)
  cors = matrix(NA, nrow = num.levels, ncol = 2)
  vars.cur = colnames(x)
  true.num.levels = num.levels
  for(i in 1:num.levels){
    fit_i = slr(
      x = x[, vars.cur, drop = FALSE], y = y, num.clusters = 2, 
      classification = classification, approx = approx, 
      amini.regularization = amini.regularization, 
      amini.regularization.parameter = amini.regularization.parameter, 
      highdegree.regularization.summary = highdegree.regularization.summary, 
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
    # recalculate full SBP
    full.sbp.est = matrix(0, nrow = p, ncol = 1)
    rownames(full.sbp.est) = colnames(x)
    full.sbp.est[vars.cur, ] = c(fit_i$sbp)
    fit_i$sbp = full.sbp.est
    # save the fit
    hslrmodels[[i]] = fit_i
    vars.cur = rownames(fit_i$sbp)[fit_i$sbp != 0]
    cors[i, ] = fit_i$cors
    Rsqs[i, ] = fit_i$Rsqs
    if(length(vars.cur) <= 2){
      true.num.levels = i
      # warning(paste0("asked for num.levels = ", num.levels, " splits, but could only make up to ", i, " splits."))
      break
    }
  }
  
  # return the slr object
  return(list(
    models = hslrmodels, 
    cors = cors, 
    Rsqs = Rsqs, 
    true.num.levels = true.num.levels,
    num.levels = num.levels
  ))
}

hslr1sc <- function(
  x, y, num.levels = 1, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  hslrmodels = list()
  Rsqs = matrix(NA, nrow = num.levels, ncol = 3)
  cors = matrix(NA, nrow = num.levels, ncol = 3)
  vars.cur = colnames(x)
  true.num.levels = num.levels
  for(i in 1:num.levels){
    fit_i = slr1sc(
      x = x[, vars.cur, drop = FALSE], 
      y = y, classification = classification, 
      approx = approx, 
      amini.regularization = amini.regularization, 
      amini.regularization.parameter = amini.regularization.parameter, 
      highdegree.regularization.summary = highdegree.regularization.summary, 
      highdegree.regularization = highdegree.regularization,
      include.leading.eigenvector = include.leading.eigenvector)
    # recalculate full SBP
    full.sbp.est = matrix(0, nrow = p, ncol = 1)
    rownames(full.sbp.est) = colnames(x)
    full.sbp.est[vars.cur, ] = c(fit_i$sbp)
    fit_i$sbp = full.sbp.est
    # save the fit
    hslrmodels[[i]] = fit_i
    vars.cur = rownames(fit_i$sbp)[fit_i$sbp != 0]
    cors[i, ] = fit_i$cors
    Rsqs[i, ] = fit_i$Rsqs
    if(length(vars.cur) <= 2){
      true.num.levels = i
      # warning(paste0("asked for num.levels = ", num.levels, " splits, but could only make up to ", i, " splits."))
      break
    }
  }
  
  # return the slr object
  return(list(
    models = hslrmodels, 
    cors = cors, 
    Rsqs = Rsqs, 
    true.num.levels = true.num.levels,
    num.levels = num.levels
  ))
}

cv.hslr = function(
  x, y, max.levels = 5, nfolds = 5, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  hslrfit = hslr(
    x = x, y = y, num.levels = max.levels, classification = classification, 
    approx = approx, 
    amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter, 
    highdegree.regularization.summary = highdegree.regularization.summary, 
    highdegree.regularization = highdegree.regularization,
    include.leading.eigenvector = include.leading.eigenvector)
  num.levels.candidates = 1:hslrfit$true.num.levels
  hslrmodels = hslrfit$models
  
  # split the data into nfolds folds
  shuffle = sample(1:n)
  idfold = (shuffle %% nfolds) + 1
  n_fold = table(idfold)
  
  # calculate lasso for each fold removed
  cvm_sqerr = matrix(NA, nfolds, hslrfit$true.num.levels) # squared error for each fold
  for(j in 1:nfolds){
    # Training data
    xtr = x[idfold != j, ]
    ytr = y[idfold != j]
    # Test data
    xte = x[idfold == j, ]
    yte = y[idfold == j]
    
    fit_j = hslr(
      x = xtr, y = ytr, num.levels = hslrfit$true.num.levels, 
      classification = classification, approx = FALSE)
    
    # for num.clusters = 1, ..., max.levels, calculate squared error
    for(m in num.levels.candidates){
      if(length(fit_j$models) >= m){
        fit_jm = fit_j$models[[m]]
        if(!classification){
          fit_jm.coefs = coefficients(fit_jm$model)
          ypred = fit_jm.coefs[1] + fit_jm.coefs[-1] * 
            balance::balance.fromSBP(x = xte, y = fit_jm$sbp)
          cvm_sqerr[j, m] = mean(crossprod(yte - ypred))
        } else{
          # how to do for classification? use auc? accuracy? see selbal.cv #################################################
          ypred = predict.glm(
            fit_jm$model,
            newdata = data.frame(
              V1 = balance::balance.fromSBP(x = xte, y = fit_jm$sbp)),
            type = "response")
          cvm_sqerr[j, m] = sum(crossprod(yte - ypred))
        }
      } else{
        cvm_sqerr[j, m] = NA
        # warning("cv.hslr(): cvm_sqerr has NA values.")
      }
    }
  }
  # Calculate CV(numclusters) and SE_CV(numclusters) for each possible 
  #   number of clusters, up to max.levels
  # scores = -sqrt(cvm_sqerr) # codacore:::findBestCutoff.CoDaBaseLearner
  cvse = apply(
    cvm_sqerr, 2, 
    function(x) stats::sd(x, na.rm = TRUE) / sqrt(length(x) - sum(is.na(x))))
  cvm = apply(
    cvm_sqerr, 2, function(x) mean(x, na.rm = TRUE))
  
  # Find numclust_min = argmin{CV(numclusters)}
  numclust_min_index = which.min(cvm)
  numclust_min = num.levels.candidates[numclust_min_index]
  
  # Find numclust_1se = maximal numclusters s.t. CV(numclusters) <= CV(nclust_min) + CV_SE(nclust_min)
  oneserule = cvm[numclust_min_index] + cvse[numclust_min_index]
  if(!is.na(oneserule)){
    which_numclust_1se = which(cvm <= oneserule)
    numclust_1se_index = which_numclust_1se[length(which_numclust_1se)]
    numclust_1se = num.levels.candidates[numclust_1se_index]
  } else{
    numclust_1se_index = numclust_min_index
    numclust_1se = numclust_min
  }
  return(
    list(
      nclusters = num.levels.candidates,
      max.levels = max.levels, 
      models = hslrmodels,
      nclusters_min = numclust_min, 
      nclusters_1se = numclust_1se, 
      cvm_sqerr = cvm_sqerr,
      cvm = cvm, 
      nclusters_min_idx = numclust_min_index,
      cvse = cvse,
      nclusters_1se_idx = numclust_1se_index
    )
  )
}

cv.hslr1sc = function(
  x, y, max.levels = 5, nfolds = 5, classification = FALSE, approx = TRUE, 
  amini.regularization = TRUE, 
  amini.regularization.parameter = 0.01, 
  highdegree.regularization.summary = "mean",
  highdegree.regularization = FALSE,
  include.leading.eigenvector = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  hslrfit = hslr1sc(
    x = x, y = y, num.levels = max.levels, classification = classification, 
    approx = approx, 
    amini.regularization = amini.regularization, 
    amini.regularization.parameter = amini.regularization.parameter, 
    highdegree.regularization.summary = highdegree.regularization.summary,
    highdegree.regularization = highdegree.regularization,
    include.leading.eigenvector = include.leading.eigenvector)
  num.levels.candidates = 1:hslrfit$true.num.levels
  hslrmodels = hslrfit$models
  
  # split the data into nfolds folds
  shuffle = sample(1:n)
  idfold = (shuffle %% nfolds) + 1
  n_fold = table(idfold)
  
  # calculate lasso for each fold removed
  cvm_sqerr = matrix(NA, nfolds, hslrfit$true.num.levels) # squared error for each fold
  for(j in 1:nfolds){
    # Training data
    xtr = x[idfold != j, ]
    ytr = y[idfold != j]
    # Test data
    xte = x[idfold == j, ]
    yte = y[idfold == j]
    
    fit_j = hslr1sc(
      x = xtr, y = ytr, num.levels = hslrfit$true.num.levels, 
      classification = classification, approx = FALSE)
    
    # for num.clusters = 1, ..., max.levels, calculate squared error
    for(m in num.levels.candidates){
      if(length(fit_j$models) >= m){
        fit_jm = fit_j$models[[m]]
        if(!classification){
          fit_jm.coefs = coefficients(fit_jm$model)
          ypred = fit_jm.coefs[1] + fit_jm.coefs[-1] * 
            balance::balance.fromSBP(x = xte, y = fit_jm$sbp)
          cvm_sqerr[j, m] = mean(crossprod(yte - ypred))
        } else{
          # how to do for classification? use auc? accuracy? see selbal.cv #################################################
          ypred = predict.glm(
            fit_jm$model,
            newdata = data.frame(
              V1 = balance::balance.fromSBP(x = xte, y = fit_jm$sbp)),
            type = "response")
          cvm_sqerr[j, m] = sum(crossprod(yte - ypred))
        }
      } else{
        cvm_sqerr[j, m] = NA
        # warning("cv.hslr(): cvm_sqerr has NA values.")
      }
    }
  }
  # Calculate CV(numclusters) and SE_CV(numclusters) for each possible 
  #   number of clusters, up to max.levels
  # scores = -sqrt(cvm_sqerr) # codacore:::findBestCutoff.CoDaBaseLearner
  cvse = apply(
    cvm_sqerr, 2, 
    function(x) stats::sd(x, na.rm = TRUE) / sqrt(length(x) - sum(is.na(x))))
  cvm = apply(
    cvm_sqerr, 2, function(x) mean(x, na.rm = TRUE))
  
  # Find numclust_min = argmin{CV(numclusters)}
  numclust_min_index = which.min(cvm)
  numclust_min = num.levels.candidates[numclust_min_index]
  
  # Find numclust_1se = maximal numclusters s.t. CV(numclusters) <= CV(nclust_min) + CV_SE(nclust_min)
  oneserule = cvm[numclust_min_index] + cvse[numclust_min_index]
  if(!is.na(oneserule)){
    which_numclust_1se = which(cvm <= oneserule)
    numclust_1se_index = which_numclust_1se[length(which_numclust_1se)]
    numclust_1se = num.levels.candidates[numclust_1se_index]
  } else{
    numclust_1se_index = numclust_min_index
    numclust_1se = numclust_min
  }
  return(
    list(
      nclusters = num.levels.candidates,
      max.levels = max.levels, 
      models = hslrmodels,
      nclusters_min = numclust_min, 
      nclusters_1se = numclust_1se, 
      cvm_sqerr = cvm_sqerr,
      cvm = cvm, 
      nclusters_min_idx = numclust_min_index,
      cvse = cvse,
      nclusters_1se_idx = numclust_1se_index
    )
  )
}
