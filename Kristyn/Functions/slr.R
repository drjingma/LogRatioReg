# slr() function with the following changes:
#   1. rank1approx option added
#   2. using balances correlations instead of effect sizes to choose active subset
#   3. spectral.clustering2() instead of spectral.clustering()
slr_old <- function(
  x, y, rank1approx = FALSE, classification = FALSE, alpha = 0.05){
  p <- ncol(x)
  
  ## Compute pairwise correlation 
  rhoMat <- slrmatrix(x = x, y = y)
  out <- list()
  out$kernel <- rhoMat
  
  ## Split into active/inactive sets
  if(rank1approx){
    rhoMat.svd <- svd(rhoMat)
    rhoMat_approx_1 <- tcrossprod(
      rhoMat.svd$u[, 1], rhoMat.svd$v[, 1]) * rhoMat.svd$d[1]
    index <- which(spectral.clustering2(rhoMat_approx_1, reindex = TRUE) == 1)
  } else{
    index <- which(spectral.clustering2(rhoMat, reindex = TRUE) == 1)
  }
  
  ## Find which set has the more predictive balance
  # if(length(index) == 1 | length(index) == p - 1){ # one valid cluster
  #   warning("slr(): Only one covariate was subsetted. No log-ratios can be made.")
  #   # fit the subset with length p - 1
  #   if(length(index) == p - 1){
  #     subset <- spectral.clustering2(rhoMat[index, index], reindex = TRUE)
  #   } else{
  #     subset <- spectral.clustering2(rhoMat[-index, -index], reindex = TRUE)
  #   }
  #   sbp.est <- matrix(0, nrow = p, ncol = 1)
  #   rownames(sbp.est) <- colnames(x)
  #   sbp.est[match(names(subset),rownames(sbp.est)), 1] <- subset
  #   est.balance <- balance::balance.fromSBP(x = x, y = sbp.est)
  #   if(!classification){
  #     refit2 <- lm(y~est.balance)
  #     pval = summary(refit2)$coefficients["est.balance", "Pr(>|t|)"]
  #   } else{
  #     refit2 = stats::glm(
  #       y~est.balance, family = binomial(link = "logit"))
  #     pval = summary(refit2)$coefficients["est.balance", "Pr(>|z|)"]
  #   }
  #   if(pval < alpha){
  #     out$index = subset
  #     refit = refit2
  #   } else{
  #     # fit only the intercept
  #     out$index = NA
  #     if(!classification){
  #       refit <- lm(y~1)
  #     } else{
  #       refit = stats::glm(y~1, family = binomial(link = "logit"))
  #     }
  #   }
  # } else{
  ## Perform spectral clustering to get the numerator/denominator groups
  ##    for the active/inactive sets
  subset1 <- spectral.clustering2(rhoMat[index, index], reindex = TRUE)
  subset2 <- spectral.clustering2(rhoMat[-index, -index], reindex = TRUE) 
  
  # Since we don't know which subset contains the active variables, 
  ## we first fit a linear model with balances obtained from both subsets. 
  sbp.est <- matrix(0, ncol = 2, nrow = p)
  rownames(sbp.est) <- colnames(x)
  sbp.est[match(names(subset1),rownames(sbp.est)), 1] <- subset1
  sbp.est[match(names(subset2),rownames(sbp.est)), 2] <- subset2
  est.balance <- balance::balance.fromSBP(x = x, y = sbp.est)
  cors = stats::cor(est.balance, y)
  
  # The correct subset should have larger coefficient. 
  ## We refit the linear model on the balance from the correct subset. 
  out$cors = cors
  
  if ( abs(cors[1, ]) > abs(cors[2, ]) ){
    # pick subset1
    out$index <- subset1
    if(!classification){
      refit <- lm(y~est.balance[, 1])
    } else{
      refit = stats::glm(
        y~est.balance[, 1], family = binomial(link = "logit"))
    }
  } else {
    # pick subset2
    out$index <- subset2
    if(!classification){
      refit <- lm(y~est.balance[, 2])
    } else{
      refit = stats::glm(
        y~est.balance[, 2], family = binomial(link = "logit"))
    }
  }
  # }
  
  out$model <- refit
  out
}

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
spectral.clustering2 <- function(W, n_eig = 2, reindex = FALSE) {
  L = graph.laplacian(W)          # 2. compute graph laplacian
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  # we will use k-means to cluster the data
  # using the leading eigenvalues in absolute values
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  if(nrow(W) == 2){
    obj = list(cluster= c(1, 2))
  } else{
    obj <- kmeans(
      ei$vectors[, 2:n_eig, drop = FALSE], centers = n_eig, nstart = 100)
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

slr <- function(
  x, y, num.clusters = 2, classification = FALSE, approx = TRUE, 
  check.significance = FALSE, alpha = 0.05
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
    clusters1 <- spectral.clustering2(rhoMat_approx, n_eig = num.clusters)
  } else{
    clusters1 <- spectral.clustering2(rhoMat, n_eig = num.clusters)
  }
  cluster.lengths = table(clusters1)
  valid.clusters = names(which(cluster.lengths >= 2)) # sets that can make log-ratios
  valid.cluster.lengths = cluster.lengths[valid.clusters]
  out$num.clusters = num.clusters
  out$num.valid.clusters = length(valid.clusters)
  
  ## Find which set has the more predictive balance
  cors = rep(NA, length(valid.clusters))
  Rsqs = rep(NA, length(valid.clusters))
  sbp.ests = matrix(0, nrow = p, ncol = length(valid.clusters))
  bal.ests = matrix(NA, nrow = n, ncol = length(valid.clusters))
  rownames(sbp.ests) <- colnames(x)
  for(i in 1:length(valid.clusters)){ # what if there is just one valid cluster? ################
    cluster.label = as.numeric(valid.clusters[i])
    index = which(clusters1 == cluster.label)
    ## Perform spectral clustering to get the numerator/denominator groups
    subset = spectral.clustering2(rhoMat[index, index], reindex = TRUE)
    ## calculate balance and its correlation with y
    sbp.ests[match(names(subset), rownames(sbp.ests)), i] = subset
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
  selected.cluster = which.max(Rsqs)
  out$index = sbp.ests[, selected.cluster]
  if(!classification){
    refit <- lm(y~bal.ests[, selected.cluster])
  } else{
    refit = stats::glm(
      y~bal.ests[, selected.cluster], family = binomial(link = "logit"))
  }
  out$model <- refit
  
  # return the full SBP vector (with entries for all p variables)
  full.sbp.est = sbp.ests[, selected.cluster, drop = FALSE]
  rownames(full.sbp.est) = colnames(x)
  out$sbp = full.sbp.est
  
  # return the slr object
  return(out)
}

slrmult = function(
  x, y, max.clusters = 5, classification = FALSE, approx = TRUE, 
  check.significance = FALSE, alpha = 0.05
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  models = list()
  cors = matrix(NA, nrow = max.clusters, ncol = max.clusters)
  Rsqs = matrix(NA, nrow = max.clusters, ncol = max.clusters)
  numclusters_candidates = 1:max.clusters
  for(i in numclusters_candidates){
    model_i = slr(
      x = x, y = y, num.clusters = i, classification = classification, 
      approx = approx, check.significance = check.significance, alpha = alpha)
    models[[i]] = model_i
    cors[i, ] = c(model_i$cors, rep(NA, max.clusters - length(model_i$cors)))
    Rsqs[i, ] = c(model_i$Rsqs, rep(NA, max.clusters - length(model_i$Rsqs)))
  }
  
  # return the slr object
  return(list(
    models = models, 
    cors = cors, 
    Rsqs = Rsqs, 
    max.cors = rowMaxs(cors, na.rm = TRUE),
    max.Rsqs = rowMaxs(Rsqs, na.rm = TRUE)
  ))
}

cv.slr = function(
  x, y, max.clusters = 5, nfolds = 5, classification = FALSE, approx = TRUE, 
  check.significance = FALSE, alpha = 0.05
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  numclusters_candidates = 1:max.clusters
  slrmultmodels = slrmult(
    x = x, y = y, max.clusters = max.clusters, classification = classification, 
    approx = approx, check.significance = check.significance, alpha = alpha)
  
  # split the data into nfolds folds
  shuffle = sample(1:n)
  idfold = (shuffle %% nfolds) + 1
  n_fold = table(idfold)
  
  # calculate lasso for each fold removed
  cvm_sqerror = matrix(NA, nfolds, max.clusters) # squared error for each fold
  cvse_errorsd =  matrix(NA, nfolds, max.clusters) # sd of error for each fold
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
        approx = approx, check.significance = check.significance, alpha = alpha)
      if(!classification){
        fit_jm.coefs = coefficients(fit_jm$model)
        ypred = fit_jm.coefs[1] + fit_jm.coefs[-1] * 
          balance::balance.fromSBP(x = xte, y = fit_jm$sbp)
        cvm_sqerror[j, m] = sum(crossprod(yte - ypred))
        cvse_errorsd[j, m] = cvm_sqerror[j, m] / n_fold[j]
      } else{
        # how to do for classification? use auc? accuracy? #################################################
        ypred = predict.glm(
          fit_jm$model,
          newdata = data.frame(balance::balance.fromSBP(x = xte, y = fit_jm$sbp)),
          type = "response")
        cvm_sqerror[j, m] = sum(crossprod(yte - ypred))
        cvse_errorsd[j, m] = cvm_sqerror[j, m] / n_fold[j]
      }
    }
  }
  # Calculate CV(numclusters) and SE_CV(numclusters) for each possible 
  #   number of clusters, up to max.clusters
  cvse = sqrt(apply(cvse_errorsd, 2, var)) / sqrt(nfolds)
  cvm = colMeans(cvm_sqerror)
  
  # Find numclust_min = argmin{CV(numclusters)}
  numclust_min_index = which.min(cvm)
  numclust_min = numclusters_candidates[numclust_min_index]
  
  # Find numclust_1se = maximal numclusters s.t. CV(numclusters) <= CV(nclust_min) + CV_SE(nclust_min)
  oneserule = cvm[numclust_min_index] + cvse[numclust_min_index]
  numclust_1se_index = which(cvm <= oneserule)[1]
  numclust_1se = numclusters_candidates[numclust_1se_index]
  return(
    list(
      nclusters = numclusters_candidates,
      max.clusters = max.clusters, 
      slrmult = slrmultmodels, 
      models = slrmultmodels$models,
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
  x, y, max.levels = 1, classification = FALSE, approx = TRUE
){
  n = nrow(x)
  p = ncol(x)
  
  fits = list()
  Rsqs = matrix(NA, nrow = max.levels, ncol = 2)
  cors = matrix(NA, nrow = max.levels, ncol = 2)
  vars.cur = colnames(x)
  num.levels = max.levels
  for(i in 1:max.levels){
    fit_i = slr(
      x = x[, vars.cur, drop = FALSE], 
      y = y, num.clusters = 2, classification = classification, 
      approx = approx)
    # recalculate full SBP
    activevars = vars.cur
    full.sbp.est = matrix(0, nrow = p, ncol = 1)
    rownames(full.sbp.est) = colnames(x)
    full.sbp.est[activevars, ] = c(fit_i$sbp)
    fit_i$sbp = full.sbp.est
    # save the fit
    fits[[i]] = fit_i
    vars.cur = names(fit_i$index)[fit_i$index != 0]
    cors[i, ] = fit_i$cors
    Rsqs[i, ] = fit_i$Rsqs
    if(length(vars.cur) <= 2){
      num.levels = i
      # warning(paste0("asked for max.levels = ", max.levels, " splits, but could only make up to ", i, " splits."))
      break
    }
  }
  
  # return the slr object
  return(list(
    models = fits, 
    cors = cors, 
    Rsqs = Rsqs, 
    max.cors = rowMaxs(cors, na.rm = TRUE),
    max.Rsqs = rowMaxs(Rsqs, na.rm = TRUE), 
    num.levels = num.levels,
    max.levels = max.levels
  ))
}

# cv.hslr = function(
#   x, y, max.levels = 5, classification = FALSE, approx = TRUE, nfolds = 5
# ){
#   n = nrow(x)
#   p = ncol(x)
#   
#   # fit slr on the original data set for each cluster size
#   fits = list()
#   numclusters_candidates = 1:max.clusters
#   for(i in numclusters_candidates){
#     fits[[i]] = slr(
#       x = x, y = y, num.clusters = i, classification = classification, 
#       approx = approx)
#   }
#   
#   # split the data into nfolds folds
#   shuffle = sample(1:n)
#   idfold = (shuffle %% nfolds) + 1
#   n_fold = table(idfold)
#   
#   # calculate lasso for each fold removed
#   cvm_sqerror = matrix(NA, nfolds, max.clusters) # squared error for each fold
#   cvse_errorsd =  matrix(NA, nfolds, max.clusters) # sd of error for each fold
#   for(j in 1:nfolds){
#     # Training data
#     xtr = x[idfold != j, ]
#     ytr = y[idfold != j]
#     # Test data
#     xte = x[idfold == j, ]
#     yte = y[idfold == j]
#     
#     # for num.clusters = 1, ..., max.clusters, calculate squared error
#     for(m in numclusters_candidates){
#       fit_jm = slr(
#         x = xtr, y = ytr, num.clusters = m, classification = classification, 
#         approx = approx)
#       if(!classification){
#         fit_jm.coefs = coefficients(fit_jm$model)
#         ypred = fit_jm.coefs[1] + fit_jm.coefs[-1] * 
#           balance::balance.fromSBP(x = xte, y = fit_jm$sbp)
#         cvm_sqerror[j, m] = sum(crossprod(yte - ypred))
#         cvse_errorsd[j, m] = cvm_sqerror[j, m] / n_fold[j]
#       } else{
#         # how to do for classification? use auc? #################################################
#         ypred = predict.glm(
#           fit_jm$model,
#           newdata = data.frame(balance::balance.fromSBP(x = xte, y = fit_jm$sbp)),
#           type = "response")
#         cvm_sqerror[j, m] = sum(crossprod(yte - ypred))
#         cvse_errorsd[j, m] = cvm_sqerror[j, m] / n_fold[j]
#       }
#     }
#   }
#   # Calculate CV(numclusters) and SE_CV(numclusters) for each possible 
#   #   number of clusters, up to max.clusters
#   cvse = sqrt(apply(cvse_errorsd, 2, var)) / sqrt(nfolds)
#   cvm = colMeans(cvm_sqerror)
#   
#   # Find numclust_min = argmin{CV(numclusters)}
#   numclust_min_index = which.min(cvm)
#   numclust_min = numclusters_candidates[numclust_min_index]
#   
#   # Find numclust_1se = maximal numclusters s.t. CV(numclusters) <= CV(nclust_min) + CV_SE(nclust_min)
#   oneserule = cvm[numclust_min_index] + cvse[numclust_min_index]
#   numclust_1se_index = which(cvm <= oneserule)
#   numclust_1se = numclusters_candidates[min(numclust_1se_index)]
#   return(
#     list(
#       nclusters = numclusters_candidates,
#       max.clusters = max.clusters, 
#       models = fits, 
#       nclusters_min = numclust_min, 
#       nclusters_1se = numclust_1se, 
#       cvm = cvm, 
#       nclusters_min_idx = numclust_min_index,
#       cvse = cvse,
#       nclusters_1se_idx = numclust_1se_index
#     )
#   )
# }