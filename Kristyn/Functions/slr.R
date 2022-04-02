# slr() function with the following changes:
#   1. rank1approx option added
#   2. using balances correlations instead of effect sizes to choose active subset
#   3. spectral.clustering2() instead of spectral.clustering()
slrT2 <- function(
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

cv.slr = function(
  x, y, max.clusters = 5, nfolds = 5, classification = FALSE, approx = TRUE, 
  check.significance = FALSE, alpha = 0.05
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  numclusters_candidates = 1:max.clusters
  slrmodels = list()
  for(i in numclusters_candidates){
    slrmodels[[i]] = slr(
      x = x, y = y, num.clusters = i, classification = classification, 
      approx = approx, check.significance = check.significance, alpha = alpha)
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
        approx = approx, check.significance = check.significance, alpha = alpha)
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
  cvse = apply(cvm_sqerr, 2, stats::sd)/sqrt(K)
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
  check.significance = FALSE, alpha = 0.05
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
      x = x[, vars.cur, drop = FALSE], 
      y = y, num.clusters = 2, classification = classification, 
      approx = approx)
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
  check.significance = FALSE, alpha = 0.05
){
  n = nrow(x)
  p = ncol(x)
  
  # fit slr on the original data set for each cluster size
  hslrfit = hslr(
    x = x, y = y, num.levels = max.levels, classification = classification, 
    approx = approx, check.significance = check.significance, alpha = alpha)
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
      classification = classification, 
      approx = FALSE, check.significance = check.significance, alpha = alpha)
    
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
  which_numclust_1se = which(cvm <= oneserule)
  numclust_1se_index = which_numclust_1se[length(which_numclust_1se)]
  numclust_1se = num.levels.candidates[numclust_1se_index]
  return(
    list(
      nclusters = num.levels.candidates,
      max.levels = max.levels, 
      models = hslrmodels,
      nclusters_min = numclust_min, 
      nclusters_1se = numclust_1se, 
      cvm = cvm, 
      nclusters_min_idx = numclust_min_index,
      cvse = cvse,
      nclusters_1se_idx = numclust_1se_index
    )
  )
}


