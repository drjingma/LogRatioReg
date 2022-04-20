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

