
# old name: fitILReta
fitBMLassoThresh = function(
  y, X, 
  W, # similarity matrix
  sbp,
  hsc_method = "kmeans", # "shimalik", "kmeans"
  force_levelMax = FALSE, 
  stopping_rule = NULL,
  #   NULL means none
  #   "natural" means force_levelMax = FALSE, 
  #   "TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm"
  lambda = NULL, nlam = 20, 
  eta = NULL, neta = 20,
  nfolds = 5,
  intercept = TRUE, standardize = TRUE
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n){
    stop("fitBMLassoThresh(): dimensions of y and X don't match!!")
  }
  if(is.null(colnames(X))){
    colnames(X) = paste("V", 1:p, sep = "")
  }
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  if(nrow(W) != p | ncol(W) != p){
    stop("fitBMLassoThresh(): W isn't pxp matrix!!")
  }
  if(is.null(rownames(sbp))){
    warning("fitBMLassoThresh(): sbp matrix doesn't have named rows!!")
  }
  
  # get lambda (if not provided)
  if(is.null(lambda) | is.null(nlam)){
    glmnet.fit = glmnet(
      x = getIlrX(X = X, sbp = sbp), y = y, 
      lambda = lambda, nlambda = nlam,
      intercept = intercept, standardize = standardize)
    lambda = glmnet.fit$lambda
    # check lambda length -- if it doesn't have length nlam, rerun
    if(nlam != length(lambda)){
      log_lambda <- log(lambda)
      lambda <- exp(seq(max(log_lambda), min(log_lambda),length.out = nlam))
    }
  }
  
  # get eta (if not provided)
  if(is.null(eta) | is.null(neta)){
    if(is.null(neta)) neta = 5
    eta = seq(0, 1, length.out = neta + 1)[2:(neta + 1)]
  } else{
    if(neta != length(eta)) stop("fitBMLassoThresh(): neta != length(eta)")
  }
  
  # thresholding with eta: Iterate solution paths along eta
  meets_threshold <- clust_thresh <- sbp_thresh <- theta0 <- theta <- list()
  num_covariates = rep(NA, neta)
  for(i in 1:neta){
    # thresholding
    # using a similarity matrix (close to 1 = highly correlated with y)
    # if there is an element that is greater than eta
    meets_threshold_i = apply(W, 1, function(row) !all(row < eta[i]))
    num_covariates[i] = sum(meets_threshold_i)
    if(sum(meets_threshold_i) <= 2){ # cannot cluster
      theta[[i]] = rep(NA, sum(meets_threshold_i))
      theta0[[i]] = NA
      meets_threshold[[i]] = meets_threshold_i
      sbp_thresh[[i]] = rep(NA, sum(meets_threshold_i))
    } else{
      W_thresh = W[meets_threshold_i, meets_threshold_i, drop = FALSE]
      X_thresh = X[, meets_threshold_i, drop = FALSE]
      # model fitting
      hsclust_thresh = HSClust(
        W = W_thresh, force_levelMax = force_levelMax, 
        stopping_rule = stopping_rule, method = hsc_method
      )
      SBP_thresh = sbp.fromHSClust(
        levels_matrix = hsclust_thresh$allLevels, 
        row_names = rownames(sbp)[meets_threshold_i])
      modelfit = cvBMLasso(
        y = y, X = X_thresh, sbp = SBP_thresh, lambda = lambda, nlam = nlam, 
        nfolds = nfolds, intercept = intercept, standardize = standardize)
      # save the model
      theta[[i]] = modelfit$theta
      theta0[[i]] = modelfit$theta0
      meets_threshold[[i]] = meets_threshold_i
      clust_thresh[[i]] = hsclust_thresh
      sbp_thresh[[i]] = SBP_thresh
    }
  }
  
  return(list(
    theta0 = theta0,
    theta = theta,
    meets_threshold = meets_threshold,
    sbp_thresh = sbp_thresh,
    clust_thresh = clust_thresh,
    num_covariates = num_covariates, 
    lambda = lambda,
    eta = eta, 
    W = W,
    sbp = sbp
  ))
}

# old name: cvILReta
cvBMLassoThresh <- function(
  y, X,
  W, # similarity matrix
  sbp,
  hsc_method = "kmeans", # "shimalik", "kmeans"
  force_levelMax = TRUE, 
  stopping_rule = NULL, 
  #   NULL means none
  #   "natural" means force_levelMax = FALSE, 
  #   "TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm"
  lambda = NULL, nlam = 20, 
  eta = NULL, neta = 5,
  nfolds = 10, foldid = NULL, 
  intercept = TRUE, standardize = TRUE,
  seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n){
    stop("cvBMLassoThresh(): dimensions of y and X don't match")
  }
  if(is.null(colnames(X))){
    colnames(X) = paste("V", 1:p, sep = "")
  }
  # check if lambda is given, assign nlam accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  if(nrow(W) != p | ncol(W) != p){
    stop("cvBMLassoThresh(): W isn't pxp matrix")
  }
  if(is.null(rownames(sbp))){
    warning("cvBMLassoThresh(): sbp matrix doesn't have named rows")
  }
  
  # get lambda (if not provided)
  if(is.null(lambda) | is.null(nlam)){
    glmnet.fit = glmnet(
      x = getIlrX(X = X, sbp = sbp), y = y, 
      lambda = lambda, nlambda = nlam,
      intercept = intercept, standardize = standardize)
    lambda = glmnet.fit$lambda
    # check lambda length -- if it doesn't have length nlam, rerun
    if(nlam != length(lambda)){
      log_lambda <- log(lambda)
      lambda <- exp(seq(max(log_lambda), min(log_lambda),length.out = nlam))
    }
  }
  
  # get eta (if not provided)
  if(is.null(eta) | is.null(neta)){
    if(is.null(neta)) neta = 5
    eta = seq(0, max(W), length.out = neta + 1)[2:(neta + 1)]
  } else{
    if(neta != length(eta)){
      stop("cvBMLassoThresh(): neta != length(eta)")
    }
  }
  
  # fit the models
  fitObj = fitBMLassoThresh(
    y = y, X = X, W = W, hsc_method = hsc_method,
    force_levelMax = force_levelMax, stopping_rule = stopping_rule, sbp = sbp, 
    lambda = lambda, nlam = nlam, eta = eta, neta = neta, nfolds = nfolds,
    intercept = intercept, standardize = standardize)
  
  # nlam <- length(fitObj$lambda)
  # nalpha <- length(fitObj$alpha)
  errs <- array(NA, dim=c(nlam, neta, nfolds))
  
  # define error function
  errfun <- function(est, truth) colMeans((est - truth)^2, na.rm = TRUE)
  
  # make folds, if necessary
  if(is.null(foldid)){
    sizes_equal <- round(n / nfolds)
    sizes <- rep(sizes_equal, nfolds)
    sizes[nfolds] <- sizes[nfolds] + n - sizes_equal * nfolds
    endpoints <- c(0, cumsum(sizes))
    perm <- sample(n)
    folds <- list()
    for (i in 1:nfolds) folds[[i]] <- perm[
      seq(endpoints[i] + 1, endpoints[i + 1])]
  } else{
    folds = list()
    for(i in 1:nfolds) folds[[i]] = which(foldid == i)
  }
  
  # Fit based on folds and compute error metric
  for (i in 1:nfolds) {
    # fit model on all but the ith fold
    fit_cv <- fitBMLassoThresh(
      y = y[-folds[[i]]], X = X[-folds[[i]], , drop = FALSE], W = W, sbp = sbp,
      hsc_method = hsc_method,
      force_levelMax = force_levelMax, stopping_rule = stopping_rule,
      lambda = lambda, nlam = nlam, eta = eta, neta = neta, nfolds = nfolds,
      intercept = intercept, standardize = standardize)
    pred_te <- lapply(1:neta, function(k) {
      if(all(is.na(fit_cv$sbp_thresh[[k]]))){
        matrix(NA, nrow = length(folds[[i]]), ncol = neta)
      } else{
        if (intercept) {
          # Xb[folds[[i]], fit_cv$meets_threshold[[k]], drop = FALSE] %*% 
          getIlrX(
            X = X[folds[[i]], fit_cv$meets_threshold[[k]], drop = FALSE], 
            sbp = fit_cv$sbp_thresh[[k]]) %*%
            fit_cv$theta[[k]] +
            rep(fit_cv$theta0[[k]], each = length(folds[[i]]))
        } else {
          # X[folds[[i]], ] %*% fit_cv$beta[[k]]
          # Xb[folds[[i]], fit_cv$meets_threshold[[k]]]
          getIlrX(
            X = X[folds[[i]], fit_cv$meets_threshold[[k]], drop = FALSE], 
            sbp = fit_cv$sbp_thresh[[k]]) %*% 
            fit_cv$theta[[k]]
        }
      }
    })
    for (k in 1:neta) errs[, k, i] <- errfun(pred_te[[k]], y[folds[[i]]])
    cat("##########################\n")
    cat(sprintf("Finished model fits for fold[%s].\n", i))
    cat("##########################\n")
  }
  m <- apply(errs, c(1, 2), function(x) mean(x, na.rm = TRUE)) # rows = lambda, cols = eta
  # se <- apply(errs, c(1, 2), function(x) stats::sd(x, na.rm = TRUE)) / sqrt(nfolds)
  ibest <- which(m == min(m, na.rm = TRUE), arr.ind = TRUE)[1, , drop = FALSE]
  
  return(list(
    theta0 = fitObj$theta0,
    theta = fitObj$theta,
    meets_threshold = fitObj$meets_threshold,
    clust_thresh = fitObj$clust_thresh,
    sbp_thresh = fitObj$sbp_thresh,
    num_covariates = fitObj$num_covariates,
    lambda = fitObj$lambda,
    eta = eta, 
    fits = fitObj,
    sbp = sbp,
    cvm = m,
    min.idx = ibest
  ))
}












# old name: fitILReta
fitBMThresh = function(
  y, X, 
  W, # similarity matrix
  sbp,
  hsc_method = "kmeans", # "shimalik", "kmeans"
  multiple_balances = TRUE,
  force_levelMax = FALSE, 
  stopping_rule = NULL,
  #   NULL means none
  #   "natural" means force_levelMax = FALSE, 
  #   "TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm"
  eta = NULL, neta = 20,
  intercept = TRUE, standardize = TRUE
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n){
    stop("fitBMThresh(): dimensions of y and X don't match!!")
  }
  if(is.null(colnames(X))){
    colnames(X) = paste("V", 1:p, sep = "")
  }
  if(nrow(W) != p | ncol(W) != p){
    stop("fitBMThresh(): W isn't pxp matrix!!")
  }
  if(is.null(rownames(sbp))){
    warning("fitBMThresh(): sbp matrix doesn't have named rows!!")
  }
  
  # get eta (if not provided)
  if(is.null(eta) | is.null(neta)){
    if(is.null(neta)) neta = 5
    eta = seq(0, 1, length.out = neta + 1)[2:(neta + 1)]
  } else{
    if(neta != length(eta)) stop("fitBMThresh(): neta != length(eta)")
  }
  
  # thresholding with eta: Iterate solution paths along eta
  meets_threshold <- clust_thresh <- sbp_thresh <- theta0 <- theta <- list()
  num_covariates = rep(NA, neta)
  for(i in 1:neta){
    # thresholding
    # using a similarity matrix (close to 1 = highly correlated with y)
    # if there is an element that is greater than eta
    meets_threshold_i = apply(W, 1, function(row) !all(row < eta[i]))
    num_covariates[i] = sum(meets_threshold_i)
    if(sum(meets_threshold_i) <= 2){ # cannot cluster
      theta[[i]] = rep(NA, sum(meets_threshold_i))
      theta0[[i]] = NA
      meets_threshold[[i]] = meets_threshold_i
      sbp_thresh[[i]] = rep(NA, sum(meets_threshold_i))
    } else{
      W_thresh = W[meets_threshold_i, meets_threshold_i, drop = FALSE]
      X_thresh = X[, meets_threshold_i, drop = FALSE]
      # model fitting
      hsclust_thresh = HSClust(
        W = W_thresh, force_levelMax = force_levelMax, 
        stopping_rule = stopping_rule, method = hsc_method
      )
      SBP_thresh = sbp.fromHSClust(
        levels_matrix = hsclust_thresh$allLevels, 
        row_names = rownames(sbp)[meets_threshold_i])
      # modelfit = cvBMLasso(
      #   y = y, X = X_thresh, sbp = SBP_thresh, 
      #   nfolds = nfolds, intercept = intercept, standardize = standardize)
      if(multiple_balances){
        ilrX_thresh = getIlrX(X = X_thresh, sbp = SBP_thresh)
      } else {
        ilrX_thresh = getIlrX(X = X_thresh, sbp = SBP_thresh[, 1, drop = FALSE])
      }
      if(intercept){
        modelfit = lm(y ~ ilrX_thresh)
        coeffs = coefficients(modelfit)
        # save the model
        theta[[i]] = coeffs[-1]
        theta0[[i]] = coeffs[1]
      } else{
        modelfit = lm(y ~ -1 + ilrX_thresh)
        coeffs = coefficients(modelfit)
        # save the model
        theta[[i]] = coefs
        theta0[[i]] = NA
      }
      # save stuff
      meets_threshold[[i]] = meets_threshold_i
      clust_thresh[[i]] = hsclust_thresh
      sbp_thresh[[i]] = SBP_thresh
    }
  }
  
  return(list(
    theta0 = theta0,
    theta = theta,
    meets_threshold = meets_threshold,
    sbp_thresh = sbp_thresh,
    clust_thresh = clust_thresh,
    num_covariates = num_covariates, 
    eta = eta, 
    W = W,
    sbp = sbp
  ))
}

# old name: cvILReta
cvBMThresh <- function(
  y, X,
  W, # similarity matrix
  sbp,
  hsc_method = "kmeans", # "shimalik", "kmeans"
  multiple_balances = TRUE, 
  force_levelMax = TRUE, 
  stopping_rule = NULL, 
  #   NULL means none
  #   "natural" means force_levelMax = FALSE, 
  #   "TooManyCells", "newmangirmanmodularity", "ngmod", "tmc", "ngm"
  eta = NULL, neta = 5,
  nfolds = 10, foldid = NULL, 
  intercept = TRUE, standardize = TRUE,
  seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n){
    stop("cvBMThresh(): dimensions of y and X don't match")
  }
  if(is.null(colnames(X))){
    colnames(X) = paste("V", 1:p, sep = "")
  }
  if(nrow(W) != p | ncol(W) != p){
    stop("cvBMThresh(): W isn't pxp matrix")
  }
  if(is.null(rownames(sbp))){
    warning("cvBMThresh(): sbp matrix doesn't have named rows")
  }
  
  # get eta (if not provided)
  if(is.null(eta) | is.null(neta)){
    if(is.null(neta)) neta = 5
    eta = seq(0, max(W), length.out = neta + 1)[2:(neta + 1)]
  } else{
    if(neta != length(eta)){
      stop("cvBMThresh(): neta != length(eta)")
    }
  }
  
  # fit the models
  fitObj = fitBMThresh(
    y = y, X = X, W = W, hsc_method = hsc_method, 
    multiple_balances = multiple_balances,
    force_levelMax = force_levelMax, stopping_rule = stopping_rule, sbp = sbp, 
    eta = eta, neta = neta, intercept = intercept, standardize = standardize)
  
  # make folds, if necessary
  if(is.null(foldid)){
    shuffle = sample(1:n)
    idfold = (shuffle %% nfolds) + 1
  } else {
    idfold = foldid
  }
  n_fold = as.vector(table(idfold))
  
  
  # Fit based on folds and compute error metric
  # calculate squared error for each fold, needed for CV(eta) calculation
  cvm_sqerror = matrix(rep(NA, nfolds * neta), nfolds, neta)
  for (i in 1:nfolds) {
    Xtrain_i = X[idfold != i, , drop = FALSE]
    Ytrain_i = y[idfold != i]
    Xtest_i = X[idfold == i, , drop = FALSE]
    Ytest_i = y[idfold == i]
    # fit model on all but the ith fold
    fit_i <- fitBMThresh(
      y = Ytrain_i, X = Xtrain_i, W = W, sbp = sbp,
      hsc_method = hsc_method, force_levelMax = force_levelMax, 
      stopping_rule = stopping_rule, eta = eta, neta = neta,
      intercept = intercept, standardize = standardize)
    # calculate predicted y on test set for each eta's model
    for(k in 1:neta){
      if(all(is.na(fit_i$sbp_thresh[[k]]))){
        ypred_test_i = rep(NA, n_fold[i])
      } else {
        ypred_test_i = getIlrX(
          X = Xtest_i[, fit_i$meets_threshold[[k]], drop = FALSE], 
          sbp = fit_i$sbp_thresh[[k]]) %*%
          fit_i$theta[[k]]
        if(intercept){
          ypred_test_i = ypred_test_i + rep(fit_i$theta0[[k]], each = n_fold[i])
        }
      }
      # prediction error on test set for ith fold and kth eta's model
      cvm_sqerror[i, k] = sum(crossprod(ypred_test_i, Ytest_i))
    }
    
    # pred_te <- lapply(1:neta, function(k) {
    #   if(all(is.na(fit_i$sbp_thresh[[k]]))){
    #     rep(NA, n_fold[i])
    #   } else{
    #     if (intercept) {
    #       getIlrX(
    #         X = Xtest_i[, fit_i$meets_threshold[[k]], drop = FALSE], 
    #         sbp = fit_i$sbp_thresh[[k]]) %*%
    #         fit_i$theta[[k]] +
    #         rep(fit_i$theta0[[k]], each = n_fold[i])
    #     } else {
    #       getIlrX(
    #         X = Xtest_i[, fit_i$meets_threshold[[k]], drop = FALSE], 
    #         sbp = fit_i$sbp_thresh[[k]]) %*% 
    #         fit_i$theta[[k]]
    #     }
    #   }
    # })
    # for (k in 1:neta) cvm_sqerror[k, i] <- sum(crossprod(pred_te[[k]], Ytest_i))
    #
    cat("##########################\n")
    cat(sprintf("Finished model fits for fold[%s].\n", i))
    cat("##########################\n")
  }
  # fields::image.plot(cvm_sqerror)
  # m <- apply(errs, 1, function(x) mean(x, na.rm = TRUE)) # rows = eta
  # ibest <- which.min(m)
  cvm = colMeans(cvm_sqerror)
  min.idx = which.min(cvm)
  return(list(
    theta0 = fitObj$theta0,
    theta = fitObj$theta,
    meets_threshold = fitObj$meets_threshold,
    clust_thresh = fitObj$clust_thresh,
    sbp_thresh = fitObj$sbp_thresh,
    num_covariates = fitObj$num_covariates,
    eta = eta, 
    fits = fitObj,
    sbp = sbp,
    cvm = cvm,
    min.idx = min.idx
  ))
}



