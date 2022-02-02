fitILReta = function(
  y, X, 
  W, # normalized similarity/distance matrix (all values between 0 & 1)
  sbp,
  clustering_method = "hsc", # "hs", "hsc"
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
  if(length(y) != n) stop("fitILReta : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  # check if lambda is given, assign nlam accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  if(nrow(W) != p | ncol(W) != p) stop("fitILReta: W isn't pxp matrix")
  if(is.null(rownames(sbp))) warning("sbp matrix doesn't have named rows")
  
  # get lambda (if not provided)
  if(is.null(lambda) | is.null(nlam)){
    glmnet.fit = glmnet(
      x = computeBalances(X, sbp = sbp), y = y, 
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
    if(neta != length(eta)) stop("neta != length(eta)")
  }
  
  # thresholding with eta: Iterate solution paths along eta
  meets_threshold <- clust_thresh <- sbp_thresh <- theta0 <- theta <- list()
  num_covariates = rep(NA, neta)
  for(i in 1:neta){
    # thresholding
    if(clustering_method == "hsc"){ 
      # using a similarity matrix (close to 1 = highly correlated with y)
      # if there is an element that is greater than eta
      meets_threshold_i = apply(W, 1, function(row) !all(row < eta[i]))
    } else if(clustering_method == "hc"){ 
      # using a distance matrix (close to 1 = not highly correlated with y)
      # if there is an element that is less than eta
      meets_threshold_i = apply(W, 1, function(row) !all(row > eta[i]))
    }
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
      if(clustering_method == "hsc"){
        hsclust_thresh = HSClust(
          W = W_thresh, force_levelMax = force_levelMax, 
          stopping_rule = stopping_rule, method = hsc_method
        )
        SBP_thresh = sbp.fromHSClust(
          levels_matrix = hsclust_thresh$allLevels, 
          row_names = rownames(sbp)[meets_threshold_i])
        modelfit = cvILR(
          y = y, X = X_thresh, sbp = SBP_thresh, lambda = lambda, nlam = nlam, 
          nfolds = nfolds, intercept = intercept, standardize = standardize)
        # save the model
        theta[[i]] = modelfit$bet
        theta0[[i]] = modelfit$int
        meets_threshold[[i]] = meets_threshold_i
        clust_thresh[[i]] = hsclust_thresh
        sbp_thresh[[i]] = SBP_thresh
      } else if(clustering_method == "hc"){
        hclust_thresh = hclust(as.dist(W_thresh),method = linkage)
        SBP_thresh = sbp.fromHclust(hclust_thresh)
        modelfit = cvILR(
          y = y, X = X_thresh, sbp = SBP_thresh, lambda = lambda, nlam = nlam, 
          nfolds = nfolds, intercept = intercept, standardize = standardize)
        # save the model
        theta[[i]] = modelfit$bet
        theta0[[i]] = modelfit$int
        meets_threshold[[i]] = meets_threshold_i
        clust_thresh[[i]] = hclust_thresh
        sbp_thresh[[i]] = SBP_thresh
      }
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
cvILReta <- function(
  y, X,
  W, # normalized similarity/distance matrix (all values between 0 & 1)
  sbp,
  clustering_method = "hsc", # "hs", "hsc"
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
  if(length(y) != n) stop("cvILReta : dimensions of y and X don't match")
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  # check if lambda is given, assign nlam accordingly
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  if(nrow(W) != p | ncol(W) != p) stop("cvILReta: W isn't pxp matrix")
  if(is.null(rownames(sbp))) warning("sbp matrix doesn't have named rows")
  
  # get lambda (if not provided)
  if(is.null(lambda) | is.null(nlam)){
    glmnet.fit = glmnet(
      x = computeBalances(X, sbp = sbp), y = y, 
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
    if(neta != length(eta)) stop("neta != length(eta)")
  }
  
  # fit the models
  fitObj = fitILReta(
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
    fit_cv <- fitILReta(
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
          computeBalances(
            X[folds[[i]], fit_cv$meets_threshold[[k]], drop = FALSE], 
            sbp = fit_cv$sbp_thresh[[k]]) %*%
            fit_cv$theta[[k]] +
            rep(fit_cv$theta0[[k]], each = length(folds[[i]]))
        } else {
          # X[folds[[i]], ] %*% fit_cv$beta[[k]]
          # Xb[folds[[i]], fit_cv$meets_threshold[[k]]]
          computeBalances(
            X[folds[[i]], fit_cv$meets_threshold[[k]], drop = FALSE], 
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

