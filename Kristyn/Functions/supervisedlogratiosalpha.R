# fitSLRalpha_old <- function(
#   y, X, A = NULL, U = NULL, btree = NULL, sbp = NULL, linkage = "complete", 
#   rho.type = "squared",
#   Q = NULL, intercept = TRUE, lambda = NULL, 
#   alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
#   rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5, maxite = 1e6, scaling = FALSE
# ){
#   # scaling doesn't work rn.
#   # browser()
#   
#   if(is.null(btree)) btree = getSupervisedTree(y, X, linkage, rho.type)
#   if(is.null(U)) U = getU(btree = btree, sbp = sbp)
#   # if(is.null(Xb)) Xb = computeBalances(X, btree = btree, sbp = sbp, U = U)
#   
#   if(is.null(A)){
#     A = Matrix(U, sparse = TRUE)
#   }
#   Z = log(X)
#   
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   # Centering and scaling #####
#   
#   Z_use <- Z <- as.matrix(Z) # X to-be-used for model fit
#   y_use <- as.vector(y) # y to-be-used for model fit
#   
#   #####
#   Z_mu = colMeans(Z_use)
#   y_mu = mean(y_use)
#   #####
#   
#   # Center X and y if intercept is to be included
#   if (intercept) {
#     Z_use <- Z - matrix(rep(1, times = n), ncol = 1) %*% 
#       matrix(colMeans(Z), nrow = 1)
#     Xb_use = Z %*% U - matrix(rep(1, times = n), ncol = 1) %*% 
#       matrix(colMeans(Z %*% U), nrow = 1)
#     y_use <- y - mean(y)
#   }
#   # scale Xb if scaling = TRUE #####
#   # note: only scale if also centering
#   if(scaling) {
#     Z_sd = apply(Z_use, 2, function(x) sqrt(sum(x^2) / n))
#     # Xb_sd = apply(Xb_use, 2, function(x) sqrt(sum(x^2) / n))
#     Z_use = apply(Z_use, 2, function(x) x / sqrt(sum(x^2) / n))
#   }
#   #####
#   
#   # Generate alpha sequence
#   if (is.null(alpha)) {
#     alpha <- seq(0, 1, len = nalpha)
#   } else {
#     if (min(alpha) < 0 | max(alpha) > 1) stop("alpha range is out of [0, 1].")
#     #alpha <- sort(alpha)
#     nalpha <- length(alpha)
#   }
#   # Generate lambda sequence
#   if(is.null(lambda)) {
#     lambda <- max(abs(t(Z_use) %*% y_use))/n * 
#       exp(seq(0, log(lam.min.ratio), len = nlam))
#   } else {
#     if (min(lambda) < 0) stop("lambda cannot be negative.")
#     #lambda <- sort(lambda)
#     nlam <- length(lambda)
#   }
#   
#   # SVD of (I:-A)
#   if (is.null(Q)) Q <- rare:::svdA(A)
#   # SVD of X
#   E <- rare:::svdX(Z_use, rho)
#   # Two implications up to this point: 1. Q and E will be stored in memory and be
#   # passed in the solver as arguments (which can be problematic when n or p is large).
#   # 2. ADMM iterates at a fixed rho. If varying rho is adopted, E should be updated as well.
#   
#   # Iterate solution paths along alpha
#   beta0 <- beta <- theta0 <- theta <- list()
#   for (i in 1:nalpha) {
#     # if (alpha[i] == 0) {
#     #   # linear log contrasts model ... but where is the constraint on beta? ##
#     #   ret <- glmnet(Z_use, y_use, family = "gaussian", lambda = lambda, 
#     #                 standardize = FALSE, intercept = FALSE, 
#     #                 thresh = min(eps1, eps2), maxit = maxite)
#     #   if(scaling){ # scale back
#     #     beta[[i]] <- as.matrix(ret$beta / Z.sd)
#     #   } else{
#     #     beta[[i]] <- as.matrix(ret$beta)
#     #   }
#     #   beta[[i]] <- as.matrix(ret$beta)
#     #   theta[[i]] <- NA
#     # } else if (alpha[i] == 1) {
#     if (alpha[i] == 1) {
#       # balance regression case!
#       ret <- glmnet(Z_use %*% A, y_use, family = "gaussian", lambda = lambda, 
#                     standardize = FALSE, intercept = FALSE, 
#                     thresh = min(eps1, eps2), maxit = maxite)
#       # ret$beta is actually theta
#       # back-scale
#       theta_tilde = as.matrix(ret$beta)
#       if(scaling){ # scale back
#         # beta[[i]] <- as.matrix(diag(1/Z_sd) %*% A %*% ret$beta)
#         # theta[[i]] <- as.matrix(diag(1 / Z_sd) %*% ret$beta)
#         theta[[i]] <- theta_tilde
#         beta[[i]] <- diag(1 / Z_sd) %*% A %*% theta_tilde
#         # all.equal(unname(as.matrix(diag(1/Z_sd) %*% A %*% ret$beta)), 
#         #           unname(A %*% theta[[i]]))
#       } else{
#         # beta[[i]] <- as.matrix(A %*% ret$beta)
#         # theta[[i]] <- as.matrix(ret$beta)
#         theta[[i]] <- theta_tilde
#         beta[[i]] <- A %*% theta_tilde
#       }
#     } else {
#       # general case when 0 < alpha < 1
#       ret <- rare:::our_solver(Z_use, as.matrix(y_use), Q, E, lambda, alpha[i], rho, 
#                                eps1, eps2, maxite)
#       if(scaling){ # scale back
#         beta[[i]] <- diag(1 / Z_sd) %*% ret$beta
#         theta[[i]] <- diag(1 / Z_sd) %*% ret$gamma # gamma is actually theta
#       } else{
#         beta[[i]] <- ret$beta
#         theta[[i]] <- ret$gamma # gamma is actually theta
#       }
#     }
#     cat(sprintf("Finished model fits for alpha[%s].\n", i))
#   }
#   # Take care of intercept
#   if (intercept) {
#     # beta0 <- lapply(beta, function(b) 
#     #   (sum(y) - c(matrix(colSums(Z), nrow = 1) %*% b))/n)
#     beta0 = lapply(beta, function(b) y_mu - Z_mu %*% b)
#     theta0 = lapply(theta, function(th) y_mu - Z_mu %*% A %*% th)
#   }
#   list(beta0 = beta0, beta = beta, theta0 = theta0, theta = theta, 
#        lambda = lambda, alpha = alpha, A = A, Q = Q, intercept = intercept, 
#        scaling = scaling, btree = btree, U = U)
# }

fitSLRalpha <- function(
  y, X, Xb = NULL, A = NULL, U = NULL, btree = NULL, sbp = NULL, 
  linkage = "complete", rho.type = "squared",
  Q = NULL, intercept = TRUE, lambda = NULL, 
  alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
  rho = 1e-2, eps1 = 1e-7, eps2 = 1e-7, maxite = 1e5, scaling = FALSE
){
  # scaling doesn't work rn.
  # browser()

  if(is.null(btree)) btree = getSupervisedTree(y, X, linkage, rho.type)
  if(is.null(U)) U = getU(btree = btree, sbp = sbp)
  if(is.null(Xb)) Xb = computeBalances(X, btree = btree, sbp = sbp, U = U)
  
  Z = log(X)
  
  if(is.null(A)){
    A = Matrix(U, sparse = TRUE)
  }
  
  n <- nrow(Xb)
  p <- ncol(X)
  
  # Generate alpha sequence
  if (is.null(alpha)) {
    alpha <- seq(0, 1, len = nalpha)
  } else {
    if (min(alpha) < 0 | max(alpha) > 1) stop("alpha range is out of [0, 1].")
    #alpha <- sort(alpha)
    nalpha <- length(alpha)
  }
  
  # Iterate solution paths along alpha
  beta0 <- beta <- theta0 <- theta <- list()
  for (i in 1:nalpha) {
    if (alpha[i] == 1) {
      # balance regression case!
      
      # Centering and scaling #####
      Xb_use <- Xb <- as.matrix(Xb) # X to-be-used for model fit
      Z_use <- Z <- as.matrix(Z)
      y_use <- as.vector(y) # y to-be-used for model fit
      # Center X and y if intercept is to be included
      if (intercept) {
        Z_mu = colMeans(Z)
        Xb_mu = colMeans(Xb_use)
        y_mu = mean(y_use)
        Z_use <- Z - matrix(rep(1, times = n), ncol = 1) %*% 
          matrix(Z_mu, nrow = 1)
        Xb_use <- Xb - matrix(rep(1, times = n), ncol = 1) %*% 
          matrix(Xb_mu, nrow = 1)
        y_use <- y - y_mu
      }
      # scale Xb if scaling = TRUE #####
      # note: only scale if also centering
      if(scaling) {
        Z_sd = apply(Z_use, 2, function(x) sqrt(sum(x^2) / n))
        Xb_sd = apply(Xb_use, 2, function(x) sqrt(sum(x^2) / n))
        Xb_use = apply(Xb_use, 2, function(x) x / sqrt(sum(x^2) / n))
      }
      #####
      
      # Generate lambda sequence
      if(is.null(lambda)) {
        lambda <- max(abs(t(Z_use) %*% y_use))/n * 
          exp(seq(0, log(lam.min.ratio), len = nlam))
      } else {
        if (min(lambda) < 0) stop("lambda cannot be negative.")
        #lambda <- sort(lambda)
        nlam <- length(lambda)
      }
      ret = fitSLRalpha1(
        y = y_use, X = X, Xb = Xb_use, A = A, U = U, btree = btree, 
        linkage = linkage, rho.type = rho.type, 
        intercept = FALSE, scaling = FALSE, # already done above
        lambda = lambda, eps1 = eps1, eps2 = eps2, maxite = maxite)
      # back-scale
      theta_tilde = as.matrix(ret$theta)
      if(scaling){ # scale back
        theta[[i]] <- diag(1 / Xb_sd) %*% theta_tilde
        beta[[i]] <- diag(1 / Z_sd) %*% A %*% theta_tilde
      } else{
        theta[[i]] <- theta_tilde
        beta[[i]] <- A %*% theta_tilde
      }
    } else {
      # general case when 0 <= alpha < 1
      Z = log(X)
      
      # Centering and scaling #####
      Z_use <- Z <- as.matrix(Z) # X to-be-used for model fit
      y_use <- as.vector(y) # y to-be-used for model fit
      # Center X and y if intercept is to be included
      if (intercept) {
        Z_mu = colMeans(Z_use)
        y_mu = mean(y_use)
        Z_use <- Z - matrix(rep(1, times = n), ncol = 1) %*% 
          matrix(Z_mu, nrow = 1)
        y_use <- y - y_mu
      }
      # scale Z if scaling = TRUE #####
      # note: only scale if also centering
      if(scaling) {
        Z_sd = apply(Z_use, 2, function(x) sqrt(sum(x^2) / n))
        Xb_sd = apply(Xb_use, 2, function(x) sqrt(sum(x^2) / n))
        Z_use = apply(Z_use, 2, function(x) x / sqrt(sum(x^2) / n))
      }
      #####
      
      # Generate lambda sequence
      if(is.null(lambda)) {
        lambda <- max(abs(t(Z_use) %*% y_use))/n * 
          exp(seq(0, log(lam.min.ratio), len = nlam))
      } else {
        if (min(lambda) < 0) stop("lambda cannot be negative.")
        #lambda <- sort(lambda)
        nlam <- length(lambda)
      }
      
      # SVD of (I:-A)
      if (is.null(Q)) Q <- rare:::svdA(A)
      # SVD of X
      E <- rare:::svdX(Z_use, rho)
      # Two implications up to this point: 1. Q and E will be stored in memory and be
      # passed in the solver as arguments (which can be problematic when n or p is large).
      # 2. ADMM iterates at a fixed rho. If varying rho is adopted, E should be updated as well.
      
      ret = fitSLRalphaNot1(
        y, X = X, A = A, U = U, btree = btree, linkage = linkage, 
        rho.type = rho.type, Q = Q, E = E, 
        intercept = FALSE, scaling = FALSE, # alrady done above
        lambda = lambda, rho = rho, eps1 = eps1, eps2 = eps2, maxite = maxite
      )
      # ret <- rare:::our_solver(
      #   Z_use, as.matrix(y_use), Q, E, lambda, alpha[i], rho, eps1, eps2, 
      #   maxite)
      if(scaling){ # scale back
        beta[[i]] <- diag(1 / Z_sd) %*% ret$beta
        theta[[i]] <- diag(1 / Xb_sd) %*% ret$gamma # gamma is actually theta
      } else{
        beta[[i]] <- ret$beta
        theta[[i]] <- ret$gamma # gamma is actually theta
      }
    }
    cat(sprintf("Finished model fits for alpha[%s].\n", i))
  }
  # Take care of intercept
  if (intercept) {
    beta0 = lapply(beta, function(b) y_mu - Z_mu %*% b)
    theta0 = lapply(theta, function(th) y_mu - Z_mu %*% A %*% th)
  }
  list(beta0 = beta0, beta = beta, theta0 = theta0, theta = theta, 
       lambda = lambda, alpha = alpha, A = A, Q = Q, intercept = intercept, 
       scaling = scaling, btree = btree, U = U)
}

fitSLRalpha1 <- function(
  y, X, Xb = NULL, A = NULL, U = NULL, btree = NULL, sbp = NULL, linkage = "complete", 
  rho.type = "squared", intercept = TRUE, lambda = NULL, nlam = 50, 
  lam.min.ratio = 1e-4, eps1 = 1e-7, eps2 = 1e-7, maxite = 1e5, 
  scaling = FALSE
){
  # browser()
  # scaling doesn't work rn.
  # browser()
  
  if(is.null(btree)) btree = getSupervisedTree(y, X, linkage, rho.type)
  if(is.null(U)) U = getU(btree = btree, sbp = sbp)
  if(is.null(Xb)) Xb = computeBalances(X, btree = btree, sbp = sbp, U = U)
  
  n <- nrow(Xb)
  p <- ncol(X)
  
  # Centering and scaling #####
  
  Xb_use <- Xb <- as.matrix(Xb) # X to-be-used for model fit
  Z_use <- Z <- as.matrix(Z)
  y_use <- as.vector(y) # y to-be-used for model fit
  
  # Center X and y if intercept is to be included
  if (intercept) {
    Z_mu = colMeans(Z)
    Xb_mu = colMeans(Xb_use)
    y_mu = mean(y_use)
    Z_use <- Z - matrix(rep(1, times = n), ncol = 1) %*% 
      matrix(Z_mu, nrow = 1)
    Xb_use <- Xb - matrix(rep(1, times = n), ncol = 1) %*% 
      matrix(Xb_mu, nrow = 1)
    y_use <- y - y_mu
  }
  # scale Xb if scaling = TRUE #####
  # note: only scale if also centering
  if(scaling) {
    Z_sd = apply(Z_use, 2, function(x) sqrt(sum(x^2) / n))
    Xb_sd = apply(Xb_use, 2, function(x) sqrt(sum(x^2) / n))
    Xb_use = apply(Xb_use, 2, function(x) x / sqrt(sum(x^2) / n))
  }
  #####
  
  # Generate lambda sequence
  if(is.null(lambda)) {
    lambda <- max(abs(t(Z_use) %*% y_use))/n * 
      exp(seq(0, log(lam.min.ratio), len = nlam))
  } else {
    if (min(lambda) < 0) stop("lambda cannot be negative.")
    #lambda <- sort(lambda)
    nlam <- length(lambda)
  }
  
  # Iterate solution paths along alpha
      # balance regression case!
      ret <- glmnet(Xb_use, y_use, family = "gaussian", lambda = lambda, 
                    standardize = FALSE, intercept = FALSE, 
                    thresh = min(eps1, eps2), maxit = maxite)
      # ret$beta is actually theta
      # back-scale
      theta_tilde = as.matrix(ret$beta)
      if(scaling){ # scale back
        theta <- diag(1 / Xb_sd) %*% theta_tilde
        beta <- A %*% theta[[i]]
      } else{
        theta <- theta_tilde
        beta <- A %*% theta_tilde
      }
  # Take care of intercept
  if (intercept) {
    # beta0 <- lapply(beta, function(b) 
    #   (sum(y) - c(matrix(colSums(Z), nrow = 1) %*% b))/n)
    beta0 = y_mu - Z_mu %*% beta
    theta0 = y_mu - Xb_mu %*% theta
  } else{
    beta0 = 0
    theta0 = 0
  }
  list(beta0 = beta0, beta = beta, theta0 = theta0, theta = theta, 
       lambda = lambda, A = A, intercept = intercept, 
       scaling = scaling, btree = btree, U = U)
}

fitSLRalphaNot1 <- function(
  y, X, A = NULL, U = NULL, btree = NULL, sbp = NULL, linkage = "complete", 
  rho.type = "squared", Q = NULL, E = NULL, intercept = TRUE, lambda = NULL,
  nlam = 50, lam.min.ratio = 1e-4, rho = 1e-2, 
  eps1 = 1e-7, eps2 = 1e-7, maxite = 1e5, scaling = FALSE
){
  # scaling doesn't work rn.
  # browser()
  
  if(is.null(btree)) btree = getSupervisedTree(y, X, linkage, rho.type)
  if(is.null(U)) U = getU(btree = btree, sbp = sbp)
  # if(is.null(Xb)) Xb = computeBalances(X, btree = btree, sbp = sbp, U = U)
  
  if(is.null(A)){
    A = Matrix(U, sparse = TRUE)
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Iterate solution paths along alpha
  beta0 <- beta <- theta0 <- theta <- list()
    
      # general case when 0 <= alpha < 1
      Z = log(X)
      
      # Centering and scaling #####
      Z_use <- Z <- as.matrix(Z) # X to-be-used for model fit
      y_use <- as.vector(y) # y to-be-used for model fit
      # Center X and y if intercept is to be included
      if (intercept) {
        Z_mu = colMeans(Z_use)
        y_mu = mean(y_use)
        Z_use <- Z - matrix(rep(1, times = n), ncol = 1) %*% 
          matrix(Z_mu, nrow = 1)
        y_use <- y - y_mu
      }
      # scale Xb if scaling = TRUE #####
      # note: only scale if also centering
      if(scaling) {
        Z_sd = apply(Z_use, 2, function(x) sqrt(sum(x^2) / n))
        # Xb_sd = apply(Xb_use, 2, function(x) sqrt(sum(x^2) / n))
        Z_use = apply(Z_use, 2, function(x) x / sqrt(sum(x^2) / n))
      }
      #####
      
      # SVD of (I:-A)
      if (is.null(Q)) Q <- rare:::svdA(A)
      # SVD of X
      E <- rare:::svdX(Z_use, rho)
      # Two implications up to this point: 1. Q and E will be stored in memory and be
      # passed in the solver as arguments (which can be problematic when n or p is large).
      # 2. ADMM iterates at a fixed rho. If varying rho is adopted, E should be updated as well.
      
      ret <- rare:::our_solver(
        Z_use, as.matrix(y_use), Q, E, lambda, alpha[i], rho, eps1, eps2, 
        maxite)
      if(scaling){ # scale back
        beta <- diag(1 / Z_sd) %*% ret$beta
        theta <- diag(1 / Z_sd) %*% ret$gamma # gamma is actually theta
      } else{
        beta <- ret$beta
        theta <- ret$gamma # gamma is actually theta
      }
  # Take care of intercept
  if (intercept) {
    beta0 = y_mu - Z_mu %*% beta
    theta0 = y_mu - Z_mu %*% A %*% theta
  } else{
    beta0 = 0
    theta0 = 0
  }
  list(beta0 = beta0, beta = beta, theta0 = theta0, theta = theta, 
       lambda = lambda, alpha = alpha, A = A, Q = Q, intercept = intercept, 
       scaling = scaling, btree = btree, U = U)
}












cvSLRalpha <- function(
  y, X, A = NULL, U = NULL, linkage = "complete", rho.type = "squared",
  Q = NULL, intercept = TRUE, lambda = NULL, 
  alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
  rho = 1e-2, eps1 = 1e-7, eps2 = 1e-7, maxite = 1e5, nfolds = 5, 
  foldid = NULL, scaling = FALSE
){
  # scaling doesn't work rn.
  # browser()
  
  # fit the models
  fitObj = fitSLRalpha(
    y = y, X = X, A = A, U = U, linkage = linkage, rho.type = rho.type, 
    Q = Q, intercept = intercept, lambda = lambda, alpha = alpha, nlam = nlam, 
    lam.min.ratio = lam.min.ratio, nalpha = nalpha, rho = rho, 
    eps1 = eps1, eps2 = eps2, maxite = maxite, scaling = scaling)
  btree = fitObj$btree
  errtype = "mean-squared-error"
  Xb = computeBalances(X, btree)
  
  n <- length(y)
  nlam <- length(fitObj$lambda)
  nalpha <- length(fitObj$alpha)
  errs <- array(NA, dim=c(nlam, nalpha, nfolds))
  
  # define error function
  errfun <- function(est, truth) colMeans((est - truth)^2)
  # if (errtype == "mean-absolute-error") {
  #   errfun <- function(est, truth) colMeans(abs(est - truth))
  # } else if (errtype != "mean-squared-error") {
  #   stop("The error function needs to be either mean squared error or mean absolute error.")
  # }
  
  # make folds, if necessary
  if(is.null(foldid)){
    nn <- round(n / nfolds)
    sizes <- rep(nn, nfolds)
    sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
    b <- c(0, cumsum(sizes))
    # set.seed(100) # set.seed for random number generator
    ii <- sample(n)
    folds <- list()
    for (i in seq(nfolds)) folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  } else{
    folds = list()
    for(i in 1:nfolds) folds[[i]] = which(foldid == i)
  }
  
  # Fit based on folds and compute error metric
  for (i in seq(nfolds)) {
    # fit model on all but the ith fold
    fit_cv <- fitSLRalpha(
      y = y[-folds[[i]]], X = X[-folds[[i]], ], Xb[-folds[[i]], ],
      A = fitObj$A, U = fitObj$U, 
      btree = fitObj$btree, 
      linkage = linkage, rho.type = rho.type, 
      Q = fitObj$Q, intercept = fitObj$intercept, lambda = fitObj$lambda, 
      alpha = fitObj$alpha, nlam = nlam, lam.min.ratio = lam.min.ratio, 
      nalpha = nalpha, rho = rho, eps1 = eps1, eps2 = eps2, maxite = maxite,
      scaling = fitObj$scaling)
    pred_te <- lapply(seq(nalpha), function(k) {
      if (fitObj$intercept) {
        # log(X)[folds[[i]], ] %*% fit_cv$beta[[k]] +
        #   rep(fit_cv$beta0[[k]], each = length(folds[[i]]))
        Xb[folds[[i]], ] %*% fit_cv$theta[[k]] +
          rep(fit_cv$theta0[[k]], each = length(folds[[i]]))
      } else {
        # X[folds[[i]], ] %*% fit_cv$beta[[k]]
        Xb[folds[[i]], ] %*% fit_cv$theta[[k]]
      }
    })
    for (k in seq(nalpha)) errs[, k, i] <- errfun(pred_te[[k]], y[folds[[i]]])
    cat("##########################\n")
    cat(sprintf("Finished model fits for fold[%s].\n", i))
    cat("##########################\n")
  }
  m <- apply(errs, c(1, 2), mean)
  se <- apply(errs, c(1, 2), stats::sd) / sqrt(nfolds)
  ibest <- which(m == min(m), arr.ind = TRUE)[1, , drop = FALSE]
  
  # list (folds = folds, errs = errs, m = m, se = se, ibest = ibest,
  #       lambda.best = fitObj$lambda[ibest[1]], alpha.best = fitObj$alpha[ibest[2]])
  return(list(
    beta0 = fitObj$beta0,
    beta = fitObj$beta,
    that0 = fitObj$theta0,
    theta = fitObj$theta,
    lambda = fitObj$lambda,
    fits = fitObj,
    btree = btree,
    U = U,
    cvm = m,
    min.idx = ibest
  ))
}
# 
# rarefit.predict <- function(fitObj, cvObj, newx) {
#   ibest.lambda <- cvObj$ibest[1]
#   ibest.alpha <- cvObj$ibest[2]
#   if (fitObj$intercept) {
#     as.vector(newx %*% fitObj$beta[[ibest.alpha]][, ibest.lambda] + 
#                 fitObj$beta0[[ibest.alpha]][ibest.lambda])
#   } else {
#     as.vector(newx %*% fitObj$beta[[ibest.alpha]][, ibest.lambda])
#   }
# }