### scale, unlike the original
### reminder: scaling doesn't affect S.hat
# fitSLR2new <- function(
#   y, X, A = NULL, U = NULL, linkage = "complete", rho.type = "squared",
#   Q = NULL, intercept = TRUE, lambda = NULL, 
#   alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
#   rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5, maxite = 1e6
# ){
#   btree = getSupervisedTree(y, X, linkage, rho.type)
#   U = getU(btree = btree)
#   if(is.null(A)){
#     A = Matrix(U, sparse = TRUE)
#   }
#   Z = log(X)
#   
#   n <- nrow(X)
#   p <- ncol(X)
#   
#   Z_use <- Z <- as.matrix(Z) # X to-be-used for model fit
#   y_use <- as.vector(y) # y to-be-used for model fit
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
#   beta0 <- beta <- theta <- list()
#   for (i in 1:nalpha) {
#     if (alpha[i] == 0) {
#       # linear log contrasts model ... but where is the constraint on beta? ###########################
#       ret <- glmnet(Z_use, y_use, family = "gaussian", lambda = lambda, 
#                     standardize = TRUE, intercept = intercept, 
#                     thresh = min(eps1, eps2), maxit = maxite)
#       beta0[[i]] = ret$a0
#       beta[[i]] <- as.matrix(ret$beta)
#       theta[[i]] <- NA
#     } else if (alpha[i] == 1) {
#       # balance regression case!
#       ret <- glmnet(Z_use %*% U, y_use, family = "gaussian", lambda = lambda, 
#                     standardize = TRUE, intercept = intercept, 
#                     # penalty.factor = c(rep(1, p-1), 0), # not needed here
#                     thresh = min(eps1, eps2), maxit = maxite)
#       beta0[[i]] = ret$a0
#       beta[[i]] <- as.matrix(A %*% ret$beta) # ret$beta is actually theta
#       theta[[i]] <- as.matrix(ret$beta)
#     } else {
#       # general case when 0 < alpha < 1
#       ret <- rare:::our_solver(Z_use, as.matrix(y_use), Q, E, lambda, alpha[i], rho, 
#                                eps1, eps2, maxite)
#       # beta0[[i]] = a0 ???????????????????
#       beta[[i]] <- ret$beta
#       theta[[i]] <- ret$gamma # gamma is actually theta
#     }
#     cat(sprintf("Finished model fits for alpha[%s].\n", i))
#   }
#   list(a0 = beta0, beta = beta, theta = theta, lambda = lambda, 
#        alpha = alpha, A = A, Q = Q, intercept = intercept, btree = btree, U = U)
# }
# 
# cvSLR2new <- function(
#   y, X, A = NULL, U = NULL, linkage = "complete", rho.type = "squared",
#   Q = NULL, intercept = TRUE, lambda = NULL, 
#   alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
#   rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5, maxite = 1e6, nfolds = 5, 
#   foldid = NULL
# ){
#   
#   fitObj = fitSLR2new(
#     y, X, A = A, U = U, linkage = linkage, rho.type = rho.type, 
#     Q = Q, intercept = intercept, lambda = lambda, 
#     alpha = alpha, nlam = nlam, lam.min.ratio = lam.min.ratio, nalpha = nalpha,
#     rho = rho, eps1 = eps1, eps2 = eps2, maxite = maxite)
#   btree = fitObj$btree
#   errtype = "mean-squared-error"
#   
#   n <- length(y)
#   nlam <- length(fitObj$lambda)
#   nalpha <- length(fitObj$alpha)
#   errs <- array(NA, dim=c(nlam, nalpha, nfolds))
#   
#   # define error function
#   errfun <- function(est, truth) colMeans((est - truth)^2)
#   if (errtype == "mean-absolute-error") {
#     errfun <- function(est, truth) colMeans(abs(est - truth))
#   } else if (errtype != "mean-squared-error") {
#     stop("The error function needs to be either mean squared error or mean absolute error.")
#   }
#   
#   # make folds
#   if(is.null(foldid)){
#     nn <- round(n / nfolds)
#     sizes <- rep(nn, nfolds)
#     sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
#     b <- c(0, cumsum(sizes))
#     # set.seed(100) # set.seed for random number generator
#     ii <- sample(n)
#     folds <- list()
#     for (i in seq(nfolds)) folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
#   } else{
#     folds = list()
#     for(i in 1:nfolds) folds[[i]] = which(foldid == i)
#   }
#   
#   # Fit based on folds and compute error metric
#   for (i in seq(nfolds)) {
#     # fit model on all but the ith fold
#     fit_cv <- fitSLR2(y = y[-folds[[i]]], X = X[-folds[[i]], ], A = fitObj$A, Q = fitObj$Q,
#                       intercept = fitObj$intercept, lambda = fitObj$lambda, alpha = fitObj$alpha, rho.type = rho.type)
#     pred_te <- lapply(seq(nalpha), function(k) {
#       if (fitObj$intercept) {
#         X[folds[[i]], ] %*% fit_cv$beta[[k]] + rep(fit_cv$a0[[k]], each = length(folds[[i]]))
#       } else {
#         X[folds[[i]], ] %*% fit_cv$beta[[k]]
#       }
#     })
#     for (k in seq(nalpha)) errs[, k, i] <- errfun(pred_te[[k]], y[folds[[i]]])
#     cat("##########################\n")
#     cat(sprintf("Finished model fits for fold[%s].\n", i))
#     cat("##########################\n")
#   }
#   m <- apply(errs, c(1, 2), mean)
#   se <- apply(errs, c(1, 2), stats::sd) / sqrt(nfolds)
#   ibest <- which(m == min(m), arr.ind = TRUE)[1, , drop = FALSE]
#   
#   # list (folds = folds, errs = errs, m = m, se = se, ibest = ibest,
#   #       lambda.best = fitObj$lambda[ibest[1]], alpha.best = fitObj$alpha[ibest[2]])
#   return(list(
#     int = fitObj$a0,
#     bet = fitObj$beta,
#     thet = fitObj$theta,
#     lambda = fitObj$lambda,
#     fits = fitObj,
#     btree = btree,
#     U = U,
#     cvm = m,
#     min.idx = ibest
#   ))
# }










fitSLRalpha <- function(
  y, X, A = NULL, U = NULL, linkage = "complete", rho.type = "squared",
  Q = NULL, intercept = TRUE, lambda = NULL, 
  alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
  rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5, maxite = 1e6, scaling = FALSE
){
  btree = getSupervisedTree(y, X, linkage, rho.type)
  U = getU(btree = btree)
  if(is.null(A)){
    A = Matrix(U, sparse = TRUE)
  }
  Z = log(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  Z_use <- Z <- as.matrix(Z) # X to-be-used for model fit
  y_use <- as.vector(y) # y to-be-used for model fit
  # Center X and y if intercept is to be included
  if (intercept) {
    Z_use <- Z - matrix(rep(1, times = n), ncol = 1) %*% 
      matrix(colMeans(Z), nrow = 1)
    # scale X if scaling = TRUE
    if(scaling) {
      Z.sd <- apply(Z_use, 2, sd)
      Z_use <- scale(Z_use, center=FALSE, scale = Z.sd)
    }
    y_use <- y - mean(y)
  }
  
  # Generate alpha sequence
  if (is.null(alpha)) {
    alpha <- seq(0, 1, len = nalpha)
  } else {
    if (min(alpha) < 0 | max(alpha) > 1) stop("alpha range is out of [0, 1].")
    #alpha <- sort(alpha)
    nalpha <- length(alpha)
  }
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
  
  # Iterate solution paths along alpha
  beta0 <- beta <- theta <- list()
  for (i in 1:nalpha) {
    if (alpha[i] == 0) {
      # linear log contrasts model ... but where is the constraint on beta? ###########################
      ret <- glmnet(Z_use, y_use, family = "gaussian", lambda = lambda, 
                    standardize = FALSE, intercept = FALSE, 
                    thresh = min(eps1, eps2), maxit = maxite)
      if(scaling){ # scale back
        beta[[i]] <- as.matrix(ret$beta / Z.sd)
      } else{
        beta[[i]] <- as.matrix(ret$beta)
      }
      beta[[i]] <- as.matrix(ret$beta)
      theta[[i]] <- NA
    } else if (alpha[i] == 1) {
      # balance regression case!
      ret <- glmnet(Z_use %*% U, y_use, family = "gaussian", lambda = lambda, 
                    standardize = FALSE, intercept = FALSE, 
                    # penalty.factor = c(rep(1, p-1), 0), # not needed here
                    thresh = min(eps1, eps2), maxit = maxite)
      # ret$beta is actually theta
      if(scaling){ # scale back
        beta[[i]] <- as.matrix(A %*% (ret$beta / Z.sd))
      } else{
        beta[[i]] <- as.matrix(A %*% ret$beta)
      }
      theta[[i]] <- as.matrix(ret$beta)
    } else {
      # general case when 0 < alpha < 1
      ret <- rare:::our_solver(Z_use, as.matrix(y_use), Q, E, lambda, alpha[i], rho, 
                        eps1, eps2, maxite)
      if(scaling){ # scale back
        beta[[i]] <- ret$beta / Z.sd
      } else{
        beta[[i]] <- ret$beta
      }
      theta[[i]] <- ret$gamma # gamma is actually theta
    }
    cat(sprintf("Finished model fits for alpha[%s].\n", i))
  }
  # Take care of intercept
  if (intercept) {
    beta0 <- lapply(beta, function(b) 
      (sum(y) - c(matrix(colSums(Z), nrow = 1) %*% b))/n)
  }
  list(a0 = beta0, beta = beta, theta = theta, lambda = lambda, 
       alpha = alpha, A = A, Q = Q, intercept = intercept, scaling = scaling,
       btree = btree, U = U)
}

cvSLRalpha <- function(
  y, X, A = NULL, U = NULL, linkage = "complete", rho.type = "squared",
  Q = NULL, intercept = TRUE, lambda = NULL, 
  alpha = NULL, nlam = 50, lam.min.ratio = 1e-4, nalpha = 10,
  rho = 1e-2, eps1 = 1e-6, eps2 = 1e-5, maxite = 1e6, nfolds = 5, 
  foldid = NULL, scaling = FALSE
){
  
  fitObj = fitSLRalpha(
    y, X, A = A, U = U, linkage = linkage, rho.type = rho.type, 
    Q = Q, intercept = intercept, lambda = lambda, 
    alpha = alpha, nlam = nlam, lam.min.ratio = lam.min.ratio, nalpha = nalpha,
    rho = rho, eps1 = eps1, eps2 = eps2, maxite = maxite, scaling = scaling)
  btree = fitObj$btree
  errtype = "mean-squared-error"
  
  n <- length(y)
  nlam <- length(fitObj$lambda)
  nalpha <- length(fitObj$alpha)
  errs <- array(NA, dim=c(nlam, nalpha, nfolds))
  
  # define error function
  errfun <- function(est, truth) colMeans((est - truth)^2)
  if (errtype == "mean-absolute-error") {
    errfun <- function(est, truth) colMeans(abs(est - truth))
  } else if (errtype != "mean-squared-error") {
    stop("The error function needs to be either mean squared error or mean absolute error.")
  }
  
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
    fit_cv <- fitSLRalpha(y = y[-folds[[i]]], X = X[-folds[[i]], ], A = fitObj$A, Q = fitObj$Q,
                      intercept = fitObj$intercept, lambda = fitObj$lambda, alpha = fitObj$alpha, 
                      scaling = scaling, rho.type = rho.type)
    pred_te <- lapply(seq(nalpha), function(k) {
      if (fitObj$intercept) {
        X[folds[[i]], ] %*% fit_cv$beta[[k]] + rep(fit_cv$a0[[k]], each = length(folds[[i]]))
      } else {
        X[folds[[i]], ] %*% fit_cv$beta[[k]]
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
    int = fitObj$a0,
    bet = fitObj$beta,
    thet = fitObj$theta,
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