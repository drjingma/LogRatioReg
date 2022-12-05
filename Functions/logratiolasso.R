# log-ratio lasso method, code issues
# fix k error
cv_two_stage <- function(z, y, family = "gaussian", lambda_1 = NULL, k_max = 10,
                         n_folds = 10, nlambda = 20, folds = NULL, gamma = 0) {
  p <- ncol(z)
  #returns list of length lambda of vectors of size k_max with the prediction error.
  two_step_obj <- two_stage(z, y, family = family, lambda_1 = lambda_1, 
                            k_max = k_max, nlambda = nlambda)
  lambda_1 <- two_step_obj$lambda
  nlambda = length(lambda_1) ############################################################
  
  #Create equally size folds
  if(is.null(folds)) {
    folds <- sample(cut(seq(1,nrow(z)),breaks=n_folds,labels=FALSE))
  } else {
    n_folds = max(folds) - min(folds) + 1
  }
  
  #Perform k-fold cross validation
  mse <- list()
  for(i in 1:n_folds){
    print(paste0("Starting CV fold ", i))
    #Segement your data by fold using the which() function
    test_indices <- which(folds==i,arr.ind=TRUE)
    out <- predict_two_stage(
      z = z[-test_indices, ], y = y[-test_indices], new_z = z[test_indices, ],
      lambda_1 = lambda_1, k_max = k_max, nlambda = nlambda, family = family )
    
    if(family == "gaussian") {
      mse_fun <- function(y_pred) {
        if(!is.null(ncol(y_pred))) {
          if(ncol(y_pred > 1)) {
            return(apply(y_pred,2,function(yp){sum((yp - y[test_indices])**2)}))
          }
        }
        sum((yp - y[test_indices])**2)
      }
    } else if (family == "binomial") {
      mse_fun <- function(y_pred) {
        if(!is.null(ncol(y_pred))) {
          if(ncol(y_pred > 1)) {
            return(apply(y_pred,2,function(yp){sum(log(1 + exp(-(2*y[test_indices] - 1) * y_pred)))}))
          }
        }
        sum(log(1 + exp(-(2*y[test_indices] - 1) * y_pred)))
      }
    }
    mse[[i]] <- lapply(out$y_pred, mse_fun)
    
    #glmnet may return fewer than nlambda entries if convergence fails
    while(length(mse[[i]]) < nlambda) {
      ##########################################################################
      mse[[i]][[length(mse[[i]]) + 1]] = rep(Inf, k_max) # rep(Inf, k), k undefined
      ##########################################################################
    }
  }
  
  # # old code
  # mse_full <- lapply(mse[[1]], function(x){x / n_folds})
  # #print(length(mse_full))
  # for(i in 2:n_folds) { # this looks like it's averaging over the folds, but why like this?
  #   for(j in 1:length(mse[[1]])) {
  #     mse_full[[j]] <- mse_full[[j]] + mse[[i]][[j]] / n_folds
  #   }
  # }
  # # fix min k, use 1se rule to choose lamdba?
  # # for each k, find the lambda that minimizes (or 1se from minimizer of) mse
  # mse_full <- simplify2array(mse_full)
  # best <- which(mse_full == min(mse_full), arr.ind = TRUE) # minimizing lambda and number of log-ratios (1, ..., k_max)
  # lambda_min <- best[1,2]
  # k_min <- best[1,1]
  # beta_min <- out_to_beta(two_step_obj$coef[[lambda_min]], k_max, p)[, k_min]
  
  # 1se rule
  mse_full0 = array(NA, dim = c(k_max, nlambda, n_folds))
  for(i in 1:n_folds){
    mse_full0[, , i] = matrix(unlist(mse[[i]]), nrow = k_max, ncol = nlambda)
  }
  mse_full_mean = apply(
    mse_full0, c(1, 2), function(x) base::mean(x, na.rm = TRUE))
  mse_full_se = apply(
    mse_full0, c(1, 2), 
    function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  best = which(
    mse_full_mean <= min(mse_full_mean, na.rm = TRUE) + 
      gamma * 
      mse_full_se[which(
        mse_full_mean == min(mse_full_mean, na.rm = TRUE), arr.ind = TRUE)][1], 
    arr.ind = TRUE)
  
  # min
  best0 = which(mse_full_mean == min(mse_full_mean), arr.ind = TRUE)
  lambda_min <- best0[1,2]
  k_min <- best0[1,1]
  # pick the largest lambda with this k_min value (has the 1st index, due to ordering lambda largest -> smallest)
  k_gammase = k_min
  lambda_gammase = best[best[, "row"] == k_gammase, , drop = FALSE][1, 2]
  
  beta_min <- out_to_beta(two_step_obj$coef[[lambda_min]], k_max, p)[, k_gammase]
  beta_gammase <- out_to_beta(two_step_obj$coef[[lambda_gammase]], k_max, p)[, k_gammase]
  
  # list(mse = mse_full, best_params = best, lambda_min = lambda_min, k_min = k_min,
  #      two_step_obj = two_step_obj, beta_min = beta_min, lambda = lambda_1)
  list(
    mse = mse_full_mean, min_params = best0, gammase_params = best, 
    lambda_min = lambda_min, k_min = k_min,
    two_step_obj = two_step_obj, 
    beta_min = beta_min, 
    beta_gammase = beta_gammase, 
    lambda = lambda_1)
}
