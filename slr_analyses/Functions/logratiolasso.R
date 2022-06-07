# log-ratio lasso method, code issues
# fix k error
cv_two_stage <- function(z, y, family = "gaussian", lambda_1 = NULL, k_max = 10,
                         n_folds = 10, nlambda = 20, folds = NULL) {
  p <- ncol(z)
  #returns list of length lambda of vectors of size k_max with the prediction error.
  two_step_obj <- two_stage(z, y, family = family, lambda_1 = lambda_1, 
                            k_max = k_max, nlambda = nlambda)
  lambda_1 <- two_step_obj$lambda
  
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
  
  mse_full <- lapply(mse[[1]], function(x){x / n_folds})
  #print(length(mse_full))
  for(i in 2:n_folds) {
    for(j in 1:length(mse[[1]])) {
      mse_full[[j]] <- mse_full[[j]] + mse[[i]][[j]] / n_folds
    }
  }
  mse_full <- simplify2array(mse_full)
  best <- which(mse_full == min(mse_full), arr.ind = TRUE)
  lambda_min <- best[1,2]
  k_min <- best[1,1]
  beta_min <- out_to_beta(two_step_obj$coef[[lambda_min]], k_max, p)[, k_min]
  
  list(mse = mse_full, best_params = best, lambda_min = lambda_min, k_min = k_min,
       two_step_obj = two_step_obj, beta_min = beta_min, lambda = lambda_1)
}