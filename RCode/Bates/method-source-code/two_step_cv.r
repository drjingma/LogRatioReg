short_to_long <- function(i, j, p) {
  #take pairwise coordinates and return long-form index
  k <- 1
  for(a in 1:p) {
    for(b in (a+1):p) {
      if(a == i & b == j) {
        return(k)
      }
      k <- k + 1
    }
  }
}

custom_fs <- function(expanded_set, resp, k_max, selected_vars, p) {
  if(k_max <= 0) {
    return(0)
  }
  subs <- myfs(expanded_set, resp, nsteps = k_max, center = FALSE)

  list(theta_vals = subs$beta, theta_ind = subs$ind, subset = which(selected_vars))
}


two_step_generic <- function(z, y, k_max = NULL, lambda_1 = lambda_1, second.stage = c("y", "yhat")) {
  p <- ncol(x)
  
  constrained_fit <- glmnet.constr(z, y, family = "gaussian", lambda = lambda_1)
  betas <- constrained_fit$beta
  #print(betas)
  #print(lambda_1)
  selected_vars <- apply(betas, 2, function(x){(abs(x) > 0)})

  output <- list()
  for (i in 1:length(lambda_1)) {
    if(sum(selected_vars[, i] > 0)) {
      #print(paste0("which: i",i))
      #print(selected_vars)
      expanded_set <- small_to_big_z(z, which(selected_vars[, i]))
      y_hat <- z %*% betas[, i]
      resp <- y  
      if(second.stage == "yhat") {
        resp <- y_hat
      }
      
      d <- sum(selected_vars[, i]) - 1
      k <- min(d, k_max)
  
      output[[i]] <- custom_fs(expanded_set, resp, k, selected_vars[, i], p)
      #print(output[[i]])
    } else {
      output[[i]] <- NA
    }
  }
  
  list(coef = output, lambda = lambda_1, selected_vars = selected_vars)
}


out_to_beta <- function(obj, k, p) {
  #takes output of "two_step_generic$coef[[i]]" and returns p by k matrix of beta coefficients
  betas <- matrix(0, p, k)
  
  if(is.na(obj)) {
    return(betas)
  } else if (class(obj) != "list") {
    return(betas)
  }
  #print(paste0("obj", obj))
  d <- min(k, nrow(obj$theta_ind))
  if(d == 0) {
    return(betas)
  }
  
  # print(obj$subset)
  # print(obj$theta_ind)
  # print(paste0("d:", d))
  # print(dim(betas))
  # print(paste0("p:", p))
  
  for(i in 1:d) {
    beta <- rep(0, p)
    
    # print(i)
    # print(obj$subset[obj$theta_ind[i,1]])
    # print(obj$subset[obj$theta_ind[i,2]])
    # print(obj$theta_vals[[i]])

    for(j in 1:i){ 
      beta[obj$subset[obj$theta_ind[j,1]]] <- beta[obj$subset[obj$theta_ind[j,1]]] + 
        obj$theta_vals[[i]][j]
      beta[obj$subset[obj$theta_ind[j,2]]] <- beta[obj$subset[obj$theta_ind[j,2]]] - 
        obj$theta_vals[[i]][j]
    }
    
    # print(paste0("i completed"))
    betas[, i] <- beta
  }
  #fill in remaining entries with same beta
  if(d < k) {
    for(i in (d+1):k) {
      betas[, i] <- betas[, d]
    }
  }
  betas
}


predict_two_step_generic <- function(z, y, new_z, lambda_1 = NULL, k_max = 5, ...) {
  p <- ncol(z)
  fit <- two_step_generic(z, y, lambda_1 = lambda_1, k_max = k_max, ...)
  betas <- lapply(fit$coef, function(obj){out_to_beta(obj, k_max, p)})
  
  predictor <- function(betas) {
    if(!is.null(ncol(betas))) {
      return(apply(betas,2,function(b){new_z %*% b}))
    }
    return(matrix(rep(new_z %*% betas, k_max), ncol = k_max, byrow = FALSE))
  }
  
  y_pred <- lapply(betas, predictor)
  #list of n_pred by k matrix of predictions
  
  out <- list(two_step_obj = fit, y_pred = y_pred)
  out
}


cv_two_step_generic <- function(z, y, lambda_1 = NULL, k_max = 5, n_folds = 10, ...) {
  #returns list of length lambda of vectors of size k_max with the prediction error.
  two_step_obj <- two_step_generic(z, y, lambda_1 = lambda_1, k_max = k_max, ...)
  
  #Create 10 equally size folds
  folds <- sample(cut(seq(1,nrow(z2)),breaks=n_folds,labels=FALSE))

  #Perform 10 fold cross validation
  mse <- list()
  for(i in 1:n_folds){
      #print(paste0("entering fold", i))
      #Segement your data by fold using the which() function 
      testIndexes <- which(folds==i,arr.ind=TRUE)
      out <- predict_two_step_generic(z[-testIndexes, ], y[-testIndexes], new_z = z[testIndexes, ],
                                      lambda_1 = lambda_1, k_max = k_max, ...)
      
      mse_fun <- function(y_pred) {
        if(!is.null(ncol(y_pred))) {
          if(ncol(y_pred > 1)) {
            return(apply(y_pred,2,function(yp){sum((yp - y[testIndexes])**2)}))
          }
        }
        sum((y_pred - y[testIndexes])**2)
      }
      mse[[i]] <- lapply(out$y_pred, mse_fun)
  }
  
  
  mse_full <- lapply(mse[[1]], function(x){x / n_folds})
  #print(length(mse_full))
  for(i in 1:n_folds) {
    for(j in 1:length(mse[[1]])) {
      mse_full[[j]] <- mse_full[[j]] + mse[[i]][[j]] / n_folds
    }
  }
  mse_full <- simplify2array(mse_full)
  best <- which(mse_full == min(mse_full), arr.ind = TRUE)
  lambda_min <- best[1,2]
  k_min <- best[1,1]
  beta_min <- out_to_beta(two_step_obj$coef[[lambda_min]], k_max, p)[, k_min]
  
  list(mse = mse_full, best = best, lambda_min = lambda_min, k_min = k_min, 
       two_step_obj = two_step_obj, beta_min = beta_min)
}
