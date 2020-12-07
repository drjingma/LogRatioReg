approximate_fs <- function(x, y, k_max = 10) {
  # Implements approximate forward stepwise algorithm
  # does not assume centered or scaled response
  
  yc <- y - mean(y)
  x_norm <- scale(x)
  p <- ncol(x)
  n <- nrow(x)
  
  #store the ratios
  ratios <- matrix(0, nrow = 2, ncol = k_max) #
  thetas <- matrix(0, nrow = k_max, ncol = k_max) #beta for each ratio at each step
  intercepts <- rep(0, k_max)

  x_current <- rep(1, n) #current covariate matrix. initially just an intercept
  
  for(k in 1:k_max) {
    #this loop can be optimized to do half as many matrix inversions.
    
    #find the best ratio
    cors <- t(x_norm) %*% yc 
    big <- which.max(cors)
    small <- which.min(cors)
    ratios[, k] <- c(big, small)
    #regress out earlier features
    x_new_norm <- (x_norm[, big] - x_norm[, small])
    if(k == 1){
      x_current_norm <- x_new_norm
    } else {
      x_current_norm <- cbind(x_current_norm, x_new_norm)
    }
    # full hat matrix
    hat_matrix <- x_current_norm %*% solve(t(x_current_norm) %*% x_current_norm) %*% t(x_current_norm)
    yc <- yc - hat_matrix %*% yc
    #hat matrix of most recent feature
    hat_small <- x_new_norm %*% solve(t(x_new_norm) %*% x_new_norm) %*% t(x_new_norm) 
    x_norm <- x_norm - hat_small %*% x_norm
      
    #add this feature to the model
    x_new <- x[, big] - x[, small]
    x_current <- cbind(x_current, x_new)
    
    #find regression coefficients on orginal scale (as currently written, this is a lot redundant computation)
    theta_current <- solve(t(x_current) %*% x_current) %*% t(x_current) %*% y
    thetas[1:k, k] <- theta_current[-1] 
    intercepts[k] <- theta_current[1]
  }
  
  #represent beta in the normal `beta` representation
  betas <- matrix(0, nrow = p, ncol = k_max)
  for(k in 1:k_max) {
    for(i in 1:k) {
      betas[ratios[1, i], k] <- betas[ratios[1, i], k] + thetas[i, k]
      betas[ratios[2, i], k] <- betas[ratios[2, i], k] - thetas[i, k]
    }
  }
  
  list(theta = thetas, ratios = ratios, beta = betas, intercept = intercepts)
}

# test case
#n <- 100
#p <- 20
#set.seed(314159)
#x <- matrix(rnorm(n *p), nrow = n, ncol = p)
#y <- .2 + x[, 1] - x[, 2] + .5 * x[, 3] - .5 * x[, 4] + rnorm(n, 0, .1)
#out <- approximate_fs(x, y, k_max = 5)



predict.approximate_fs <- function(fs_model, x) {
  #arguments: fs_model: output from approximate_fs
  #x: matrix of covariates
  
  predictions <- x %*% fs_model$beta + fs_model$intercept
}
#out2 <- predict.approximate_fs(out, x[1:20, ])
#apply(out2, 2, function(x){x - y[1:20]}) #check the residuals



cv.approximate_fs <- function(x, y, k_max = 5, n_folds = 10, fold_id = NULL, ...) {
  # performs k-fold CV to determine optimal number of steps for approximate fs
  # fold_id must be integers 1,2,...,k if used
  
  if(is.null(fold_id)) {
    fold_id <- sample(1:length((y)) %/% n_folds + 1, length(y), replace = FALSE)
  } else {
    n_folds <- max(fold_id)
  }
  
  cvm <- rep(0, k_max)
  for(i in 1:n_folds) {
    model <- approximate_fs(x[fold_id != i, ], y[fold_id != i], k_max, ...)
    predictions <- predict.approximate_fs(model, x[fold_id == i, ])
    cvm <-  cvm + apply(predictions, 2, function(x){ mean((x - y[fold_id == i])^2) })
  }
  
  cvm <- cvm / n_folds
  best <- which.min(cvm)

  full_model <- approximate_fs(x, y, k_max = best)
  beta_best <- full_model$beta[,best]

  return(list(cvm = cvm, beta = beta_best, model = full_model))
}
#out3 <- cv.approximate_fs(x, y, k_max = 10 ) 
#out3$cvm
#which.min(out3$cvm)


cv.fs <- function(x, y, k_max = 5, n_folds = 10, fold_id = NULL, ...) {
  # performs k-fold CV to determine optimal number of steps for approximate fs
  # fold_id must be integers 1,2,...,k if used
  
  if(is.null(fold_id)) {
    fold_id <- sample(1:length((y)) %/% n_folds + 1, length(y), replace = FALSE)
  } else {
    n_folds <- max(fold_id)
  }
  
  cvm <- rep(0, k_max)
  for(i in 1:n_folds) {
    model <- approximate_fs(x[fold_id != i, ], y[fold_id != i], k_max, ...)
    predictions <- predict.approximate_fs(model, x[fold_id == i, ])
    cvm <-  cvm + apply(predictions, 2, function(x){ mean((x - y[fold_id == i])^2) })
  }
  
  cvm <- cvm / n_folds
  best <- which.min(cvm)

  full_model <- approximate_fs(x, y, k_max = best)
  beta_best <- full_model$beta[,best]

  return(list(cvm = cvm, beta = beta_best, model = full_model))
}
  