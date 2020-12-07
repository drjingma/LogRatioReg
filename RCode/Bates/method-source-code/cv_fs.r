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
    for(k in 1:k_max) {
      model <- regsubsets(x[fold_id != i, ], y[fold_id != i], method = "forward", nvmax = k_max, intercept = FALSE)
      predictions <- predict.regsubsets(model, x[fold_id == i, ], id = k)
      cvm[k] <-  cvm[k] + apply(predictions, 2, function(x){ mean((x - y[fold_id == i])^2) })
    }
  }
  
  cvm <- cvm / n_folds
  best <- which.min(cvm)
  
  #extract best beta
  full_model <- regsubsets(x, y, method = "forward", nvmax = best, intercept = FALSE)
  beta_best <- 1:p * 0
  coefi=coef(full_model, id=best)
  beta_best[as.numeric(names(coefi))] <- coefi
  
  return(list(cvm = cvm, beta = beta_best, model = full_model))
}

predict.regsubsets = function (object , newdata , id,...){
  coefi=coef(object ,id=id)
  xvars=names(coefi)
  beta_temp <- 1:p * 0
  beta_temp[as.numeric(xvars)] <- coefi
  
  newdata %*% beta_temp
}


# test code 
# n <- 100
# p <- 30
# 
# c <- 0
# sigma <- autoRegressiveCorr(p, c)
# beta <- c(2,-2,1,-1, rep(0, p-4))
# X  <- abs(mvrnorm(n = n, rep(0, p), sigma))
# x <- X
# z <- log(X)
# y <- z %*% beta * sig + rnorm(n)
# y <- y - mean(y)
# z <- scale(z, center = TRUE, scale = FALSE)
# 
# full_model <- regsubsets(z, y, method = "forward", nvmax = 10, intercept = FALSE)
# predict.regsubsets(full_model, z, 6)
# out <- cv.fs(z, y, k_max = 10, n_folds = 10)
