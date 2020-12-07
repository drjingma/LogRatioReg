wide_sparsenet <- function(z, y) {
  p <- ncol(z)
  big_z <- small_to_big_z(z, 1:p)
  s_fit <- cv.sparsenet(big_z, y)

  theta_to_beta_v2(coef(s_fit)[-1], 1:p, p)
}

two_step <- function(x, y) {
  #two step procedure for gaussian case
  #step 1: constrained lasso
  #step 2: sparsenet

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
	
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- y
  

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  stage2_coef <- coef(stage2_fit)[-1]
  
  theta_to_beta_v2(stage2_coef, selected_vars, p)
}

two_step_yhat <- function(x, y) {
  #two step procedure for gaussian case
  #step 1: constrained lasso
  #step 2: sparsenet

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- z %*% constr_beta
  

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  stage2_coef <- coef(stage2_fit)[-1]
  
  theta_to_beta_v2(stage2_coef, selected_vars, p)
}

two_step2 <- function(x, y) {
  #two step procedure for gaussian case
  #step 1: constrained lasso
  #step 2: sparsenet

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- y
  

  stage2_fit <- sparsenet(expanded_set, resp)
  stage2_coef <- coef(stage2_fit)$g9[-1,15]  

  theta_to_beta_v2(stage2_coef, selected_vars, p)
}


two_step4 <- function(x, y) {
  #fit L_0 by cv

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- y
  

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  best_lam <- max(which.min(stage2_fit$cvm[,9]) - 5, 1)
  stage2_coef <- coef(stage2_fit$sparsenet.fit)$g9[-1, best_lam]  

  theta_to_beta_v2(stage2_coef, selected_vars, p)
}

two_step5 <- function(x, y) {
  #fit L_0 by cv

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- z %*% constr_beta

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  stage2_coef <- coef(stage2_fit$sparsenet.fit)$g1[-1, 50]  

  theta_to_beta_v2(stage2_coef, selected_vars, p)
}

two_step6 <- function(x, y) {
  #fit L_0 by cv

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- y

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  stage2_coef <- coef(stage2_fit, which = "parms.1se") 

  theta_to_beta_v2(stage2_coef, selected_vars, p)
}

post_proc <- function(x, y) {
  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  if(length(selected_vars) == 0) {return(rep(0,p))}
  expanded_set <- small_to_big_z(z, selected_vars)

  resp <- z %*% constr_beta
  d <- sum(constr_beta) - 1
  if(d <= 0) {
    return(constr_beta)
  }

  subs <- regsubsets(expanded_set, resp, method = "forward", nvmax = d)
  r2 <- subs$rss / subs$nullrss
  cut <- which.max(r2 < .05)
  if(cut == 1) {return(rep(0,p))}
  stage2_coef <- rep(0, ncol(expanded_set))
  vars <- as.numeric(names(coef(subs, cut)[-1]))
  if(sum(is.na(vars)) > 0) {return(rep(0,p))}
  stage2_coef[vars] <- coef(subs, cut)[-1]
  
  theta_to_beta_v2(stage2_coef, selected_vars, p)
}


small_to_big_z <- function(z, indices) {
  d <- length(indices)
  if(d <= 1) {
    return(NULL)
  }

  big_z <- matrix(0, nrow = nrow(z), ncol = d*(d - 1) / 2)
  
  col <- 1
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      big_z[, col] <- z[, indices[i]] - z[, indices[j]]
      col <- col + 1
    }
  }
  
  big_z
}

spanning_basis <- function(z, indices) {
  d <- length(indices)
  big_z <- matrix(0, nrow = nrow(z), ncol = d-1)
  
  col <- 1
  for(j in 2:d) {
    big_z[, col] <- z[, indices[1]] - z[, indices[j]]
    col <- col + 1
  }
  
  big_z
}

theta_to_beta_v2 <- function(theta, indices, p) {
  beta <- rep(0,p)
  indices <- which(indices)
  d <- length(indices)
  
  k <- 1
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      beta[indices[i]] <- beta[indices[i]] + theta[k]
      beta[indices[j]] <- beta[indices[j]] - theta[k]
      k <- k + 1
    }
  }
  
  beta
}

two_step_cvm <- function(x, y) {
  #two step procedure for gaussian case
  #step 1: constrained lasso
  #step 2: sparsenet

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  resp <- y
  

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  
  stage2_fit$cvm[stage2_fit$which.min[1], stage2_fit$which.min[2]]
}


two_step_debiased <- function(x, y) {
  #two step procedure for gaussian case
  #step 1: constrained lasso
  #step 2: sparsenet

  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian")
  cv_fit <- cv.glmnet.constr(constrained_fit, z, y)
  best <- which.min(cv_fit$cvm)
  constr_beta <- constrained_fit$b[,best]
  
  selected_vars <- which(abs(constr_beta) > 0)
  expanded_set <- small_to_big_z(z, selected_vars)
  
  basis <- spanning_basis(z, selected_vars)
  model <- lm(y ~ basis - 1)
  y_hat_debiased <- predict(model)

  resp <- y_hat_debiased
  

  stage2_fit <- cv.sparsenet(expanded_set, resp)
  stage2_coef <- coef(stage2_fit)[-1]
  
  theta_to_beta_v2(stage2_coef, selected_vars, p)
}

two_step_generic <- function(x, y, lambda_1 = NULL, k_max = NULL) {
  p <- ncol(x)
  
  y <- y - mean(y)
  z <- scale(log(x), center = TRUE, scale = FALSE)

  constrained_fit <- glmnet.constr(z, y, family = "gaussian", lambda = lambda_1)
  betas <- constrained_fit$beta
  selected_vars <- apply(betas, 2, function(x){(abs(x) > 0)})

  output <- list
  for (i in 1:length(lambda_1)) {
    expanded_set <- small_to_big_z(z, selected_vars[,i])
    y_hat <- z %*% betas[,i]
    resp <- y  

    d <- sum(selected_vars) - 1
    k <- min(d, k_max)

    output[[i]] <- custom_fs(expanded_set, resp, k_max, selected_vars, p)
  }
}

custom_fs <- function(expanded_set, resp, k_max, selected_vars, p) {
  if(k_max <= 0) {
    return(0)
  }
  subs <- myfs(expanded_set, resp, nsteps = k_max, center = FALSE)

  list(theta_vals = subs$bhat, theta_ind = subs$ind, subset = which(selected_vars))
}

theta_long_to_wide <- function(theta) {
  p <- (1 + sqrt(1 + 8*length(theta))) / 2
  theta_mat <- matrix(0, p, p)

  k <- 1
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      theta_mat[i,j] <- theta[k]
      k <- k + 1
    }
  }

  theta_mat
}

theta_sel_to_theta <- function(theta, vars, indices, p) {
  beta <- rep(0,p)
  d <- length(indices)
  theta_out <- rep(0, p * (p-1) / 2)

  k <- 1
  k_inner <- 1
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      if(i %in% indices[vars] & j %in% indices[vars]) {
        theta_out[k_inner] <- theta[k]
        k <- k + 1
        print(i)

      }
      k_inner <- k_inner + 1
    }
  }
  
  theta_out
}
# theta_sel_to_theta(c(1), c(2,3), 3)
# theta_sel_to_theta(c(1), c(1,3), 3)
# theta_sel_to_theta(c(1,2,3), c(1,2,3), 3)

theta_wide_to_long <- function(theta_mat) {
  p <- ncol(theta_mat)

  k <- 1
  theta <- rep(0, p * (p - 1) / 2)
  for(i in 1:(p-1)) {
    for(j in (i+1):p) {
      theta[k] <- theta_mat[i, j]
      k <- k + 1
    }
  }

  theta
}

autoRegressiveCorr <- function(dim, corr) {
  #returns autoregressive correlation matrix with specified correlation
  mat <- matrix(0, nrow = dim, ncol = dim)
  for(row in 1:dim) {
    for(col in 1:dim) {
      mat[row, col] <- corr ** abs(row - col)
    }
  }
  
  mat
}
