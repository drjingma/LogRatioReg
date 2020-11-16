# Calculates objective function of constrained lasso for compositional data
#   with subcompositional coherence,
#   given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
compositionalLasso = function(Xtilde, Ytilde, mu, beta, gamma, lambda){
  n = dim(Xtilde)[1]
  objective = sum(crossprod(Ytilde - Xtilde %*% beta)) / (2 * n) + 
    lambda * sum(abs(beta)) + # typical Lasso objective up until this point
    gamma * sum(beta) + mu * sum(beta)^2 / 2
  return(objective)
}

# Fit LASSO on (standardized?) data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta_start - p vector, an optiona starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.0001
fitCompositionalLASSOstandardized = function(
  Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001, mu = 1){
  
  n = dim(Xtilde)[1]
  p = dim(Xtilde)[2]
  
  # Starting value for alpha = gamma / mu
  alpha = 0
  
  # Check that n is the same between Xtilde and Ytilde
  if(length(Ytilde) != n) stop("Y does not have same number of observations as X (n observations).")
  
  # Check that lambda is non-negative
  if(lambda < 0) stop("lambda is negative, when it should be non-negative.")
  
  # Check for starting point beta_start. If none supplied, initialize with a vector of zeros. If supplied, check for compatibility with Xtilde in terms of p
  if(is.null(beta_start)){
    beta_start = rep(0, p)
  } else{
    # check for compatibility
    if(length(beta_start) != p) stop("beta_start does not have same number of parameters as X (p parameters).")
  }
  
  # Coordinate-descent implementation. Stop when the difference between objective functions is less than eps.
  old_objective = compositionalLasso(Xtilde, Ytilde, mu, beta_start, gamma = alpha * mu, lambda)
  objectives = c(old_objective) # for testing ##################################
  difference = 2 * eps
  beta = beta_start
  difference_beta = 2 * eps
  R = Ytilde - Xtilde %*% matrix(beta)
  while(difference_beta >= eps){
    old_beta = beta
    for (j in 1:p){
      old_betaj = beta[j]
      A_j = (1 / n) * Xtilde[ , j] %*% R + old_betaj # unconstrained LASSO softthresh arg, Y - X_{-j} beta_{-j}^k, partial residual
      A2_j = A_j - mu * (sum(beta[-j]) + alpha) # constrained LASSO softthresh arg
      constS = 1 / ((crossprod(Xtilde[, j]) / n) + mu) # constrained LASSO constant in front of softthresh
      beta[j] = constS * soft(A2_j, lambda)
      R = R + Xtilde[ , j] * (old_betaj - beta[j]) # update residual
    }
    objective = compositionalLasso(Xtilde, Ytilde, mu, beta, gamma = alpha * mu, lambda)
    objectives = c(objectives, objective) # for testing ########################
    difference = sqrt((old_objective - objective)^2) # should be positive, but just in case.
    old_objective = objective
    # update alpha
    alpha = alpha + sum(beta)
    difference_beta = sqrt(sum((old_beta - beta)^2))
    # print(paste0("difference in beta: ", difference_beta))
    # print(paste0("difference in loss: ", difference))
  }
  fmin = old_objective
  return(list(beta = beta, fmin = fmin, objectives = objectives))
}


# Fit LASSO on standardized data for a sequence of lambda values. Sequence version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.0001
fitCompositionalLASSOstandardized_seq = function(
  Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 50, eps = 0.0001, mu = 1){
  
  # Check that n is the same between Xtilde and Ytilde
  n = dim(Xtilde)[1]
  p = dim(Xtilde)[2]
  
  # Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >=0, and make sure the values are sorted from largest to smallest. If none of the supplied values satisfy the requirement, print the warning message and proceed as if the values were not supplied.
  # If lambda_seq is not supplied, calculate lambda_max (the minimal value of lambda that gives zero solution), and create a sequence of length n_lambda as
  lambda_seq = sort(lambda_seq[which(lambda_seq >= 0)], decreasing = TRUE)
  
  if(is.null(lambda_seq) | length(lambda_seq) == 0){
    lambda_max = max(abs(crossprod(Xtilde, Ytilde) / n))
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  } else{
    n_lambda = length(lambda_seq)
  }
  
  # Apply fitLASSOstandardized going from largest to smallest lambda (make sure supplied eps is carried over).
  if(length(lambda_seq) != n_lambda) stop("lambda_seq does not have length equal to n_lambda")
  
  fitLASSOstd_init = fitCompositionalLASSOstandardized(
    Xtilde, Ytilde, lambda_seq[1], beta_start = NULL, eps, mu)
  beta_mat = matrix(rep(fitLASSOstd_init$beta, n_lambda), p, n_lambda)
  fmin_vec = rep(fitLASSOstd_init$fmin, n_lambda)
  # Use warm starts strategy for setting the starting values.
  if(n_lambda > 1){
    for(i in 2:n_lambda){
      fitLASSOstd = fitCompositionalLASSOstandardized(
        Xtilde, Ytilde, lambda_seq[i], beta_start = beta_mat[ , (i - 1)], eps, mu)
      beta_mat[ , i] = fitLASSOstd$beta
      fmin_vec[i] = fitLASSOstd$fmin
    }
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.0001
fitCompositionalLASSO = function(
  X ,Y, lambda_seq = NULL, n_lambda = 50, eps = 0.0001, mu = 1){
  n = dim(X)[1]
  p = dim(X)[2]
  if(!is.null(lambda_seq)) n_lambda = length(lambda_seq)
  
  # Center and standardize X,Y based on standardizeXY function
  stdXY = standardizeXY(X, Y)
  Xtilde = stdXY$Xtilde
  Ytilde = stdXY$Ytilde
  Ymean = stdXY$Ymean
  Xmeans = stdXY$Xmeans
  weights = stdXY$weights
  
  # Fit Lasso on a sequence of values using fitLASSOstandardized_seq (make sure the parameters carry over)
  fitLASSO_seq = fitCompositionalLASSOstandardized_seq(
    Xtilde, Ytilde, lambda_seq, n_lambda, eps, mu)
  lambda_seq = fitLASSO_seq$lambda_seq
  
  # Perform back scaling and centering to get original intercept and coefficient vector for each lambda
  beta_mat = matrix(rep(NA, p * n_lambda), p, n_lambda)
  beta0_vec = rep(NA, n_lambda)
  for(i in 1:n_lambda){
    tildebeta = fitLASSO_seq$beta_mat[ , i]
    beta_mat[ , i] = diag(weights) %*% tildebeta
    beta0_vec[i] = Ymean - crossprod(Xmeans, beta_mat[ , i])
  }
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# eps - precision level for convergence assessment, default 0.0001
cvCompositionalLASSO = function(
  X ,Y, lambda_seq = NULL, n_lambda = 50, k = 10, eps = 0.0001, mu = 1){
  n = dim(X)[1]
  if(!is.null(lambda_seq)) n_lambda = length(lambda_seq)
  
  # Fit Lasso on original data using fitLASSO
  Lasso = fitCompositionalLASSO(
    X ,Y, lambda_seq, n_lambda, eps, mu)
  beta_mat = Lasso$beta_mat
  beta0_vec = Lasso$beta0_vec
  lambda_seq = Lasso$lambda_seq
  
  # Split the data into K folds
  shuffle = sample(1:n)
  idfold = (shuffle %% k) + 1
  n_fold = as.vector(table(idfold))
  
  cvm = rep(NA, n_lambda) # want to have CV(lambda)
  cvm_sqerror = matrix(rep(NA, k * n_lambda), k, n_lambda) # calculate squared error for each fold, needed for CV(lambda) calculation
  cvse = rep(NA, n_lambda) # want to have SE_CV(lambda)
  cvse_errorsd =  matrix(rep(NA, k * n_lambda), k, n_lambda) # calculate sd of error for each fold, needed for SE_CV(lambda) calculation
  # Calculate Lasso for each fold removed
  for (j in 1:k){
    # Training data
    Xtrain = X[idfold != j, ]
    Ytrain = Y[idfold != j]
    # Test data
    Xtest = X[idfold == j, ]
    Ytest = Y[idfold == j]
    
    # Calculate LASSO on that fold using fitLASSO
    Lasso_j = fitCompositionalLASSO(
      Xtrain, Ytrain, lambda_seq, n_lambda, eps, mu)
    
    # Any additional calculations that are needed for calculating CV and SE_CV(lambda)
    for(m in 1:n_lambda){
      Ypred = Lasso_j$beta0_vec[m] + Xtest %*% Lasso_j$beta_mat[ , m]
      cvm_sqerror[j, m] = sum(crossprod(Ytest - Ypred))
      cvse_errorsd[j, m] = cvm_sqerror[j, m] / n_fold[j]
    }
    
    
  }
  
  # Calculate CV(lambda) and SE_CV(lambda) for each value of lambda
  cvse = sqrt(apply(cvse_errorsd, 2, var)) / sqrt(k)
  cvm = colMeans(cvm_sqerror)
  
  # Find lambda_min = argmin{CV(lambda)}
  lambda_min_index = which.min(cvm)
  lambda_min = lambda_seq[lambda_min_index]
  
  # Find lambda_1se = maximal lambda s.t. CV(lambda) <= CV(lambda_min) + CV_SE(lambda_min)
  bound = cvm[lambda_min_index] + cvse[lambda_min_index]
  lambda_1se_index = which(cvm <= bound)
  lambda_1se = lambda_seq[min(lambda_1se_index)]
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # Other output
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(
    list(
      lambda_seq = lambda_seq, 
      beta_mat = beta_mat, 
      beta0_vec = beta0_vec, 
      lambda_min = lambda_min, 
      lambda_1se = lambda_1se, 
      cvm = cvm, 
      cvm_idx = lambda_min_index,
      cvse = cvse,
      cvse_idx = lambda_1se_index
    ))
}