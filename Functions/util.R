
##### help generate data #######################################################

sigmoid = function(x) 1 / (1 + exp(-x))

alrinv <- function(y) {
  x <- cbind(exp(y), 1)
  x / rowSums(x)
}

geometric.mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

rgExpDecay <- function(
    p,  # number of variables
    rho, # the decaying correlation
    const = 0.01
){
  Sigma = matrix(rho, p, p)
  index = matrix(rep(1:p, p), p, p) - t(matrix(rep(1:p, p), p, p))
  Sigma = Sigma^(abs(index))
  Omega = solve(Sigma)
  Omega <- (Omega + t(Omega))/2
  Adj <- (abs(Omega)>1e-08) - diag(rep(1,p))
  # weights <- matrix(rho, p, p)
  # Omega <- Adj*weights
  # diag(Omega) <- ifelse(min(eigen(Omega)$val) > 0, 0, - min(eigen(Omega)$val)) + const
  # Sigma <- solve(Omega)
  # Sigma <- cov2cor(Sigma)
  # Omega <- solve(Sigma)
  return(list(Sigma = Sigma, Omega = Omega, Adj = Adj))
}


##### preprocess data for selbal ###############################################

getSelbalData = function(X = X, y = y, classification = FALSE, levels, labels){
  if(is.null(rownames(X)) || is.null(colnames(X))){
    X.slbl = X
    if(is.null(rownames(X))){
      rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
    }
    if(is.null(colnames(X))){
      colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
    }
  }
  if(classification && !is.factor(y)){
    if(is.null(levels) & is.null(labels)){
      stop("getSelbalData(): classification = TRUE, but levels and labels were not provided. Will let factor() choose them.")
    }
    y.slbl = factor(y, levels = levels, labels = labels)
  } else{
    y.slbl = y
  }
  return(list(X = X.slbl, y = y.slbl))
}



##### extract balance regression coefficients ##################################

getCoefsBM = function(coefs, sbp){
  # balance model defined in manuscript, i.e., without normalizing balance
  if(names(coefs)[1] == "(Intercept)"){ # model with intercept
    a0 = coefs[1]
    if(length(coefs) == 1){
      bm.coefs = 0
    } else{
      bm.coefs = coefs[-1]
    }
  } else{ # model without intercept
    a0 = NA
    bm.coefs = coefs
  }
  llc.coefs = getGammaFromTheta(bm.coefs, sbp)
  return(list(
    a0 = a0, 
    bm.coefs = bm.coefs,
    llc.coefs = llc.coefs
  ))
}

getSBPSelbal = function(X, selbal.fit){
  p = ncol(X)
  
  # U (transformation) matrix
  sbp = rep(0, p)
  names(sbp) = colnames(X)
  pba.pos = unlist(subset(
    selbal.fit$global.balance, subset = Group == "NUM", select = Taxa))
  num.pos = length(pba.pos)
  pba.neg = unlist(subset(
    selbal.fit$global.balance, subset = Group == "DEN", select = Taxa))
  num.neg = length(pba.neg)
  sbp[pba.pos] = 1
  sbp[pba.neg] = -1
  sbp = matrix(sbp)
  rownames(sbp) = colnames(X)
  return(sbp)
}

getCoefsSelbal = function(
    X, y, selbal.fit, classification = FALSE, check = TRUE, covar = NULL){
  p = ncol(X)
  
  # U (transformation) matrix
  sbp = getSBPSelbal(X = X, selbal.fit = selbal.fit)
  ilrU = getIlrTrans(sbp = sbp)
  if(check){
    if(is.factor(y)) y = as.numeric(y) - 1
    data = data.frame(cbind(y, covar, log(as.matrix(X)) %*% matrix(ilrU)))
    colnames(data)[ncol(data)] <- "V1"
    if(!classification){
      supposed.fit = lm(y ~ ., data = data)
    } else{
      supposed.fit = glm(y ~ ., data = data, family = binomial(link = "logit"))
    }
    if(!isTRUE(all.equal(as.numeric(coefficients(selbal.fit$glm)), 
                         as.numeric(coefficients(supposed.fit))))){
      warning("getCoefsSelbal(): coefficients of given model and supposed model don't match!")
    }
  }
  a0 = coefficients(selbal.fit$glm)["(Intecept)"]
  bm.coefs = coefficients(selbal.fit$glm)["V1"]
  llc.coefs = ilrU %*% as.matrix(bm.coefs)
  rownames(llc.coefs) = colnames(X)
  return(list(
    a0 = a0, 
    bm.coefs = bm.coefs,
    llc.coefs = llc.coefs, 
    sbp = sbp
  ))
}



##### transform coefficients for linear log contrasts model ####################

# getBetaFromTheta = function(
#     theta, sbp = NULL
# ){
#   # get U
#   U = getIlrTrans(sbp = sbp)
#   # calculate beta
#   beta = U %*% theta
#   return(beta)
# }

# getBetaFromCodacore = function(SBP_codacore, coeffs_codacore, p){
#   if(!is.matrix(SBP_codacore)) SBP_codacore = matrix(SBP_codacore)
#   if(ncol(SBP_codacore) != length(coeffs_codacore)){
#     stop("getBetaFromCodacore: SBP and coefficients don't match")
#   }
#   kplus = apply(SBP_codacore, 2, function(col) sum(col == 1))
#   kminus = apply(SBP_codacore, 2, function(col) sum(col == -1))
#   reciprocals = matrix(
#     0, nrow = nrow(SBP_codacore), ncol = ncol(SBP_codacore))
#   for(i in 1:ncol(SBP_codacore)){
#     reciprocals[SBP_codacore[, i] == 1, i] = 1 / kplus[i]
#     reciprocals[SBP_codacore[, i] == -1, i] = -1 / kminus[i]
#   }
#   ReciprocalstimesCoeffs = matrix(
#     NA, nrow = nrow(SBP_codacore), ncol = ncol(SBP_codacore))
#   for(i in 1:ncol(ReciprocalstimesCoeffs)){
#     ReciprocalstimesCoeffs[, i] = reciprocals[, i] * coeffs_codacore[i]
#   }
#   beta = rowSums(ReciprocalstimesCoeffs)
#   return(beta)
# }

getGammaFromTheta = function(theta, sbp){
  # theta defined in manuscript, i.e., coefficient of non-normalized balance
  if(!is.matrix(sbp)) sbp = matrix(sbp)
  if(ncol(sbp) != length(theta)){
    stop("getGammaFromTheta: SBP and coefficients don't match")
  }
  kplus = apply(sbp, 2, function(col) sum(col == 1))
  kminus = apply(sbp, 2, function(col) sum(col == -1))
  reciprocals = matrix(
    0, nrow = nrow(sbp), ncol = ncol(sbp))
  for(i in 1:ncol(sbp)){
    reciprocals[sbp[, i] == 1, i] = 1 / kplus[i]
    reciprocals[sbp[, i] == -1, i] = -1 / kminus[i]
  }
  ReciprocalstimesCoeffs = matrix(
    NA, nrow = nrow(sbp), ncol = ncol(sbp))
  for(i in 1:ncol(ReciprocalstimesCoeffs)){
    ReciprocalstimesCoeffs[, i] = reciprocals[, i] * theta[i]
  }
  beta = rowSums(ReciprocalstimesCoeffs)
  names(beta) = rownames(sbp)
  return(beta)
}

# slr.fromSBP <- function(x, y){
#   
#   if(!identical(colnames(x), rownames(y))){
#     
#     stop("Component names for data matrix and balance matrix do not match.")
#   }
#   
#   x <- as.matrix(x)
#   
#   if(any(x == 0)){
#     
#     message("Alert: Replacing 0s with next smallest value to calculate balances.")
#     zeros <- x == 0
#     x[zeros] <- min(x[!zeros])
#   }
#   
#   res <- apply(y, 2, function(z) slr.fromContrast(x, z))
#   rownames(res) <- as.character(1:nrow(res))
#   return(res)
# }
slr.fromContrast <- function(x, contrast){
  
  if(length(contrast) != ncol(x)) stop("Contrast must have length ncol(x) = D.")
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")
  
  # lpos <- sum(contrast == 1)
  # lneg <- sum(contrast == -1)
  # const <- sqrt((lpos*lneg)/(lpos+lneg))
  logX <- log(x)
  ipos <- rowMeans(logX[, contrast == 1, drop = FALSE])
  ineg <- rowMeans(logX[, contrast == -1, drop = FALSE])
  
  ipos - ineg
}

##### balance regression model tools ###########################################

# Compute the transformation vector/matrix from SBP vector/matrix
getIlrTrans = function(sbp = NULL, detailed = FALSE){
  # check that sbp is in matrix form, even if it's just a p-vector
  # normalize SBP
  if(!all(sbp == 0)){
    pos.vec = colSums(sbp == 1)
    neg.vec = colSums(sbp == -1)
    U = sbp
    for(i in 1:ncol(U)){
      # divide 1s by pos.vec
      U[(U[, i] == 1), i] = 1 / pos.vec[i]
      # divide -1s by neg.vec
      U[(U[, i] == -1), i] = -1 / neg.vec[i]
    }
    U.unscaled = U
    # multiply by sqrt(rs / (r + s))
    norm.const = sqrt((pos.vec * neg.vec) / (pos.vec + neg.vec))
    U = sweep(U, MARGIN = 2, STATS = norm.const, FUN = "*")
  } else{
    U = sbp
    norm.const = NA
    pos.vec = 0
    neg.vec = 0
  }
  if(detailed){
    return(
      list(
        ilr.trans = U, 
        ilr.trans.unscaled = U.unscaled,
        const = norm.const, 
        pos.vec = pos.vec, 
        neg.vec = neg.vec
      )
    )
  } else{
    return(U)
  }
}

# getIlrX = function(X, sbp = NULL, U = NULL){
#   p = dim(X)[2]
#   
#   # checks
#   if(!("matrix" %in% class(X))){
#     if("numeric" %in% class(X)){
#       X = matrix(X, ncol = length(X))
#     } else{
#       warning("getIlrX: X is neither of class matrix nor numeric!!")
#     }
#   }
#   if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
#   
#   # get U
#   if(is.null(U)){
#     U = getIlrTrans(sbp = sbp)
#   }
#   
#   # calculate balances from U
#   ilr_balances = log(X) %*% U
#   
#   return(ilr_balances)
# }



##### metrics ##################################################################

randidx = function(sbp1, sbp2, adjusted = FALSE){
  sbp1 = as.numeric(sbp1 + 2)
  sbp2 = as.numeric(sbp2 + 2)
  if(adjusted){
    fossil::adj.rand.index(sbp1, sbp2)
  } else{
    fossil::rand.index(sbp1, sbp2)
  }
}

# getMSEyhat = function(
#     y, n, betahat0, betahat, predictors){
#     yhat = betahat0 + predictors %*% betahat
#     mse = as.vector(crossprod(y - yhat) / n)
#   return(mse)
# }

getEstimationAccuracy = function(
    true.llc.coefs, est.llc.coefs, relative = FALSE
  ){
  # L1 distance
  EA1 = sum(abs(est.llc.coefs - true.llc.coefs))
  if(relative){
    EA1 = EA1 / sum(abs(true.llc.coefs))
  }
  # L2 distance
  EA2 = as.vector(sqrt(crossprod(est.llc.coefs - true.llc.coefs)))
  if(relative){
    EA2 = EA2 / sqrt(crossprod(true.llc.coefs))
  }
  # LInfty distance
  EAInfty = max(abs(est.llc.coefs - true.llc.coefs))
  if(relative){
    EAInfty = EAInfty / max(abs(true.llc.coefs))
  }
  return(list(
    EA1 = EA1, 
    EA2 = EA2, 
    EAInfty = EAInfty
  ))
}
getSelectionAccuracy = function(non0.true.llc.coefs, non0.est.llc.coefs){
  is0.true.llc.coefs = !non0.true.llc.coefs
  is0.est.llc.coefs = !non0.est.llc.coefs
  FP = sum(non0.est.llc.coefs & is0.true.llc.coefs)
  FN = sum(is0.est.llc.coefs & non0.true.llc.coefs)
  TP = sum(non0.est.llc.coefs & non0.true.llc.coefs)
  TPR = TP / sum(non0.true.llc.coefs)
  FPR = FP / sum(is0.true.llc.coefs)
  # F-score = precision / recall
  # precision = # true positive results / # of positive results
  #   (including those not identified correctly)
  precision = TP / sum(non0.est.llc.coefs) 
  # recall = sensitivity = TPR = # true positive results / # of true positives
  Fscore = 2 * precision * TPR / (precision + TPR)
  return(list(
    TPR = TPR, 
    FPR = FPR,
    Fscore = Fscore
  ))
}
# package metrics
getMetricsLLC = function(
    # y.train = NULL, y.test = NULL, logX.train = NULL, logX.test = NULL, 
    # n.train = NULL, n.test = NULL,
    # betahat0 = NULL, 
    est.llc.coefs, 
    true.sbp = NULL, true.llc.coefs = NULL, non0.true.llc.coefs,
    metrics = c("estimation", "selection"), 
    relative = FALSE
){
  result = list()
  
  # if("prediction" %in% metrics){
  #   if(is.null(y.train) | is.null(y.test) | 
  #      is.null(logX.train) | is.null(logX.test) | 
  #      is.null(n.train) | is.null(n.test) | is.null(betahat0)) {
  #     stop("getMetricsLLC() prediction: one or more of these are missing -- y.train, y.test, logX.train, logX.test, n.train, n.test, betahat0")
  #   }
  #   # 1. prediction error #
  #   # 1a. on training set #
  #   result$PEtr = getMSEyhat(
  #     y = y.train, n = n.train, betahat0 = betahat0, betahat = betahat, 
  #     predictors = logX.train)
  #   # 1b. on test set #
  #   result$PEte = getMSEyhat(
  #     y = y.test, n = n.test, betahat0 = betahat0, betahat = betahat, 
  #     predictors = logX.test)
  # }
  
  if("estimation" %in% metrics){
    if(is.null(true.llc.coefs)) {
      stop("getMetricsLLC() estimation: true.llc.coefs arg is missing")
    }
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    EA = getEstimationAccuracy(true.llc.coefs = true.llc.coefs, est.llc.coefs = est.llc.coefs, relative = relative)
    # save
    result$EA1 = EA$EA1
    result$EA2 = EA$EA2
    result$EAInfty = EA$EAInfty
  }
  
  if("selection" %in% metrics){
    if(is.null(true.llc.coefs) & is.null(true.sbp)) {
      stop("getMetricsLLC() selection: true.llc.coefs and true.sbp args are both missing, need at least one of them")
    }
    if(is.null(true.llc.coefs) & !is.null(true.sbp)){
      true.llc.coefs = as.numeric(apply(true.sbp, 1, function(row) any(row != 0)))
    }
    
    # 3. selection accuracy #
    non0.est.llc.coefs = abs(est.llc.coefs) > 10e-8
    SA = getSelectionAccuracy(
      non0.true.llc.coefs = non0.true.llc.coefs, 
      non0.est.llc.coefs = non0.est.llc.coefs)
    # save
    result$TPR = SA$TPR
    result$FPR = SA$FPR
    result$Fscore = SA$Fscore
  }
  return(unlist(result))
}
getMetricsBM = function(
    # y.train = NULL, y.test = NULL, ilrX.train = NULL, ilrX.test = NULL, 
    # n.train = NULL, n.test = NULL,
    # thetahat0 = NULL, thetahat, 
    est.llc.coefs, true.sbp = NULL,
    true.llc.coefs = NULL, non0.true.llc.coefs,
    metrics = c("estimation", "selection"), 
    relative = FALSE
){
  result = list()
  
  # if("prediction" %in% metrics){
  #   if(is.null(y.train) | is.null(y.test) | 
  #      is.null(ilrX.train) | is.null(ilrX.test) | 
  #      is.null(n.train) | is.null(n.test) | is.null(thetahat0)) {
  #     stop("getMetricsBM() prediction: one or more of these are missing -- y.train, y.test, ilrX.train, ilrX.test, n.train, n.test, thetahat0")
  #   }
  #   # 1. prediction error #
  #   # 1a. on training set #
  #   result$PEtr = getMSEyhat(
  #     y = y.train, n = n.train, betahat0 = thetahat0, betahat = thetahat,
  #     predictors = ilrX.train)
  #   # 1b. on test set #
  #   result$PEte = getMSEyhat(
  #     y = y.test, n = n.test, betahat0 = thetahat0, betahat = thetahat,
  #     predictors = ilrX.test)
  # }
  
  if("estimation" %in% metrics){
    if(is.null(true.llc.coefs)) {
      stop("getMetricsBM() estimation: true.llc.coefs arg is missing")
    }
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    EA = getEstimationAccuracy(true.llc.coefs = true.llc.coefs, est.llc.coefs = est.llc.coefs, relative = relative)
    # save
    result$EA1 = EA$EA1
    result$EA2 = EA$EA2
    result$EAInfty = EA$EAInfty
  }
  
  if("selection" %in% metrics){
    if(is.null(true.llc.coefs) & is.null(true.sbp)) {
      stop("getMetricsBM() selection: true.llc.coefs and true.sbp args are both missing, need at least one of them")
    }
    if(is.null(true.llc.coefs) & !is.null(true.sbp)){
      true.llc.coefs = as.numeric(apply(true.sbp, 1, function(row) any(row != 0)))
    }
    # 3. selection accuracy #
    non0.est.llc.coefs = abs(est.llc.coefs) > 10e-8
    SA = getSelectionAccuracy(
      non0.true.llc.coefs = non0.true.llc.coefs,
      non0.est.llc.coefs = non0.est.llc.coefs)
    # save
    result$TPR = SA$TPR
    result$FPR = SA$FPR
    result$Fscore = SA$Fscore
  }
  return(unlist(result))
}


# population correlation matrix ################################################
popGammajk = function(
    alpha1j, alpha1k, beta1, var_epsilon, var_epsilonj, var_epsilonk, U){
  varU = stats::var(U) #(1 / 12) * (0.5 - (-0.5))
  corrjk = ((alpha1j - alpha1k) * beta1 * varU) / 
    sqrt((beta1^2 * varU + var_epsilon) * 
           ((alpha1j - alpha1k)^2 * varU + var_epsilonj + var_epsilonk))
  return(abs(corrjk))
}
popGamma = function(
    alpha1, beta1, var_epsilon, var_epsilon2, U
){
  p = length(alpha1)
  if(length(var_epsilon2) == 1) var_epsilon2 = rep(var_epsilon2, p)
  
  rhoMat = matrix(0, p, p)
  for (j in 1:p){
    for (k in 1:p){
      if (k==j){next}
      else {
        rhoMat[j, k] = popGammajk(
          alpha1j = alpha1[j], alpha1k = alpha1[k], beta1 = beta1, 
          var_epsilon = var_epsilon, var_epsilonj = var_epsilon2[j], 
          var_epsilonk = var_epsilon2[k], U = U)
      }
    }
  }
  return(rhoMat)
}
