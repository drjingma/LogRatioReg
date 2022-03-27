

sigmoid = function(x) 1 / (1 + exp(-x))

alrinv <- function(y) {
  x <- cbind(exp(y), 1)
  x / rowSums(x)
}



##### extract coefficients #####################################################

getCoefsBM = function(coefs, sbp){
  if(names(coefs)[1] == "(Intercept)"){
    a0 = coefs[1]
    bm.coefs = coefs[-1]
  } else{
    a0 = NA
    bm.coefs = coefs
  }
  llc.coefs = getBetaFromTheta(bm.coefs, sbp)
  return(list(
    a0 = a0, 
    bm.coefs = bm.coefs,
    llc.coefs = llc.coefs
  ))
}

getCoefsSelbal = function(
  X, y, selbal.fit, classification = FALSE, check = TRUE){
  p = ncol(X)
  
  # U (transformation) matrix
  ilrU = rep(0, p)
  names(ilrU) = colnames(X)
  pba.pos = unlist(subset(
    selbal.fit$global.balance, subset = Group == "NUM", select = Taxa))
  num.pos = length(pba.pos)
  pba.neg = unlist(subset(
    selbal.fit$global.balance, subset = Group == "DEN", select = Taxa))
  num.neg = length(pba.neg)
  ilrU[pba.pos] = 1 / num.pos
  ilrU[pba.neg] = -1 / num.neg
  norm.const = sqrt((num.pos * num.neg) / (num.pos + num.neg))
  ilrU = norm.const * ilrU
  if(check){
    if(!classification){
      supposed.fit = lm(y ~ log(as.matrix(X)) %*% as.matrix(ilrU))
    } else{
      supposed.fit = glm(
        y ~ log(as.matrix(X)) %*% matrix(ilrU), family = binomial(link = "logit"))
    }
    if(!isTRUE(all.equal(as.numeric(coefficients(selbal.fit$glm)), 
                  as.numeric(coefficients(supposed.fit))))){
      warning("getCoefsSelbal(): coefficients of given model and supposed model don't match!")
    }
  }
  a0 = coefficients(selbal.fit$glm)[1]
  bm.coefs = coefficients(selbal.fit$glm)[-1]
  llc.coefs = ilrU %*% as.matrix(bm.coefs)
  rownames(llc.coefs) = colnames(X)
  return(list(
    a0 = a0, 
    bm.coefs = bm.coefs,
    llc.coefs = llc.coefs
  ))
}


##### transform coefficients for linear log contrasts model ####################

getBetaFromTheta = function(
  theta, sbp = NULL
){
  # get U
  U = getIlrTrans(sbp = sbp)
  # calculate beta
  beta = U %*% theta
  return(beta)
}

getBetaFromCodacore = function(SBP_codacore, coeffs_codacore, p){
  if(!is.matrix(SBP_codacore)) SBP_codacore = matrix(SBP_codacore)
  if(ncol(SBP_codacore) != length(coeffs_codacore)){
    stop("getBetaFromCodacore: SBP and coefficients don't match")
  }
  kplus = apply(SBP_codacore, 2, function(col) sum(col == 1))
  kminus = apply(SBP_codacore, 2, function(col) sum(col == -1))
  reciprocals = matrix(
    0, nrow = nrow(SBP_codacore), ncol = ncol(SBP_codacore))
  for(i in 1:ncol(SBP_codacore)){
    reciprocals[SBP_codacore[, i] == 1, i] = 1 / kplus[i]
    reciprocals[SBP_codacore[, i] == -1, i] = -1 / kminus[i]
  }
  ReciprocalstimesCoeffs = matrix(
    NA, nrow = nrow(SBP_codacore), ncol = ncol(SBP_codacore))
  for(i in 1:ncol(ReciprocalstimesCoeffs)){
    ReciprocalstimesCoeffs[, i] = reciprocals[, i] * coeffs_codacore[i]
  }
  beta = rowSums(ReciprocalstimesCoeffs)
  return(beta)
}



##### other selbal helper functions ############################################

getSelbalData = function(X = X, y = y, classification = TRUE, levels, labels){
  if(is.null(rownames(X)) && is.null(colnames(X))){
    X.slbl = X
    rownames(X.slbl) = paste("Sample", 1:nrow(X.slbl), sep = "_")
    colnames(X.slbl) = paste("V", 1:ncol(X.slbl), sep = "")
  }
  if(classificatin && !is.factor(y)){
    if(is.null(levels) | is.null(labels)){
      stop("getSelbalData(): classification = TRUE, but levels and labels were not provided. Will let factor() choose them.")
    }
    y.slbl = factor(y, levels = levels, labels = labels)
  } else{
    y.slbl = y
  }
  return(list(X = X.slbl, y = y.slbl))
}




##### balance regression model tools ###########################################

# Compute the transformation matrix U
#   using a hierarchical tree or SBP matrix
# old name: getU
getIlrTrans = function(sbp = NULL, detailed = FALSE){
  # normalize SBP matrix to get U
  pos.vec = colSums(sbp == 1)
  neg.vec = colSums(sbp == -1)
  U = sbp
  for(i in 1:ncol(U)){
    # divide 1s by r
    U[(U[, i] == 1), i] = 1 / pos.vec[i]
    # divide -1s by s
    U[(U[, i] == -1), i] = -1 / neg.vec[i]
  }
  # multiply by sqrt(rs / (r + s))
  norm.const = sqrt((pos.vec * neg.vec) / (pos.vec + neg.vec))
  U = sweep(U, MARGIN = 2, STATS = norm.const, FUN = "*")
  # U2 = t(t(U) * norm.const) # equal to U
  if(detailed){
    return(
      list(
        ilr.trans = U, 
        const = norm.const, 
        pos.vec = pos.vec, 
        neg.vec = neg.vec
      )
    )
  } else{
    return(U)
  }
}

getIlrX = function(X, sbp = NULL, U = NULL){
  p = dim(X)[2]
  
  # checks
  if(!("matrix" %in% class(X))){
    if("numeric" %in% class(X)){
      X = matrix(X, ncol = length(X))
    } else{
      warning("getIlrX: X is neither of class matrix nor numeric!!")
    }
  }
  if(is.null(colnames(X))) colnames(X) = paste("V", 1:p, sep = "")
  
  # get U
  if(is.null(U)){
    U = getIlrTrans(sbp = sbp)
  }
  
  # calculate balances from U
  ilr_balances = log(X) %*% U
  
  return(ilr_balances)
}

# Fit any balance regression model to compositional data X and response y
#   using the balances' associated binary tree, SBP matrix, or 
#   transformation matrix U
# old name: fitILR
fitBMLasso = function(
  y = NULL, X, sbp = NULL, lambda = NULL, nlam = 20, 
  intercept = TRUE, standardize = TRUE, classification = FALSE
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n){
    stop("fitBMLasso(): dimensions of y and X don't match!!")
  }
  if(is.null(colnames(X))){
    colnames(X) = paste("V", 1:p, sep = "")
  }
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
  # compute balances
  ilrX = getIlrX(X = X, sbp = sbp)
  
  # run glmnet
  if(!classification){
    glmnet.fit = glmnet(
      x = ilrX, y = y, lambda = lambda, nlambda = nlam, 
      intercept = intercept, standardize = standardize)
  } else{
    glmnet.fit = glmnet(
      x = ilrX, y = y, lambda = lambda, nlambda = nlam, 
      intercept = intercept, standardize = standardize, 
      family = binomial(link = "logit"))
  }
  
  # check lambda length
  if(nlam != length(glmnet.fit$lambda)){
    lambda <- log(glmnet.fit$lambda)
    lambda_new <- exp(seq(max(lambda), min(lambda),length.out = nlam))
    if(!classification){
      glmnet.fit = glmnet(
        x = ilrX, y = y, lambda = lambda_new, 
        intercept = intercept, standardize = standardize)
    } else{
      glmnet.fit = glmnet(
        x = ilrX, y = y, lambda = lambda_new, 
        intercept = intercept, standardize = standardize, 
        family = binomial(link = "logit"))
    }
  }
  return(list(
    theta0 = glmnet.fit$a0, 
    theta = glmnet.fit$beta, 
    lambda = glmnet.fit$lambda,
    glmnet = glmnet.fit, 
    sbp = sbp
  ))
}

# old name: cvILR
cvBMLasso = function(
  y = NULL, X, sbp = NULL, lambda = NULL, nlam = 20, 
  nfolds = 10, foldid = NULL, intercept = TRUE, standardize = TRUE, 
  keep = FALSE, classification = FALSE
){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n){
    stop("cvBMLasso(): dimensions of y and X don't match!!")
  }
  if(is.null(colnames(X))){
    colnames(X) = paste("V", 1:p, sep = "")
  }
  if(!is.null(lambda)){ # lambda is given
    nlam = length(lambda)
  }
  
  # compute balances
  ilrX = getIlrX(X = X, sbp = sbp)
  
  # run cv.glmnet
  if(!classification){
    cv_exact = cv.glmnet(
      x = ilrX, y = y, lambda = lambda, nlambda = nlam, nfolds = nfolds, 
      foldid = foldid, intercept = intercept, type.measure = c("mse"), 
      keep = keep, standardize = standardize)
  } else{
    cv_exact = cv.glmnet(
      x = ilrX, y = y, lambda = lambda, nlambda = nlam, nfolds = nfolds, 
      foldid = foldid, intercept = intercept, type.measure = c("mse"), 
      keep = keep, standardize = standardize, family = binomial(link = "logit"))
  }
  
  # check lambda length
  if(nlam != length(cv_exact$lambda)){
    lambda <- log(cv_exact$lambda)
    lambda_new <- exp(seq(max(lambda), min(lambda),length.out = nlam))
    if(!classification){
      cv_exact = cv.glmnet(
        x = ilrX, y = y, lambda = lambda_new, nfolds = nfolds, 
        foldid = foldid, intercept = intercept, type.measure = c("mse"), 
        keep = keep, standardize = standardize)
    } else{
      cv_exact = cv.glmnet(
        x = ilrX, y = y, lambda = lambda_new, nfolds = nfolds, 
        foldid = foldid, intercept = intercept, type.measure = c("mse"), 
        keep = keep, standardize = standardize, 
        family = binomial(link = "logit"))
    }
  }
  return(list(
    theta0 = cv_exact$glmnet.fit$a0,
    theta = cv_exact$glmnet.fit$beta,
    lambda = cv_exact$lambda,
    cv.glmnet = cv_exact,
    glmnet = cv_exact$glmnet.fit,
    cvm = cv_exact$cvm,
    cvsd = cv_exact$cvsd,
    sbp = sbp
  ))
}






##### metrics ##################################################################

getMSEyhat = function(y, n, betahat0, betahat, predMat, classification = FALSE){
  if(!classification){
    yhat = betahat0 + predMat %*% betahat
    mse = as.vector(crossprod(y - yhat) / n)
  } else{
    yhat = NULL
    mse = NULL
  }
  
  return(mse)
}
getEstimationAccuracy = function(true.beta, betahat){
  EA1 = sum(abs(betahat - true.beta))
  EA2 = as.vector(sqrt(crossprod(betahat - true.beta)))
  EAInfty = max(abs(betahat - true.beta))
  return(list(
    EA1 = EA1, 
    EA2 = EA2, 
    EAInfty = EAInfty
  ))
}
getSelectionAccuracy = function(is0.true.beta, non0.true.beta, non0.betahat){
  FP = sum(is0.true.beta & non0.betahat)
  FN = sum((non0.true.beta != non0.betahat) & non0.true.beta)
  TPR = sum((non0.true.beta == non0.betahat) & non0.betahat) / 
    sum(non0.true.beta)
  # F-score = precision / recall
  # precision = # true positive results / # of positive results
  #   (including those not identified correctly)
  precision = sum((non0.true.beta == non0.betahat) & non0.betahat) / 
    sum(non0.betahat) 
  # recall = sensitivity = TPR = # true positive results / # of true positives
  Fscore = 2 * precision * TPR / (precision + TPR)
  return(list(
    FP = FP, 
    FN = FN, 
    TPR = TPR, 
    precision = precision, 
    Fscore = Fscore
  ))
}
# package metrics
getMetricsLLC = function(
  y.train = NULL, y.test = NULL, logX.train = NULL, logX.test = NULL, 
  n.train = NULL, n.test = NULL,
  betahat0 = NULL, betahat, 
  true.sbp = NULL, true.beta = NULL, is0.true.beta, non0.true.beta,
  metrics = c("prediction", "betaestimation", "selection"), 
  classification = FALSE
){
  result = list()
  
  if("prediction" %in% metrics){
    if(is.null(y.train) | is.null(y.test) | 
       is.null(logX.train) | is.null(logX.test) | 
       is.null(n.train) | is.null(n.test) | is.null(betahat0)) {
      stop("getMetricsBalanceReg() prediction: one or more of these are missing -- y.train, y.test, logX.train, logX.test, n.train, n.test, betahat0")
    }
    # 1. prediction error #
    # 1a. on training set #
    result$PEtr = getMSEyhat(
      y = y.train, n = n.train, betahat0 = betahat0, betahat = betahat, 
      predMat = logX.train, classification = classification)
    # 1b. on test set #
    result$PEte = getMSEyhat(
      y = y.test, n = n.test, betahat0 = betahat0, betahat = betahat, 
      predMat = logX.test, classification = classification)
  }
  
  if("betaestimation" %in% metrics){
    if(is.null(true.beta)) {
      stop("getMetricsBalanceReg() betaestimation: true.beta arg is missing")
    }
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    EA = getEstimationAccuracy(true.beta = true.beta, betahat = betahat)
    # 2b. estimation accuracy for active set
    EA.active = getEstimationAccuracy(
      true.beta = true.beta[non0.true.beta], betahat = betahat[non0.true.beta])
    # 2c. estimation accuracy for inactive set
    EA.inactive = getEstimationAccuracy(
      true.beta = true.beta[is0.true.beta], betahat = betahat[is0.true.beta])
    # save
    result$EA1 = EA$EA1
    result$EA2 = EA$EA2
    result$EAInfty = EA$EAInfty
    result$EA1Active = EA.active$EA1
    result$EA2Active = EA.active$EA2
    result$EAInftyActive = EA.active$EAInfty
    result$EA1Inactive = EA.inactive$EA1
    result$EA2Inactive = EA.inactive$EA2
    result$EAInftyInactive = EA.inactive$EAInfty
  }
  
  if("selection" %in% metrics){
    if(is.null(true.beta) & is.null(true.sbp)) {
      stop("getMetricsBalanceReg() selection: true.beta and true.sbp args are both missing, need at least one of them")
    }
    if(is.null(true.beta) & !is.null(true.sbp)){
      true.beta = as.numeric(apply(true.sbp, 1, function(row) any(row != 0)))
    }
    
    
    # 3. selection accuracy #
    non0.betahat = abs(betahat) > 10e-8
    SA = getSelectionAccuracy(
      is0.true.beta = is0.true.beta, 
      non0.true.beta = non0.true.beta, 
      non0.betahat = non0.betahat)
    SA.pos = getSelectionAccuracy(
      is0.true.beta = is0.true.beta[true.beta > 0], 
      non0.true.beta = non0.true.beta[true.beta > 0], 
      non0.betahat = non0.betahat[true.beta > 0])
    SA.neg = getSelectionAccuracy(
      is0.true.beta = is0.true.beta[true.beta < 0], 
      non0.true.beta = non0.true.beta[true.beta < 0], 
      non0.betahat = non0.betahat[true.beta < 0])
    # save
    result$FP = SA$FP
    result$FN = SA$FN
    result$TPR = SA$TPR
    result$precision = SA$precision
    result$Fscore = SA$Fscore
    result$"FP+" = SA.pos$FP
    result$"FN+" = SA.pos$FN
    result$"TPR+" = SA.pos$TPR
    result$"FP-" = SA.neg$FP
    result$"FN-" = SA.neg$FN
    result$"TPR-" = SA.neg$TPR
  }
  return(unlist(result))
}
getMetricsBM = function(
  y.train = NULL, y.test = NULL, ilrX.train = NULL, ilrX.test = NULL, 
  n.train = NULL, n.test = NULL,
  thetahat0 = NULL, thetahat, betahat, true.sbp = NULL,
  true.beta = NULL, is0.true.beta, non0.true.beta,
  metrics = c("prediction", "betaestimation", "selection"), 
  classification = FALSE
){
  result = list()
  
  if("prediction" %in% metrics){
    if(is.null(y.train) | is.null(y.test) | 
       is.null(ilrX.train) | is.null(ilrX.test) | 
       is.null(n.train) | is.null(n.test) | is.null(thetahat0)) {
      stop("getMetricsBalanceReg() prediction: one or more of these are missing -- y.train, y.test, ilrX.train, ilrX.test, n.train, n.test, thetahat0")
    }
    # 1. prediction error #
    # 1a. on training set #
    result$PEtr = getMSEyhat(
      y = y.train, n = n.train, betahat0 = thetahat0, betahat = thetahat,
      predMat = ilrX.train, classification = classification)
    # 1b. on test set #
    result$PEte = getMSEyhat(
      y = y.test, n = n.test, betahat0 = thetahat0, betahat = thetahat,
      predMat = ilrX.test, classification = classification)
  }
  
  if("betaestimation" %in% metrics){
    if(is.null(true.beta)) {
      stop("getMetricsBalanceReg() betaestimation: true.beta arg is missing")
    }
    # 2. estimation accuracy #
    # 2a. estimation of beta #
    EA = getEstimationAccuracy(true.beta = true.beta, betahat = betahat)
    # 2b. estimation accuracy for active set
    EA.active = getEstimationAccuracy(
      true.beta = true.beta[non0.true.beta], betahat = betahat[non0.true.beta])
    # 2c. estimation accuracy for inactive set
    EA.inactive = getEstimationAccuracy(
      true.beta = true.beta[is0.true.beta], betahat = betahat[is0.true.beta])
    # save
    result$EA1 = EA$EA1
    result$EA2 = EA$EA2
    result$EAInfty = EA$EAInfty
    result$EA1Active = EA.active$EA1
    result$EA2Active = EA.active$EA2
    result$EAInftyActive = EA.active$EAInfty
    result$EA1Inactive = EA.inactive$EA1
    result$EA2Inactive = EA.inactive$EA2
    result$EAInftyInactive = EA.inactive$EAInfty
  }
  
  if("selection" %in% metrics){
    if(is.null(true.beta) & is.null(true.sbp)) {
      stop("getMetricsBalanceReg() selection: true.beta and true.sbp args are both missing, need at least one of them")
    }
    if(is.null(true.beta) & !is.null(true.sbp)){
      true.beta = as.numeric(apply(true.sbp, 1, function(row) any(row != 0)))
    }
    # 3. selection accuracy #
    non0.betahat = abs(betahat) > 10e-8
    SA = getSelectionAccuracy(
      is0.true.beta = is0.true.beta,
      non0.true.beta = non0.true.beta,
      non0.betahat = non0.betahat)
    SA.pos = getSelectionAccuracy(
      is0.true.beta = is0.true.beta[true.beta > 0],
      non0.true.beta = non0.true.beta[true.beta > 0],
      non0.betahat = non0.betahat[true.beta > 0])
    SA.neg = getSelectionAccuracy(
      is0.true.beta = is0.true.beta[true.beta < 0],
      non0.true.beta = non0.true.beta[true.beta < 0],
      non0.betahat = non0.betahat[true.beta < 0])
    # save
    result$FP = SA$FP
    result$FN = SA$FN
    result$TPR = SA$TPR
    result$precision = SA$precision
    result$Fscore = SA$Fscore
    result$"FP+" = SA.pos$FP
    result$"FN+" = SA.pos$FN
    result$"TPR+" = SA.pos$TPR
    result$"FP-" = SA.neg$FP
    result$"FN-" = SA.neg$FN
    result$"TPR-" = SA.neg$TPR
  }
  return(unlist(result))
}

# # solution paths
# roc.for.coef <- function(beta_hat, beta, eps = 1e-08){
#   TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
#   FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
#   FP = sum((abs(beta_hat) > eps) * (abs(beta) <= eps))
#   tpr <- TP/(TP + FN)
#   S_hat <- sum((abs(beta_hat) > eps))
#   prec = TP/(TP + FP)
#   # out <- c(S_hat,tpr)
#   # names(out) <- c('S_hat','tpr')
#   out = c("S_hat" = S_hat, "tpr" = tpr, "TP" = TP, "precision" = prec)
#   return(out)
# }
# roc.for.coef.LR <- function(beta_hat,beta,sbp,eps=1e-08){
#   if (is.null(sbp)){
#     stop('A sequential binary partition tree is needed for roc evaluation!')
#   }
#   # first identify the variable at the LR scale
#   index <- which(abs(beta_hat) > eps)
#   # map to original variable
#   if (length(index)==0){
#     S_hat <- NULL
#   } else  if (length(index)==1){
#     S_hat <- names(which(abs(sbp[,index])>0))
#   } else {
#     S_hat <- names(which(rowSums(abs(sbp[,index]))>0))
#   }
#   S0 <- names(which((abs(beta) > eps)))
#   TP <- intersect(S_hat, S0)
#   tpr <- length(TP)/length(S0) # prob that an actual positive will test positive
#   prec = length(TP)/length(S_hat)
#   # out <- c(length(S_hat),tpr)
#   # names(out) <- c('S_hat','tpr')
#   out = c(
#     "S_hat" = length(S_hat), "tpr" = tpr, "TP" = length(TP), 
#     "precision" = prec)
#   return(out)
# }
