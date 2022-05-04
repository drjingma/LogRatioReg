
##### help generate data #######################################################

sigmoid = function(x) 1 / (1 + exp(-x))

alrinv <- function(y) {
  x <- cbind(exp(y), 1)
  x / rowSums(x)
}

geometric.mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}


##### preprocess data for selbal ###############################################

getSelbalData = function(X = X, y = y, classification = TRUE, levels, labels){
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
    if(is.null(levels) | is.null(labels)){
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

getSBPSelbal = function(
    p, selbal.fit){
  
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
  if(!is.null(rownames(X))) rownames(sbp) = rownames(X)
  return(matrix(sbp))
}

getCoefsSelbal = function(
  X, y, selbal.fit, classification = FALSE, check = TRUE){
  p = ncol(X)
  
  # U (transformation) matrix
  sbp = getSBPSelbal(X = X, selbal.fit = selbal.fit)
  # ilrU = rep(0, p)
  # names(ilrU) = colnames(X)
  # pba.pos = unlist(subset(
  #   selbal.fit$global.balance, subset = Group == "NUM", select = Taxa))
  # num.pos = length(pba.pos)
  # pba.neg = unlist(subset(
  #   selbal.fit$global.balance, subset = Group == "DEN", select = Taxa))
  # num.neg = length(pba.neg)
  # ilrU[pba.pos] = 1 / num.pos
  # ilrU[pba.neg] = -1 / num.neg
  is.pos = (sbp == 1)
  is.neg = (sbp == -1)
  num.pos = sum(is.pos)
  num.neg = sum(is.neg)
  ilrU = sbp
  ilrU[is.pos] = 1 / num.pos
  ilrU[is.neg] = -1 / num.neg
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



##### balance regression model tools ###########################################

# Compute the transformation vector/matrix from SBP vector/matrix
getIlrTrans = function(sbp = NULL, detailed = FALSE){
  # check that sbp is in matrix form, even if it's just a p-vector
  # normalize SBP
  pos.vec = colSums(sbp == 1)
  neg.vec = colSums(sbp == -1)
  U = sbp
  for(i in 1:ncol(U)){
    # divide 1s by pos.vec
    U[(U[, i] == 1), i] = 1 / pos.vec[i]
    # divide -1s by neg.vec
    U[(U[, i] == -1), i] = -1 / neg.vec[i]
  }
  # multiply by sqrt(rs / (r + s))
  norm.const = sqrt((pos.vec * neg.vec) / (pos.vec + neg.vec))
  U = sweep(U, MARGIN = 2, STATS = norm.const, FUN = "*")
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



##### plot hierarchical slr ####################################################

# plot the hslr tree of active/inactive variables
#   (not balances! that fluctuates at each level)
#   up to the last tree, which makes a balance
hslrCompleteBTree = function(hslr_fit = NULL, cvhslr_fit = NULL){
  if(is.null(hslr_fit) && is.null(cvhslr_fit)){
    stop("hslrTree(): must provide either hslr_fit or cvhslr_fit.")
  }
  if(!is.null(cvhslr_fit)){
    if(!is.null(hslr_fit)){
      warning("hslrTree(): both hslr_fit and cvhslr_fit were provided -- will overwrite hslr_fit = cvhslr_fit$models")
    }
    hslr_fit = cvhslr_fit$models
  }
  if(is.null(hslr_fit[[1]])){
    stop("hslrTree(): hslr_fit is NULL.")
  }
  component_names = rownames(hslr_fit[[1]]$sbp)
  p = length(component_names)
  cors = t(sapply(hslr_fit, function(level) level$cors))
  Rsqs = t(sapply(hslr_fit, function(level) level$cors))
  sbps = sapply(hslr_fit, function(level) level$sbp)
  rownames(sbps) = component_names
  colnames(sbps) <- rownames(cors) <- rownames(Rsqs) <- paste("T", 1:length(hslr_fit), sep = "")
  
  # fill in the rest of the sbp matrix to get a p x (p - 1) matrix
  sbp_full = matrix(0, nrow = p, ncol = 1)
  rownames(sbp_full) = rownames(sbps)
  colname_counter = 1
  # go through the hslr sbps to fill in the rows in-between
  for(i in 1:ncol(sbps)){
    # current active variables
    active_vars_tmp = as.numeric(which(sbps[, i] != 0))
    # current inactive variables (relative to previous sbp column)
    if(i == 1){
      prev_active_vars.tmp = 1:p
    } else{
      prev_active_vars.tmp = as.numeric(which(sbps[, i - 1] != 0))
    }
    inactive_vars_tmp = setdiff(prev_active_vars.tmp, active_vars_tmp)
    if(length(inactive_vars_tmp) == 1){ # 1 invalid cluster
      # if just 1 variable has been left out, it wasn't a valid cluster
      # let sbp vector indicate which variable was left out
      sbp_tmp = matrix(0, nrow = p, ncol = 1)
      sbp_tmp[inactive_vars_tmp, 1] = -1
      sbp_tmp[active_vars_tmp, 1] = 1
      colnames(sbp_tmp) = colnames(sbps)[i]
      colname_counter = colname_counter + 1
    } else{
      sbp_tmp = matrix(0, nrow = p, ncol = length(inactive_vars_tmp))
      # one column is the active/inactive variables
      sbp_tmp[inactive_vars_tmp, 1] = -1
      sbp_tmp[active_vars_tmp, 1] = 1
      # the rest are the inactive variables, 1 left out at a time
      for(j in 2:ncol(sbp_tmp)){
        sbp_tmp[inactive_vars_tmp[j - 1], j] = 1
        sbp_tmp[inactive_vars_tmp[j:length(inactive_vars_tmp)], j] = -1
      }
      colnames(sbp_tmp) = c(
        colnames(sbps)[i], 
        paste(
          "z", colname_counter + 0:(ncol(sbp_tmp) - 2), sep = ""))
      colname_counter = colname_counter + length(inactive_vars_tmp) - 1
    }
    sbp_full = cbind(sbp_full, sbp_tmp)
    if(i == 1){
      sbp_full = sbp_full[, -1]
    }
    if(i == ncol(sbps) && length(active_vars_tmp) == 2){
      sbp_tmp = matrix(0, nrow = p, ncol = 1)
      sbp_tmp[active_vars_tmp[1], 1] = 1
      sbp_tmp[active_vars_tmp[2], 1] = -1
      colnames(sbp_tmp) = colnames(sbps)[i]
      sbp_full = cbind(sbp_full, sbp_tmp)
    }
  }
  return(list(sbp_complete = sbp_full, Rsqs = Rsqs, cors = cors))
}

plotSBP = function(
  sbp = NULL, edges = NULL, title = NULL, 
  nodes_types = NULL, 
  # a data frame of the nodes/vertices (column = "name") labeled as 
  #   "balance", "covariate", "significant covariate" (column = "type")
  text_size = 2
){
  # get edges data frame, if one isn't given
  if(is.null(sbp) & is.null(edges)){
    stop("plotSBP: provide either sbp or edges arguments!!")
  }
  if(is.null(edges) & !is.null(sbp)){
    edges = getEdgesFromSBP(sbp)
  }
  
  # make the graph input for ggraph plot
  if(!is.null(nodes_types)){
    # check nodes_types data frame
    if(!all(unname(nodes_types[, 1]) %in% unique(unlist(edges))) ||
       !all(unique(unlist(edges)) %in% unname(nodes_types[, 1]))){
      stop("names of nodes in nodes_types data frame don't match the row names of sbp matrix")
    } 
    mygraph <- as_tbl_graph(graph_from_data_frame(
      d = edges, 
      vertices = nodes_types
    ),
    directed = TRUE)
  } else{
    mygraph <- as_tbl_graph(graph_from_data_frame(
      d = edges, vertices = data.frame(nodes = unique(unlist(edges)))),
      directed = TRUE)
  }
  
  # make the plot
  plt = ggraph(mygraph, layout = 'igraph', algorithm = 'tree') +
    geom_edge_diagonal() +
    geom_node_point()
  if(!is.null(title)){
    plt = plt + ggtitle(title)
  }
  if(!is.null(nodes_types)){
    plt = plt + 
      geom_node_label(aes(label = name, fill = type), size = text_size)
  } else{
    plt = plt + 
      geom_node_label(aes(label = name), size = text_size)
  }
  plt = plt + theme_void()
  return(plt)
}

getEdgesFromSBP = function(sbp){
  # make sure the rows (leaf nodes) and & columns (inner nodes) are named
  if(is.null(rownames(sbp))){
    rownames(sbp) = paste("s", 1:nrow(sbp), sep = "")
  }
  if(is.null(colnames(sbp))){
    colnames(sbp) = paste("z", 1:ncol(sbp), sep = "")
  }
  leafs = rownames(sbp)
  inners = colnames(sbp)
  
  # make edges data frame
  edges.df = data.frame(from = character(), to = character())
  for(j in 1:ncol(sbp)){
    from.tmp = rep(inners[j], 2) # because binary tree
    # left child - "1"'s
    if(sum(sbp[, j] == 1) == 1){ # if one "1", it's a leaf node
      left.tmp = leafs[which(sbp[, j] == 1)]
    } else{ # if multiple "1"'s, it's an inner node
      is.left.tmp = rep(FALSE, ncol(sbp))
      for(k in 1:ncol(sbp)){
        # there is only one column with the same locations of non-zeroes as 
        #   there are "1"'s in sbp[, j] -- that column has the left child node
        if(isTRUE(all.equal((sbp[, k] != 0), (sbp[, j] == 1)))){
          is.left.tmp[k] = TRUE
        }
      }
      left.tmp = inners[is.left.tmp]
    }
    # right child - "-1"'s
    if(sum(sbp[, j] == -1) == 1){ 
      right.tmp = leafs[which(sbp[, j] == -1)]
    } else{ 
      is.right.tmp = rep(FALSE, ncol(sbp))
      for(k in 1:ncol(sbp)){
        if(isTRUE(all.equal((sbp[, k] != 0), (sbp[, j] == -1)))){
          is.right.tmp[k] = TRUE
        }
      }
      right.tmp = inners[is.right.tmp]
    }
    to.tmp = c(left.tmp, right.tmp) # left & right children -- swap for ggraph
    if(length(to.tmp) == 2 & length(from.tmp) == 2){
      edges.df = rbind(edges.df, data.frame(from = from.tmp, to = to.tmp))
    }
  }
  return(edges.df)
}





