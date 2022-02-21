# metrics
getMSEyhat = function(y, n, betahat0, betahat, predMat){
  yhat = betahat0 + predMat %*% betahat
  mse = as.vector(crossprod(y - yhat) / n)
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
getMetricsBalanceReg = function(
  y.train, y.test, ilrX.train, ilrX.test, n.train, n.test, 
  thetahat0, thetahat, betahat, sbp, 
  true.beta, is0.true.beta, non0.true.beta
){
  # 1. prediction error #
  # 1a. on training set #
  PE.train = getMSEyhat(
    y = y.train, n = n.train, betahat0 = thetahat0, betahat = thetahat, 
    predMat = ilrX.train)
  # 1b. on test set #
  PE.test = getMSEyhat(
    y = y.test, n = n.test, betahat0 = thetahat0, betahat = thetahat, 
    predMat = ilrX.test)
  
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  EA = getEstimationAccuracy(true.beta = true.beta, betahat = betahat)
  # 2b. estimation accuracy for active set
  EA.active = getEstimationAccuracy(
    true.beta = true.beta[non0.true.beta], betahat = betahat[non0.true.beta])
  # 2c. estimation accuracy for inactive set
  EA.inactive = getEstimationAccuracy(
    true.beta = true.beta[is0.true.beta], betahat = betahat[is0.true.beta])
  
  # 3. selection accuracy #
  non0.thetahat = (thetahat != 0)
  # sel.cols.SBP = sbp[, non0.thetahat, drop = FALSE]
  # non0.betahat = apply(sel.cols.SBP, 1, function(row) any(row != 0))
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
  
  return(
    c(
      "PEtr" = PE.train, 
      "PEte" = PE.test, 
      "EA1" = EA$EA1, 
      "EA2" = EA$EA2, 
      "EAInfty" = EA$EAInfty, 
      "EA1Active" = EA.active$EA1, 
      "EA2Active" = EA.active$EA2, 
      "EAInftyActive" = EA.active$EAInfty, 
      "EA1Inactive" = EA.inactive$EA1, 
      "EA2Inactive" = EA.inactive$EA2, 
      "EAInftyInactive" = EA.inactive$EAInfty, 
      "FP" = SA$FP, 
      "FN" = SA$FN, 
      "TPR" = SA$TPR, 
      "precision" = SA$precision, 
      "Fscore" = SA$Fscore, 
      "FP+" = SA.pos$FP, 
      "FN+" = SA.pos$FN, 
      "TPR+" = SA.pos$TPR, 
      "precision+" = SA.pos$precision, 
      "Fscore+" = SA.pos$Fscore, 
      "FP-" = SA.neg$FP, 
      "FN-" = SA.neg$FN, 
      "TPR-" = SA.neg$TPR, 
      "precision-" = SA.neg$precision, 
      "Fscore-" = SA.neg$Fscore
    )
  )
}
# getMetricsBalanceReg = function(
#   y.train, y.test, ilrX.train, ilrX.test, n.train, n.test, 
#   thetahat0, thetahat, betahat, sbp, 
#   true.beta, is0.true.beta, non0.true.beta, 
#   metrics = c("prediction", "betaestimation", "selection")
# ){
#   result = list()
#   
#   if("prediction" %in% metrics){
#     # 1. prediction error #
#     # 1a. on training set #
#     result$PEtr = getMSEyhat(
#       y = y.train, n = n.train, betahat0 = thetahat0, betahat = thetahat, 
#       predMat = ilrX.train)
#     # 1b. on test set #
#     result$PEte = getMSEyhat(
#       y = y.test, n = n.test, betahat0 = thetahat0, betahat = thetahat, 
#       predMat = ilrX.test)
#   }
#   
#   if("betaestimation" %in% metrics){
#     # 2. estimation accuracy #
#     # 2a. estimation of beta #
#     EA = getEstimationAccuracy(true.beta = true.beta, betahat = betahat)
#     # 2b. estimation accuracy for active set
#     EA.active = getEstimationAccuracy(
#       true.beta = true.beta[non0.true.beta], betahat = betahat[non0.true.beta])
#     # 2c. estimation accuracy for inactive set
#     EA.inactive = getEstimationAccuracy(
#       true.beta = true.beta[is0.true.beta], betahat = betahat[is0.true.beta])
#     # save
#     result$EA1 = EA$EA1
#     result$EA2 = EA$EA2
#     result$EAInfty = EA$EAInfty
#     result$EA1Active = EA.active$EA1
#     result$EA2Active = EA.active$EA2
#     result$EAInftyActive = EA.active$EAInfty
#     result$EA1Inactive = EA.inactive$EA1
#     result$EA2Inactive = EA.inactive$EA2
#     result$EAInftyInactive = EA.inactive$EAInfty
#   }
#   
#   if("selection" %in% metrics){
#     # 3. selection accuracy #
#     non0.thetahat = (thetahat != 0)
#     # sel.cols.SBP = sbp[, non0.thetahat, drop = FALSE]
#     # non0.betahat = apply(sel.cols.SBP, 1, function(row) any(row != 0))
#     non0.betahat = abs(betahat) > 10e-8
#     SA = getSelectionAccuracy(
#       is0.true.beta = is0.true.beta, 
#       non0.true.beta = non0.true.beta, 
#       non0.betahat = non0.betahat)
#     SA.pos = getSelectionAccuracy(
#       is0.true.beta = is0.true.beta[true.beta > 0], 
#       non0.true.beta = non0.true.beta[true.beta > 0], 
#       non0.betahat = non0.betahat[true.beta > 0])
#     SA.neg = getSelectionAccuracy(
#       is0.true.beta = is0.true.beta[true.beta < 0], 
#       non0.true.beta = non0.true.beta[true.beta < 0], 
#       non0.betahat = non0.betahat[true.beta < 0])
#     result$FP = SA$FP
#     result$FN = SA$FN
#     result$TPR = SA$TPR
#     result$precision = SA$precision
#     result$Fscore = SA$Fscore
#     result$"FP+" = SA.pos$FP
#     result$"FN+" = SA.pos$FN
#     result$"TPR+" = SA.pos$TPR
#     result$"precision+" = SA.pos$precision
#     result$"Fscore+" = SA.pos$Fscore
#     result$"FP-" = SA.neg$FP
#     result$"FN-" = SA.neg$FN
#     result$"TPR-" = SA.neg$TPR
#     result$"precision-" = SA.neg$precision
#     result$"Fscore-" = SA.neg$Fscore
#   }
#   
#   return(unlist(result))
# }
getMetricsLLC = function(
  y.train, y.test, logX.train, logX.test, n.train, n.test, 
  betahat0, betahat, true.beta, is0.true.beta, non0.true.beta
){
  # 1. prediction error #
  # 1a. on training set #
  PE.train = getMSEyhat(
    y = y.train, n = n.train, betahat0 = betahat0, betahat = betahat, 
    predMat = logX.train)
  # 1b. on test set #
  PE.test = getMSEyhat(
    y = y.test, n = n.test, betahat0 = betahat0, betahat = betahat, 
    predMat = logX.test)
  
  # 2. estimation accuracy #
  # 2a. estimation of beta #
  EA = getEstimationAccuracy(true.beta = true.beta, betahat = betahat)
  # 2b. estimation accuracy for active set
  EA.active = getEstimationAccuracy(
    true.beta = true.beta[non0.true.beta], betahat = betahat[non0.true.beta])
  # 2c. estimation accuracy for inactive set
  EA.inactive = getEstimationAccuracy(
    true.beta = true.beta[is0.true.beta], betahat = betahat[is0.true.beta])
  
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
  
  return(
    c(
      "PEtr" = PE.train, 
      "PEte" = PE.test, 
      "EA1" = EA$EA1, 
      "EA2" = EA$EA2, 
      "EAInfty" = EA$EAInfty, 
      "EA1Active" = EA.active$EA1, 
      "EA2Active" = EA.active$EA2, 
      "EAInftyActive" = EA.active$EAInfty, 
      "EA1Inactive" = EA.inactive$EA1, 
      "EA2Inactive" = EA.inactive$EA2, 
      "EAInftyInactive" = EA.inactive$EAInfty, 
      "FP" = SA$FP, 
      "FN" = SA$FN, 
      "TPR" = SA$TPR, 
      "precision" = SA$precision, 
      "Fscore" = SA$Fscore, 
      "FP+" = SA.pos$FP, 
      "FN+" = SA.pos$FN, 
      "TPR+" = SA.pos$TPR, 
      "precision+" = SA.pos$precision, 
      "Fscore+" = SA.pos$Fscore, 
      "FP-" = SA.neg$FP, 
      "FN-" = SA.neg$FN, 
      "TPR-" = SA.neg$TPR, 
      "precision-" = SA.neg$precision, 
      "Fscore-" = SA.neg$Fscore
    )
  )
}

# solution paths
roc.for.coef <- function(beta_hat, beta, eps = 1e-08){
  TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
  FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
  FP = sum((abs(beta_hat) > eps) * (abs(beta) <= eps))
  tpr <- TP/(TP + FN)
  S_hat <- sum((abs(beta_hat) > eps))
  prec = TP/(TP + FP)
  # out <- c(S_hat,tpr)
  # names(out) <- c('S_hat','tpr')
  out = c("S_hat" = S_hat, "tpr" = tpr, "TP" = TP, "precision" = prec)
  return(out)
}
roc.for.coef.LR <- function(beta_hat,beta,sbp,eps=1e-08){
  if (is.null(sbp)){
    stop('A sequential binary partition tree is needed for roc evaluation!')
  }
  # first identify the variable at the LR scale
  index <- which(abs(beta_hat) > eps)
  # map to original variable
  if (length(index)==0){
    S_hat <- NULL
  } else  if (length(index)==1){
    S_hat <- names(which(abs(sbp[,index])>0))
  } else {
    S_hat <- names(which(rowSums(abs(sbp[,index]))>0))
  }
  S0 <- names(which((abs(beta) > eps)))
  TP <- intersect(S_hat, S0)
  tpr <- length(TP)/length(S0) # prob that an actual positive will test positive
  prec = length(TP)/length(S_hat)
  # out <- c(length(S_hat),tpr)
  # names(out) <- c('S_hat','tpr')
  out = c(
    "S_hat" = length(S_hat), "tpr" = tpr, "TP" = length(TP), 
    "precision" = prec)
  return(out)
}
