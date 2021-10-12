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
roc.for.coef <- function(beta_hat, beta, eps = 1e-08){
  TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
  FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
  tpr <- TP/(TP + FN)
  S_hat <- sum((abs(beta_hat) > eps))
  # out <- c(S_hat,tpr)
  # names(out) <- c('S_hat','tpr')
  out = c("S_hat" = S_hat, "tpr" = tpr, "TP" = TP)
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
  tpr <- length(TP)/length(S_hat)
  # out <- c(length(S_hat),tpr)
  # names(out) <- c('S_hat','tpr')
  out = c("S_hat" = length(S_hat), "tpr" = tpr, "TP" = length(TP))
  return(out)
}
