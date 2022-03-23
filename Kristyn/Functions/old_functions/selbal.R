fitselbal = function(X, y){
  # fit selbal - computes balances internally?
  fit = selbal::selbal(X, y, draw = FALSE)
  design = cbind(rep(1, nrow(X)), fit[[1]])
  balances_lm <- lm(y ~ design) # for predict()
  
  # get contrast (sbp matrix, to compute balances)
  contrast_selbal = matrix(0, ncol(X), 1)
  rownames(contrast_selbal) = colnames(X)
  contrast_selbal[match(fit[[2]], rownames(contrast_selbal)), 1] = 1 # positive
  contrast_selbal[match(fit[[3]], rownames(contrast_selbal)), 1] = -1 # negative
  
  # compute balances
  selbal.predict = function(x){
    design[, 2] = balance.fromSBP(x, contrast_selbal)
    predict(balances_lm, newx = design)
  }
  return(selbal.predict)
}
