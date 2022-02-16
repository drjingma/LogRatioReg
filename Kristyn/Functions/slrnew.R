slr <- function(x,y){
  
  p <- ncol(x)
  ## Compute pairwise correlation 
  rhoMat <- matrix(0,p,p)
  rownames(rhoMat) <- colnames(rhoMat) <- colnames(x)
  testmat <- rhoMat
  for (j in 1:p){
    for (k in 1:p){
      if (k==j){next}
      else {
        # testmat[j,k] <- cor.test(log(x[,j])-log(x[,k]),y)$p.value
        rhoMat[j,k] <- abs(stats::cor(log(x[,j])-log(x[,k]),y))
      }
    }
  }
  
  # rank-1 approximation and spectral clustering
  rhoMat.svd <- svd(rhoMat)
  rhoMat_approx_1 <-  tcrossprod(rhoMat.svd$u[,1], rhoMat.svd$v[,1]) * rhoMat.svd$d[1]
  index <- which(spectral.clustering(rhoMat_approx_1)==1)
  
  if(length(index) %in% c(1, p - 1)){
    if(length(index) == 1){
      subset1 <- spectral.clustering(rhoMat[-index,-index]) 
    } else if(length(index) == p - 1){
      subset1 <- spectral.clustering(rhoMat[index,index])
    }
    sbp.est = matrix(0, ncol = 1, nrow = p)
    rownames(sbp.est) <- colnames(x)
    sbp.est[match(names(subset1),rownames(sbp.est)),1] <- subset1
    est.balance <- balance::balance.fromSBP(x=x,y=sbp.est)
    coeff <- coefficients(lm(y~est.balance))
    
    # The correct subset should have larger coefficient. 
    ## We refit the linear model on the balance from the correct subset. 
    out <- list()
    out$kernel <- rhoMat
    
    out$index <- subset1
    refit <- lm(y~est.balance[,1])
    
    out$model <- refit
    return(out)
  } else{
    ## Perform spectral clustering on each subset of variables using the original correlation matrix
    subset1 <- spectral.clustering(rhoMat[index,index])
    subset2 <- spectral.clustering(rhoMat[-index,-index]) 
    
    # Since we don't know which subset contains the active variables, 
    ## we first fit a linear model with balances obtained from both subsets. 
    sbp.est <- matrix(0,ncol=2,nrow=p)
    rownames(sbp.est) <- colnames(x)
    sbp.est[match(names(subset1),rownames(sbp.est)),1] <- subset1
    sbp.est[match(names(subset2),rownames(sbp.est)),2] <- subset2
    est.balance <- balance::balance.fromSBP(x=x,y=sbp.est)
    coeff <- coefficients(lm(y~est.balance))
    
    # The correct subset should have larger coefficient. 
    ## We refit the linear model on the balance from the correct subset. 
    out <- list()
    out$kernel <- rhoMat
    
    if ( abs(coeff[2]) > abs(coeff[3]) ){
      # pick subset1
      out$index <- subset1
      refit <- lm(y~est.balance[,1])
    } else {
      # pick subset2
      out$index <- subset2
      refit <- lm(y~est.balance[,2])
    }
    out$model <- refit
    return(out)
  }
}