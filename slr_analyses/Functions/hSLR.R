# slr using hierarchical clustering on Aitchison variation matrix A
#   instead of spectral clustering on Aitchison similarity matrix max(A) - A
hslr = function( ###
    x, 
    y, 
    method=c('correlation','wald'),
    response.type=c('survival','continuous','binary'),
    threshold,
    s0.perc=0,
    zeta=0
){
  this.call <- match.call()
  method <- match.arg(method)
  response.type <- match.arg(response.type)
  
  n <- length(y)
  
  ## Compute univariate coefficient 
  xclr <- apply(x,1,function(a) log(a) - mean(log(a)))
  if (method=='wald'){
    xclr.centered <- scale(t(xclr),center=TRUE,scale=TRUE)
    if (response.type=='continuous'){
      y.centered <- y - mean(y)
      sxx <- diag(crossprod(xclr.centered))
      sxy <- crossprod(xclr.centered,y.centered)
      syy <- sum(y.centered^2)
      numer <- sxy/sxx
      sd <- sqrt((syy/sxx - numer^2)/(n - 2)) # variance estimate
      
      if (is.null(s0.perc)) {
        fudge <- median(sd)
      }
      if (!is.null(s0.perc)) {
        if (s0.perc >= 0) {
          fudge <- quantile(sd, s0.perc)
        }
        if (s0.perc < 0) {
          fudge <- 0
        }
      }
      feature.scores <- numer/(sd + fudge) # this is the wald-statistic 
      feature.scores <- stats::pt(abs(feature.scores),df=n-2) # t-distribution with n-2 df
    } else if (response.type=='binary'){
      feature.scores <- rep(0,ncol(xclr.centered))
      for (j in 1:ncol(xclr.centered)){
        fit <- glm(y~x,data=data.frame(x=xclr.centered[,j],y=as.factor(y)),family=binomial(link='logit'))
        feature.scores[j] <- coef(summary(fit))[2,3]
      }
      feature.scores <- stats::pnorm(abs(feature.scores))
    }
  }
  if (method=='correlation'){
    feature.scores <- cor(t(xclr),y)
  }
  
  which.features <- (abs(feature.scores) >= threshold)
  if (sum(which.features)<2){
    # Fit an intercept only regression model 
    if (response.type=='binary'){
      model.train <- glm(y~.,data=data.frame(y=as.factor(y)),family=binomial(link='logit'))
    } else if (response.type=='continuous'){
      model.train <- lm(y~.,data=data.frame(y=y))
    }
    object <- list(sbp=NULL)
  } else {
    x.reduced <- x[,which.features] # reduced data matrix
    p <- ncol(x.reduced)
    Aitchison.variation <- matrix(0,p,p)
    rownames(Aitchison.variation) <- colnames(Aitchison.variation) <- colnames(x.reduced)
    for (j in 1:p){
      for (k in 1:p){
        if (k==j){next}
        else {
          Aitchison.variation[j,k] <- stats::var(log(x.reduced[,j])-log(x.reduced[,k])) # Aitchison variation
        }
      }
    }
    ## Perform hierarchical clustering
    htree.est <- hclust(Aitchison.variation) ###
    sbp.est <- balance::sbp.fromHclust(htree.est) # grab 1st partition ###
    balance <- slr.fromContrast(x.reduced,sbp.est)
    # model fitting
    if (response.type=='binary'){
      model.train <- glm(y~balance,data=data.frame(balance=balance,y=as.factor(y)),family=binomial(link='logit'))
      
    } else if (response.type=='continuous'){
      model.train <- lm(y~balance,data=data.frame(balance=balance,y=y))
    }
    object <- list(sbp = sbp.est)
  }
  object$feature.scores <- feature.scores 
  object$theta <- as.numeric(coefficients(model.train))
  object$fit <- model.train
  
  class(object) <- 'hslr' ###
  return(object)
}

cv.hslr <- function(x,y,method=c('correlation','wald'),
                   response.type=c('survival','continuous','binary'),
                   threshold=NULL,s0.perc=NULL,zeta=0,
                   nfolds=10,foldid=NULL,weights=NULL,
                   type.measure = c("default", "mse", "deviance", "class", "auc", "mae", "C"),
                   parallel=FALSE,scale=FALSE,trace.it=TRUE){
  type.measure = match.arg(type.measure)
  N <- nrow(x)
  p <- ncol(x)
  
  if (is.null(weights))
    weights <- rep(1,N)
  
  if (is.null(threshold)) {
    xclr <- apply(x,1,function(a) log(a) - mean(log(a)))
    xclr.centered <- base::scale(t(xclr),center=TRUE, scale=TRUE)
    # determine threshold based on univariate score statistics or correlations    
    if (method=='correlation'){
      threshold = sort(abs(cor(xclr.centered,as.numeric(y))))
    } else if (method=='wald'){
      if (response.type=='continuous'){
        y.centered <- y - mean(y)
        sxx <- diag(crossprod(xclr.centered))
        sxy <- crossprod(xclr.centered,y.centered)
        syy <- sum(y.centered^2)
        numer <- sxy/sxx
        sd <- sqrt((syy/sxx - numer^2)/(N - 2)) # variance estimate
        
        if (is.null(s0.perc)) {
          fudge <- median(sd)
        }
        if (!is.null(s0.perc)) {
          if (s0.perc >= 0) {
            fudge <- quantile(sd, s0.perc)
          }
          if (s0.perc < 0) {
            fudge <- 0
          }
        }
        feature.scores <- numer/(sd + fudge) # this is the wald-statistic
        threshold <- sort(stats::pt(abs(feature.scores),df=N-2))
      } else if (response.type=='binary'){
        feature.scores <- rep(0,p)
        for (j in 1:p){
          fit <- glm(y~x,data=data.frame(x=xclr.centered[,j],y=as.factor(y)),family=binomial(link='logit'))
          feature.scores[j] <- coef(summary(fit))[2,3]
        }
        threshold <- sort(stats::pnorm(abs(feature.scores)))
      }
    }
  }
  
  if (is.null(foldid)) 
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  
  if (trace.it) cat("Training\n")
  
  outlist = as.list(seq(nfolds))
  if (parallel) {
    # require(foreach)
    # cl <- makeCluster(detectCores())
    # 
    # registerDoParallel(cl)
    # outlist = foreach(i = seq(nfolds), .packages = c("MASS")) %dopar%
    #   {
    #     which = foldid == i
    #     x_in = lapply(x,function(a) a[which, ,drop=FALSE])
    #     x_sub <- lapply(x,function(a) a[!which, ,drop=FALSE])
    #     y_sub <- y[!which]
    #     
    #     lapply(threshold, function(l) slr(x_sub,y_sub,type=type,response.type = response.type,s0.perc=s0.perc,l))
    #   }# end foreach
    # stopCluster(cl)
  }
  else {
    for (i in seq(nfolds)) {
      if (trace.it) cat(sprintf("Fold: %d/%d\n", i, nfolds))
      which = foldid == i
      x_in <- x[which, ,drop=FALSE]
      x_sub <- x[!which, ,drop=FALSE]
      y_sub <- y[!which]
      outlist[[i]] <- lapply(threshold, function(l) hslr(x_sub,y_sub,method=method,response.type = response.type,threshold=l,s0.perc=s0.perc,zeta = zeta))
    }# end loop i
  }
  
  # collect all out-of-sample predicted values
  predmat <- buildPredmat(outlist, threshold, x, foldid, response.type = response.type)
  
  if (response.type=='continuous'){
    mse <- (as.numeric(y)-predmat)^2
  } else if (response.type=='binary'){
    mse <- (predmat != y)
  }
  cvm <- apply(mse, 2, weighted.mean, w=weights, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(mse, cvm, FALSE)^2, 2, weighted.mean, w = weights, na.rm = TRUE)/(N - 1))
  
  out <- list(threshold=threshold,cvm=cvm,cvsd=cvsd,fit.preval = predmat, foldid = foldid)
  
  lamin <- with(out,getOptcv(threshold, cvm, cvsd))
  
  obj = c(out, as.list(lamin))
  class(obj) = "cv.hslr" ###
  obj
}


