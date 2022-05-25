# normalize <- function(contrast){
#   if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")
#   
#   lpos <- sum(contrast == 1)
#   lneg <- sum(contrast == -1)
#   const <- sqrt((lpos*lneg)/(lpos+lneg))
#   contrast[contrast==1] = 1/lpos
#   contrast[contrast==-1] = -1/lneg
#   
#   const * contrast
# }

# copied from func_libs #############################################################
#' Transform Samples with the ilr of a Balance
#'
#' @param x A relative abundance matrix with rows as samples (N) and columns as components (D).
#' @param contrast A vector. One column of a serial binary partition matrix
#'  with values [-1, 0, 1] describing D components.
#'
#' @return A transformation of samples for the balance provided.
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

# copied from func_libs #############################################################
spectral.clustering <- function(W, n_eig = 2, zeta = 0) {
  L = graph.laplacian(W,zeta = zeta) # compute graph Laplacian
  ei = eigen(L, symmetric = TRUE)    # Compute the eigenvectors and values of L
  # we will use k-means to cluster the eigenvectors corresponding to
  # the leading smallest eigenvalues
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  obj <- kmeans(ei$vectors[, 1:n_eig], centers = n_eig, nstart = 100, algorithm = "Lloyd")
  if (n_eig==2){
    cl <- 2*(obj$cluster - 1) - 1 
  } else {
    cl <- obj$cluster
  }
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl) 
}

## Potential to speed up the SLR function: do not need to evaluate the feature score each time


#' Updated: this version first screens variables and then constructs the balance on the Aitchison similarity. 
#' @param x a sample by feature data matrix. Note x cannot have zero values.
#' @param y an n by 1 response variable: continuous or binary
#' @param method method of variable screening: could be correlation-based or wald score-based.
#' @param response.type type of response variable: could be survival, continuous, or binary. Currently only continuous and binary responses are allowed.
#' @param threshold a nonnegative constant between 0 and 1. If NULL, then no variable screening is performed.
#' @param s0.perc Factor for denominator of score statistic, between 0 and 1: the percentile of standard deviation values added to the denominator. Default is 0.
slr.screen = function(
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
    Aitchison.similarity <- max(Aitchison.variation) - Aitchison.variation 
    ## Perform spectral clustering
    sbp.est <- spectral.clustering(Aitchison.similarity,zeta = zeta)
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
  
  class(object) <- 'slr.screen'
  return(object)
}

#' @param newdata x n by p of relative abundances
predict.slr.screen <- function(object,newdata = NULL,response.type=c('survival','continuous','binary')){
  # prediction will be based on the canonical space
  if (missing(newdata) || is.null(newdata)) {
    stop('No new data provided!')
  } else {
    if (is.null(object$sbp)){
      predictor <- rep(1,nrow(newdata)) * object$theta
    } else {
      newdata.reduced <- newdata[,colnames(newdata) %in% names(object$sbp)]
      new.balance <- slr.fromContrast(newdata.reduced,object$sbp)
      if (response.type=='binary'){
        fitted.results <- predict(object$fit,newdata=data.frame(balance=new.balance),type='response')
        predictor <- ifelse(fitted.results > 0.5,1,0)
      } else if (response.type=='continuous'){
        predictor <- cbind(1,new.balance) %*% object$theta
      }
    }
    as.numeric(predictor)
  }
}

buildPredmat <- function(outlist,threshold,x,foldid,response.type){
  nfolds = max(foldid)
  predlist = as.list(seq(nfolds))
  predmat = matrix(NA, length(foldid), length(threshold))
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    x_in = x[which, , drop=FALSE]
    predlist[[i]] = sapply(fitobj, function(a) predict.slr.screen(a,newdata=x_in,response.type=response.type))
    predmat[which,] <- predlist[[i]]
  }
  predmat
}

getOptcv <- function(threshold, cvm, cvsd){
  # if(match(cvname,c("AUC","C-index"),0))cvm=-cvm
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  threshold.min = max(threshold[idmin], na.rm = TRUE)
  idmin = match(threshold.min, threshold)
  semin = (cvm + cvsd)[idmin]
  id1se = cvm <= semin
  threshold.1se = max(threshold[id1se], na.rm = TRUE)
  id1se = match(threshold.1se, threshold)
  index=matrix(c(idmin,id1se),2,1,dimnames=list(c("min","1se"),"threshold"))
  list(threshold.min = threshold.min, threshold.1se = threshold.1se, index = index)
}

cv.slr.screen <- function(x,y,method=c('correlation','wald'),
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
        threshold <- sort(stats::pt(abs(feature.scores),df=n-2))
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
    #     lapply(threshold, function(l) slr.screen(x_sub,y_sub,type=type,response.type = response.type,s0.perc=s0.perc,l))
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
      outlist[[i]] <- lapply(threshold, function(l) slr.screen(x_sub,y_sub,method=method,response.type = response.type,threshold=l,s0.perc=s0.perc,zeta = zeta))
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
  class(obj) = "cv.slr.screen"
  obj
}


