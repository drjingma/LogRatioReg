# Copied from RCode/SLR.R
#   added hierarchical clustering option
#   added option to include extra non-compositional covariate(s)
# Date: 6/19/2022


## Potential to speed up the SLR function: do not need to evaluate the feature score each time

#' Updated: this version first screens variables and then constructs the balance on the Aitchison similarity. 
#' @param x a sample by feature data matrix. Note x cannot have zero values.
#' @param y an n by 1 response variable: continuous or binary
#' @param screen.method method of variable screening: could be correlation based or wald score based.
#' @param cluster.method method of clustering.
#' @param response.type type of response variable: could be survival, continuous, or binary. Currently only continuous and binary responses are allowed.
#' @param s0.perc Factor for denominator of score statistic, between 0 and 1: the percentile of standard deviation values added to the denominator. Default is 0.
#' @param threshold a nonnegative constant between 0 and 1. If NULL, then no variable screening is performed.
slr = function(
    x, 
    y, 
    screen.method=c('correlation','wald'),
    cluster.method = c('spectral', 'hierarchical'),
    response.type=c('survival','continuous','binary'),
    threshold,
    s0.perc=0,
    zeta=0
){
  this.call <- match.call()
  screen.method <- match.arg(screen.method)
  response.type <- match.arg(response.type)
  
  n <- length(y)
  
  ## Compute univariate coefficient 
  xclr <- apply(x,1,function(a) log(a) - mean(log(a)))
  if (screen.method=='wald'){
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
      # this is the wald-statistic 
      feature.scores <- numer/(sd + fudge) 
      # t-distribution with n-2 df
      feature.scores <- stats::pt(abs(feature.scores),df=n-2)
    } else if (response.type=='binary'){
      feature.scores <- rep(0,ncol(xclr.centered))
      for (j in 1:ncol(xclr.centered)){
        fit <- glm(y~x,data=data.frame(
          x=xclr.centered[,j],y=as.factor(y)),family=binomial(link='logit'))
        feature.scores[j] <- coef(summary(fit))[2,3]
      }
      feature.scores <- stats::pnorm(abs(feature.scores))
    }
  }
  if (screen.method=='correlation'){
    feature.scores <- cor(t(xclr),y)
  }
  
  which.features <- (abs(feature.scores) >= threshold)
  if (sum(which.features)<2){
    # Fit an intercept only regression model 
    if (response.type=='binary'){
      model.train <- glm(y~.,data=data.frame(
        y=as.factor(y)),family=binomial(link='logit'))
    } else if (response.type=='continuous'){
      model.train <- lm(y~.,data=data.frame(y=y))
    }
    object <- list(sbp=NULL)
  } else {
    x.reduced <- x[,which.features] # reduced data matrix
    p <- ncol(x.reduced)
    Aitchison.var <- matrix(0,p,p)
    rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(x.reduced)
    for (j in 1:p){
      for (k in 1:p){
        if (k==j){next}
        else {
          Aitchison.var[j,k] <- stats::var(log(x.reduced[,j])-log(x.reduced[,k])) # Aitchison variation
        }
      }
    }
    if(cluster.method == "spectral" | nrow(Aitchison.var) == 2){
      Aitchison.sim <- max(Aitchison.var) - Aitchison.var 
      ## Perform spectral clustering
      sbp.est <- spectral.clustering(Aitchison.sim,zeta = zeta)
    } else if(cluster.method == "hierarchical"){
        ## Perform hierarchical clustering
        htree.est <- hclust(dist(Aitchison.var))
        sbp.est <- balance::sbp.fromHclust(htree.est)[, 1] # grab 1st partition
    } else{
      stop("invalid cluster.method arg was provided!!")
    }
    balance <- slr.fromContrast(x.reduced,sbp.est)
    # model fitting
    if (response.type=='binary'){
      model.train <- glm(
        y~balance,data=data.frame(balance=balance,y=as.factor(y)),
        family=binomial(link='logit'))
      
    } else if (response.type=='continuous'){
      model.train <- lm(
        y~balance,data=data.frame(balance=balance,y=y))
    }
    object <- list(sbp = sbp.est)
  }
  object$feature.scores <- feature.scores 
  object$theta <- as.numeric(coefficients(model.train))
  object$fit <- model.train
  
  class(object) <- 'slr'
  return(object)
}

#' @param newdata x n by p of relative abundances
predict.slr <- function(
  object,newdata = NULL,response.type=c('survival','continuous','binary')
){
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
        # predictor <- ifelse(fitted.results > 0.5,1,0)
        predictor = fitted.results
      } else if (response.type=='continuous'){
        predictor <- cbind(1,new.balance) %*% object$theta
      }
    }
    as.numeric(predictor)
  }
}

# # this version doesn't summarize by fold, instead it's n x n.thresh
# buildPredmat <- function(outlist,threshold,x,foldid,response.type){
#   nfolds = max(foldid)
#   predlist = as.list(seq(nfolds))
#   predmat = matrix(NA, length(foldid), length(threshold))
#   for (i in seq(nfolds)) {
#     which = foldid == i
#     fitobj = outlist[[i]]
#     x_in = x[which, , drop=FALSE]
#     predlist[[i]] = sapply(
#       fitobj, function(a) predict.slr(
#         a,newdata=x_in,response.type=response.type))
#     predmat[which,] <- predlist[[i]]
#   }
#   predmat
# }
# this version is nfolds x n.thresh, like the other methods
buildPredmat <- function(
    outlist,threshold,x,y,foldid,response.type,type.measure
){
  nfolds = max(foldid)
  # predlist = as.list(seq(nfolds))
  predmat = matrix(NA, nfolds, length(threshold))
  for(i in 1:nfolds){ # predict for each fold
    which = foldid == i
    y.i = y[which]
    fitobj = outlist[[i]]
    x_in = x[which, , drop=FALSE]
    # predlist[[i]] = sapply(
    #   fitobj, function(a) predict.slr(
    #     a,newdata=x_in,response.type=response.type))
    predy.i = sapply(
      fitobj, function(a) predict.slr(
        a,newdata=x_in,response.type=response.type))
    for(j in 1:length(threshold)){
      predy.ij = predy.i[, j]
      if (response.type=='continuous'){ # mse, to be minimized
        predmat[i, j] <- mean((as.numeric(y.i)-predy.ij)^2)
      } else if (response.type=='binary'){
        if(type.measure == "accuracy"){# accuracy, minimize the # that don't match
          predmat[i, j] <- mean((predy.ij > 0.5) != y.i) 
        } else if(type.measure == "auc"){# auc, minimize 1 - auc 
          predmat[i, j] = 1 - pROC::auc(
            y.i,predy.ij, levels = c(0, 1), direction = "<", quiet = TRUE)
        }
      }
    }
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
  list(
    threshold.min = threshold.min, 
    threshold.1se = threshold.1se, 
    index = index
  )
}

cv.slr <- function(
    x,y,screen.method=c('correlation','wald'),
    cluster.method = c('spectral', 'hierarchical'),
    response.type=c('survival','continuous','binary'),
    threshold=NULL,s0.perc=NULL,zeta=0,
    nfolds=10,foldid=NULL,weights=NULL,
    type.measure = c(
      "default", "mse", "deviance", "class", "auc", "mae", "C", "accuracy"
    ),
    parallel=FALSE,scale=FALSE,trace.it=FALSE
){
  type.measure = match.arg(type.measure)
  N <- nrow(x)
  p <- ncol(x)
  
  if (is.null(weights)){
    # weights <- rep(1,N)
    weights = rep(1, nfolds)
  }
  
  if (is.null(threshold)) {
    xclr <- apply(x,1,function(a) log(a) - mean(log(a)))
    xclr.centered <- base::scale(t(xclr),center=TRUE, scale=TRUE)
    # determine threshold based on univariate score statistics or correlations    
    if (screen.method=='correlation'){
      threshold = sort(abs(cor(xclr.centered,as.numeric(y))))
    } else if (screen.method=='wald'){
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
          fit <- glm(y~x,data=data.frame(
            x=xclr.centered[,j],y=as.factor(y)),family=binomial(link='logit'))
          feature.scores[j] <- coef(summary(fit))[2,3]
        }
        threshold <- sort(stats::pnorm(abs(feature.scores)))
      }
    }
  }
  
  if (is.null(foldid)) {
    foldid = sample(rep(seq(nfolds), length = N))
  } else {
    nfolds = max(foldid)
  }
  
  if (nfolds < 3){
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  
  if (trace.it){
    cat("Training\n")
  }
  
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
  } else {
    for (i in seq(nfolds)) {
      if (trace.it){
        cat(sprintf("Fold: %d/%d\n", i, nfolds))
      }
      which = foldid == i
      x_in <- x[which, ,drop=FALSE]
      x_sub <- x[!which, ,drop=FALSE]
      y_sub <- y[!which]
      outlist[[i]] <- lapply(threshold, function(l) slr(
        x_sub,y_sub,screen.method=screen.method,cluster.method=cluster.method,
        response.type = response.type,threshold=l,s0.perc=s0.perc,zeta = zeta))
    }# end loop i
  }
  
  # collect all out-of-sample predicted values
  predmat <- buildPredmat(
    outlist, threshold, x, y, foldid, response.type = response.type, 
    type.measure = type.measure)

  cvm <- apply(predmat, 2, weighted.mean, w=weights, na.rm = TRUE)
  # cvsd <- sqrt(apply(
  #   scale(predmat, cvm, FALSE)^2, 2, weighted.mean, w = weights, na.rm = TRUE) / 
  #     (N - 1))
  cvsd = apply(predmat, 2, stats::sd, na.rm = TRUE) / sqrt(nfolds)
  
  out <- list(
    threshold=threshold,
    cvm=cvm,cvsd=cvsd,
    fit.preval = predmat, 
    foldid = foldid
  )
  
  lamin <- with(out,getOptcv(threshold, cvm, cvsd))
  
  obj = c(out, as.list(lamin))
  class(obj) = "cv.slr"
  obj
}


