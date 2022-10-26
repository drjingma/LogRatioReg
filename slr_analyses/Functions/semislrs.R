
slr = function(
    x,
    y,
    screen.method=c('correlation','wald'),
    cluster.method = c('spectral', 'hierarchical'),
    response.type=c('survival','continuous','binary'),
    threshold,
    s0.perc=0,
    zeta=0,
    x.unlabeled = NULL, use.unlabeled = FALSE,
    positive.slope = FALSE
){
  this.call <- match.call()
  screen.method <- match.arg(screen.method)
  response.type <- match.arg(response.type)
  
  n <- length(y)
  
  feature.scores = getFeatureScores(x, y, screen.method, response.type, s0.perc)
  which.features <- (abs(feature.scores) >= threshold)
  if (sum(which.features)<2){
    # Fit an intercept only regression model
    if (response.type=='binary'){
      model.train <- glm(y~.,data=data.frame(
        y=as.factor(y)),family=binomial(link='logit'))
    } else if (response.type=='continuous'){
      model.train <- lm(y~.,data=data.frame(y=y))
    }
    object <- list(sbp=NULL, Aitchison.var = NULL, cluster.mat = NULL)
  } else {
    x.reduced <- x[,which.features] # reduced data matrix
    if(use.unlabeled){
      if(is.null(x.unlabeled)){
        stop("use.unlabeled is TRUE but x.unlabeled is missing!!")
      }
      if(!all.equal(colnames(x), colnames(x.unlabeled))){
        stop("x and x.unlabeled don't have the same components/covariates/columns! cannot use the unlabeled data.")
      }
      all.x = rbind(x, x.unlabeled)
      all.x.reduced <- all.x[,which.features]
      Aitchison.var = getAitchisonVar(all.x.reduced)
      rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(all.x.reduced)
    } else{
      Aitchison.var = getAitchisonVar(x.reduced)
      rownames(Aitchison.var) <- colnames(Aitchison.var) <- colnames(x.reduced)
    }
    if(cluster.method == "spectral" | nrow(Aitchison.var) == 2){
      Aitchison.sim <- max(Aitchison.var) - Aitchison.var
      ## Perform spectral clustering
      sbp.est <- spectral.clustering(Aitchison.sim, zeta = zeta)
      cluster.mat = Aitchison.sim
    } else if(cluster.method == "hierarchical"){
      ## Perform hierarchical clustering
      htree.est <- hclust(dist(Aitchison.var))
      sbp.est <- balance::sbp.fromHclust(htree.est)[, 1] # grab 1st partition
      cluster.mat = Aitchison.var
    } else{
      stop("invalid cluster.method arg was provided!!")
    }
    balance <- slr.fromContrast(x.reduced, sbp.est) # predict from labeled data
    # model fitting
    if (response.type=='binary'){
      model.train <- glm(
        y~balance,data=data.frame(balance=balance,y=as.factor(y)),
        family=binomial(link='logit'))
      if(positive.slope){
        if(coefficients(model.train)[2] < 0){
          sbp.est = -sbp.est
          balance <- slr.fromContrast(x.reduced, sbp.est)
          model.train <- glm(
            y~balance,data=data.frame(balance=balance,y=as.factor(y)),
            family=binomial(link='logit'))
        }
      }
    } else if (response.type=='continuous'){
      model.train <- lm(
        y~balance,data=data.frame(balance=balance,y=y))
      if(positive.slope){
        if(coefficients(model.train)[2] < 0){
          sbp.est = -sbp.est
          balance <- slr.fromContrast(x.reduced, sbp.est)
          model.train <- lm(
            y~balance,data=data.frame(balance=balance,y=y))
        }
      }
    }
    object <- list(
      sbp = sbp.est, Aitchison.var = Aitchison.var, cluster.mat = cluster.mat)
  }
  object$feature.scores <- feature.scores
  object$theta <- as.numeric(coefficients(model.train))
  object$fit <- model.train
  
  class(object) <- 'slr'
  return(object)
}



cv.slr <- function(
    x, y, screen.method = c('correlation', 'wald'),
    cluster.method = c('spectral', 'hierarchical'),
    response.type=c('survival', 'continuous', 'binary'),
    threshold = NULL,
    s0.perc = 0, zeta = 0,
    x.unlabeled = NULL, use.unlabeled = FALSE, #
    type.measure = c(
      "default", "mse", "deviance", "class", "auc", "mae", "C", "accuracy"
    ),
    scale = FALSE, nfolds = 10,
    foldid = NULL, weights = NULL, #
    fold.unlabeled = FALSE, foldid.unlabeled = NULL, #
    trace.it = FALSE #
){
  type.measure = match.arg(type.measure)
  N <- nrow(x)
  p <- ncol(x)
  N2 = nrow(x.unlabeled)
  
  if (is.null(weights)){
    weights = rep(1, nfolds)
  }
  
  if (is.null(threshold)) {
    xclr <- apply(x,1,function(a) log(a) - mean(log(a)))
    xclr.centered <- base::scale(t(xclr),center=TRUE, scale=TRUE)
    
    # determine threshold based on univariate score statistics or correlations
    threshold = sort(
      getFeatureScores(x, y, screen.method, response.type, s0.perc))
  }
  
  if (is.null(foldid)) {
    foldid = sample(rep(seq(nfolds), length = N))
  } else {
    nfolds = max(foldid)
  }
  if (nfolds < 3){
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  if(fold.unlabeled & is.null(foldid.unlabeled)){
    if(is.null(x.unlabeled)){
      stop("fold.unlabeled is TRUE but x.unlabeled is missing!!")
    }
    if(nfolds > N2){
      stop("nfolds is greater than the number of samples in x.unlabeled! cannot divide into folds.")
    }
    foldid.unlabeled = sample(rep(seq(nfolds), length = N2))
  }
  
  if (trace.it){
    cat("Training\n")
  }
  
  outlist = as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    if (trace.it){
      cat(sprintf("Fold: %d/%d\n", i, nfolds))
    }
    which.fold.i = foldid == i
    # x_in <- x[which.fold.i, ,drop=FALSE]
    x_sub <- x[!which.fold.i, ,drop=FALSE]
    y_sub <- y[!which.fold.i]
    if(use.unlabeled){
      if(is.null(x.unlabeled)){
        stop("use.unlabeled is TRUE but x.unlabeled is missing!!")
      }
      if(fold.unlabeled){
        which.fold.i.unlabeled = foldid.unlabeled == i
        x.unlab_sub = x.unlabeled[!which.fold.i.unlabeled, ,drop=FALSE]
      } else{
        x.unlab_sub = x.unlabeled
      }
    } else{
      x.unlab_sub = NULL
    }
    outlist[[i]] <- lapply(threshold, function(l) slr(
      x = x_sub, y = y_sub,
      screen.method = screen.method, cluster.method = cluster.method,
      response.type = response.type,
      threshold = l,
      s0.perc = s0.perc, zeta = zeta,
      x.unlabeled = x.unlab_sub, use.unlabeled = use.unlabeled))
  }
  
  # collect all out-of-sample predicted values
  #   with the updated code, this is more like a CV matrix
  predmat <- buildPredmat(
    outlist, threshold, x, y, foldid, response.type = response.type,
    type.measure = type.measure)
  
  cvm <- apply(predmat, 2, weighted.mean, w=weights, na.rm = TRUE)
  cvsd = apply(predmat, 2, stats::sd, na.rm = TRUE) / sqrt(nfolds)
  
  out <- list(
    threshold = threshold,
    cvm=cvm,cvsd = cvsd,
    fit.preval = predmat,
    foldid = foldid
  )
  
  lamin <- with(out, getOptcv(threshold, cvm, cvsd))
  
  obj = c(out, as.list(lamin))
  class(obj) = "cv.slr"
  obj
}
