library(gdata)
library(glmnet)
library(Rcpp)
library(dirmult)
library(MGLM)
library(MASS)
library(compositions)
library(matrixStats)

# source("COAT-master/coat.R")
# solves t(beta)%*%xx%*%beta/2 - t(xy)%*%beta + fac||beta||_1 
#or ||X%*%bet-y||_2^2/2 + fac||beta||_1 with xx=crossprod(X), xy=crossprod(x,y)
# subject to t(cmat)%*%beta=0
cppFunction('NumericVector ConstrLassoC0(NumericVector betstart, NumericMatrix xx, NumericVector xy, NumericMatrix cmat, double fac, int maxiter, double tol){
int p = xx.nrow();
int m = cmat.ncol();
NumericVector bet(p);
NumericVector bet_old(p);
NumericVector ksi(m);
NumericVector ksi_old(m);
NumericMatrix xxcc(p,p);
double tmp;
double tmp2;
double fac2;
double myabs;
int mysgn;
int iter;
int iter2;
int k;
double eps;
double eps2;
LogicalVector nonzero(p);
IntegerVector ind(p);


for (int i=0; i<p; i++){
    bet[i] = betstart[i];
    nonzero[i] = xx(i,i)>0;
    ind[i] = i;
}

for (int i=0; i<m; i++){
    ksi[i] = 0;
}

for (int i=0; i<p; i++){
    for (int j=0; j<p; j++){
        xxcc(i,j) = xx(i,j);
        for (int ij=0; ij<m; ij++){
            xxcc(i,j) += cmat(i,ij)*cmat(j,ij);
        }
    }
}

iter2 = 1;
eps2 = 1;
while (eps2>tol & iter2<maxiter){
    iter = 1;
    eps = 1;
    while (eps>tol & iter<maxiter){
        for (int i=0; i<p; i++){
            bet_old[i] = bet[i];
        }
        std::random_shuffle(ind.begin(), ind.end());
        for (int i=0; i<p; i++){
            k = ind[i];
            if(nonzero[k]){
                fac2 = fac/xxcc(k,k);
                tmp = 0;
                for (int j=0; j<p; j++){
                    tmp += xxcc(k,j)*bet[j];
                }
                tmp = tmp - xxcc(k,k)*bet[k];
                tmp2 = 0;
                for (int j=0; j<m; j++){
                    tmp2 += cmat(k,j)*ksi[j];
                }
                tmp = (xy[k] - tmp2 - tmp)/xxcc(k,k);
                if (tmp>0){
                    mysgn = 1;
                    myabs = tmp;
                }else if (tmp<0){
                    mysgn = -1;
                    myabs = -tmp;
                }else{
                    mysgn = 0;
                    myabs = 0;
                }
                if (myabs > fac2){
                    bet[k] = mysgn*(myabs - fac2);
                }else{
                    bet[k] = 0;
                }
            }else{
                 bet[k] = 0;
            }
        }
        eps = 0;
        for (int i=0; i<p; i++){
            eps += pow(bet[i] - bet_old[i], 2.0);
        }
        eps = sqrt(eps);
        iter += 1;
    }
    for (int i=0; i<m; i++){
        ksi_old[i] = ksi[i];
        for (int j=0; j<p; j++){
            ksi[i] += cmat(j,i)*bet[j];
        }
    }
    eps2 = 0;
    for (int j=0; j<m; j++){
        eps2 += pow(ksi[j] - ksi_old[j], 2.0);
    }
    eps2 = sqrt(eps2);
    iter2 += 1;  
}
if(iter2==maxiter){
    for (int i=0; i<p; i++){
        bet[i] = 0;
    }
}
return bet;
}')


## ---- Constrained Lasso ----
#' Function to perform Lasso regression subject to the constraint on the regression coefficients.
#' @param y The outcome vector
#' @param x The n by p design matrix
#' @param Cmat The constraints
#' @param lambda The tuning parameter
ConstrLasso0 <- function(y, x, Cmat, lambda, betstart=NULL, maxiter=1000, tol=1e-8){
  # solves ||y-X%*%beta||_2^2/(2n) + lam||beta||_1 subject to t(Cmat)%*%beta=0
  Cmat <- as.matrix(Cmat)
  XX <- crossprod(x)
  Xy <- crossprod(x,y)
  fac <- length(y)*lambda
  if (is.null(betstart)) betstart <- rep(0, ncol(x))
  set.seed(0)
  bet <- ConstrLassoC0(betstart, XX, Xy, Cmat, fac, maxiter, tol) 
  return(bet)
}

ConstrLasso <- function(y, x, Cmat=NULL, lambda=NULL, nlam=20, intercept=TRUE, scaling=TRUE, maxiter=1000, tol=1e-8){
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(Cmat)) Cmat <- matrix(0,p,1)
  
  if (intercept){
    y.mean <- mean(y)
    y <- y - y.mean
    x.mean <- colMeans(x)
    x <- scale(x, center=x.mean, scale=F) # centering before scaling
  }
  if (scaling){
    x.sd <- apply(x,2,sd)
    x <- scale(x, center=F, scale=x.sd)
    Cmat <- Cmat/x.sd
  }
  
  # get sequence of tuning parameter lambda
  if (is.null(lambda)){
    maxlam <- 2*max(abs(crossprod(x,y)/n))
    lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))
  }
  nlam <- length(lambda)
  
  Cmat <- as.matrix(Cmat)
  xx <- crossprod(x)
  xy <- crossprod(x,y)
  fac <- length(y)*lambda
  betstart <- rep(0, ncol(x))
  set.seed(0)
  bet <- mapply(ConstrLassoC0, fac=fac, MoreArgs=list(xy=xy, xx=xx, cmat=Cmat, betstart=betstart, maxiter=maxiter, tol=tol))
  
  if (scaling){
    bet <- bet/x.sd
  }
  int <- rep(0, nlam)
  if (intercept){
    int <- y.mean - as.vector(x.mean%*%bet)
  }
  return(list(int=int, bet=bet, lambda=lambda))
}

cv.func <- function(method="ConstrLasso", y, x, Cmat=NULL, lambda=NULL, nlam=20, intercept=TRUE, scaling=TRUE, nfolds=5, maxiter=1000, tol=1e-8, seed=0){
  if(is.na(match(method, c("ConstrLasso")))) stop("input method is wrong!")
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(Cmat)) Cmat <- matrix(0,p,1)
  
  if (intercept){
    y.mean <- mean(y)
    y <- y - y.mean
    x.mean <- colMeans(x)
    x <- scale(x, center=x.mean, scale=F) # centering before scaling
  }
  if (scaling){
    x.sd <- apply(x,2,sd)
    x <- scale(x, center=F, scale=x.sd)
    Cmat <- Cmat/x.sd
  }
  
  # get sequence of tuning parameter lambda
  if (is.null(lambda)){
    maxlam <- 2*max(abs(crossprod(x,y)/n))
    lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))
  }
  nlam <- length(lambda)
  
  
  # run the cv folds
  err <- matrix(0, nlam, nfolds)
  set.seed(seed)
  rnum <- sample.int(n)%%nfolds+1
  for (j in 1:nfolds){
    temp <- rnum==j
    print(j)
    ytrain <- y[!temp]
    xtrain <- x[!temp,]
    ytest <- y[temp]
    xtest <- x[temp,]
    res <- do.call(method, list(y=ytrain, x=xtrain, Cmat=Cmat, lambda=lambda, intercept=FALSE, scaling=FALSE, maxiter=maxiter, tol=tol))
    err[,j] <- colMeans((xtest%*%res$bet-res$int-as.vector(ytest))^2)
  }
  cvm <- rowMeans(err)
  cvsd <- apply(err,1,sd)
  # fit with all lambda
  res.fit <- do.call(method, list(y=y, x=x, Cmat=Cmat, lambda=lambda, intercept=FALSE, scaling=FALSE, maxiter=maxiter, tol=tol))
  bet <- res.fit$bet
  if (scaling){
    bet <- bet/x.sd
  }
  int <- rep(0, nlam)
  if (intercept){
    int <- y.mean - as.vector(x.mean%*%bet)
  }  
  
  return(list(lambda=lambda, int=int, bet=bet, cvm=cvm, cvsd=cvsd))
}

#' Function to generate random observations from a logistic normal distribution.
#' 
#' @param n sample size
#' @param mu mean vector
#' @param Sigma covariance matrix 
rlogisnormal <- function(n,mu,Sigma){
  d <- MASS::mvrnorm(n=n,mu=mu,Sigma=Sigma)
  d <- matrix(d, nrow=n)
  d <- cbind(d,rep(0,n))
  d <- exp(d)
  compositions::acomp(d)
}

#' Function to calculate the geometric mean of a vector.
#' @details Zero values in the input vectors are excluded in computing the geometric mean.
gm.mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

#' Function to compute the centered log-ratio transformation of a vector.
#' @details The centered log-ratio transformation is defined as log(X_j/g(X)) where g(X) refers to the geometric mean.
clr <- function(x, base=exp(1)){
  x <- log((x / gm.mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

#' Function to calculate the rand index between two cluster assignments.
#' @details The rand index is a measure to assess the consistency between the true cluster assignment (x) 
#' and the estimated cluster assignment (y). For details, see https://en.wikipedia.org/wiki/Rand_index.
rand.index <- function(x,y){
  sum <- 0
  len <- length(x)
  for(i in 1:len)
  {
    for(j in 1:len)
    {
      if(i!=j)
      {
        bool_1 <- x[i]==x[j]
        bool_2 <- y[i]==y[j]
        if((bool_1&&bool_2)||(!bool_1&&!bool_2))
          sum <- sum+1
      }
    }
  }
  return( sum/(len*(len-1)) )
}

#' Function to calculate the adjusted rand index between two cluster assignments. 
#' @details The adjusted Rand index is the corrected-for-chance version of the Rand index.
adjusted.rand.index <- function(x,y){
  x <- as.vector(x)
  y <- as.vector(y)
  xx <- outer(x, x, "==")
  yy <- outer(y, y, "==")
  upper <- row(xx) < col(xx)
  xx <- xx[upper]
  yy <- yy[upper]
  a <- sum(as.numeric(xx & yy))
  b <- sum(as.numeric(xx & !yy))
  c <- sum(as.numeric(!xx & yy))
  d <- sum(as.numeric(!xx & !yy))
  ni <- (b + a)
  nj <- (c + a)
  abcd <- a + b + c + d
  q <- (ni * nj)/abcd
  (a - q)/((ni + nj)/2 - q)
}

#' Function to calculate the graph Laplacian matrix
#' @param W The weighted adjacency matrix
#' @details https://en.wikipedia.org/wiki/Laplacian_matrix
graph.laplacian <- function(W, normalized = TRUE, zeta=0.01){
  stopifnot(nrow(W) == ncol(W)) 
  
  n = nrow(W)    # number of vertices
  # We perturb the network by adding some links with low edge weights
  W <- W + zeta * mean(colSums(W))/n * tcrossprod(rep(1,n))
  g <- colSums(W) # degrees of vertices
  
  if(normalized){
    D_half = diag(1 / sqrt(g) )
    return(D_half %*% W %*% D_half )
  } else {
    return(W)
  }
}

#' Function to perform spectral clustering from a graph. 
#' 
#' @param W The weighted adjacency matrix
#' @param n_eig Number of clusters/eigenvectors to obtain
spectral.clustering <- function(W, n_eig = 2) {
  L = graph.laplacian(W)          # 2. compute graph laplacian
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  # we will use k-means to cluster the data
  # using the leading eigenvalues in absolute values
  ei$vectors <- ei$vectors[,base::order(abs(ei$values),decreasing=TRUE)]
  obj <- kmeans(ei$vectors[, 2:n_eig], centers = n_eig, nstart = 100)
  if (n_eig==2){
    cl <- 2*(obj$cluster - 1) - 1 
  } else {
    cl <- obj$cluster
  }
  # return the cluster membership
  return(cl) 
}

#' Function to get cluster assignments from applying different methods
#' 
#' @param x An n by p data matrix of compositions
#' @param K The number of clusters
#' @return A list with the following components:
#' \itemize{
#' \item coat - The cluster assignment from hierarchical clustering of the coat correlation matrix.
#' \item cosine - The cluster assignment from spectral clustering of the cosine dissimilarity matrix.
#' \item propr - The cluster assignment from hierarchical clustering of the proportionality matrix. 
#' }
get.clusters <- function(x,K){
  
  cl <- NULL
  fit_coat <- coat(x) # coat is a method that estimates the correlation between compositional variables.
  h <- hclust(as.dist(1-fit_coat$corr),method = 'average')
  cl$coat <- cutree(h,k=K)
  # adjusted.rand.index(cl,cl0)
  
  cosSIM <- lsa::cosine(x)
  cl$cosine <- spectral.clustering(cosSIM,K)
  # adjusted.rand.index(cl,cl0)
  
  # for each cluster, construct SLRs
  pr <- propr::propr(x, metric = "phs")
  h <- hclust(as.dist(pr@matrix))
  cl$propr <- cutree(h,k=K)
  # adjusted.rand.index(cl,cl0)
  
  return(cl)
}

#' Function to generate balances and summed log-ratios from a cluster assignment.
#' 
#' @param x An n by p data matrix of compositions.
#' @param dist_mat The p by p distance matrix.
#' @param cl The cluster assignment obtained a priori.
#' 
#' @return A list with the following components:
#' \itemize{
#' \item sbp - The sequential binary partition matrix.
#' \item slr - The summed log-ratio transform corresponding to the SBP.
#' \item ba - The isometric log-ratio (aka balances) corresponding to the SBP. 
#' }
generate.sbp.for.clusters <- function(x, dist_mat, cl){
  sbp <- NULL
  K <- length(unique(cl))
  for (j in 1:K){
    # check length of cl
    if (sum(cl==j)==2){
      sbp[[j]] <- matrix(c(1,-1),ncol=1)
      rownames(sbp[[j]]) <- rownames(dist_mat[cl==j,cl==j])
    } else {
      h <- hclust(as.dist(dist_mat[cl==j,cl==j]))
      phyloa <- ape::as.phylo(h)
      sbp[[j]] <- philr::phylo2sbp(phyloa)
    }
  }
  tmp_sbp <- as.matrix(bdiag(sbp))
  rownames(tmp_sbp) <- unlist(lapply(sbp,rownames))
  sbp <- tmp_sbp[match(paste0('s',1:p),rownames(tmp_sbp)),]
  colnames(sbp) <- paste0("z", 1:ncol(sbp))
  slr <- slr.fromSBP(x,sbp)
  ba <- balance.fromSBP(x,sbp)
  
  return(list(sbp=sbp,slr=slr,ba=ba))
}

#' Function to generate an exponential decay covariance matrix. 
#' The inverse of this matrix is sparse and defines the AR(1) graph.
#' 
#' @param p  The number of variables.
#' @param rho The correlation (non-negative).
#' @return A list with the following components:
#' \itemize{
#' \item Sigma - The covariance matrix.
#' \item Omega - The inverse covariance matrix. 
#' \item Adj -The 0-1 adjacency matrix.
#' }
rgExpDecay <- function(
  p,  
  rho, 
  const = 0.01
){
  Sigma = matrix(rho, p, p)
  index = matrix(rep(1:p, p), p, p) - t(matrix(rep(1:p, p), p, p))
  Sigma = Sigma^(abs(index))
  Omega = solve(Sigma)
  Omega <- (Omega + t(Omega))/2
  Adj <- (abs(Omega)>1e-08) - diag(rep(1,p))
  return(list(Sigma = Sigma, Omega = Omega, Adj = Adj))
}


#' Function to Transform Samples with the ilr of a Contrast
#'
#' @param x A matrix with rows as samples (n) and columns as components (p).
#' @param contrast A vector. One column of a serial binary partition matrix
#'  with values [-1, 0, 1] describing p components.
#'
#' @return A transformation of samples for the contrast provided.
#' 
slr.fromContrast <- function(x, contrast){
  
  if(length(contrast) != ncol(x)) stop("Contrast must have length ncol(x) = p.")
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")
  
  lpos <- sum(contrast == 1)
  lneg <- sum(contrast == -1)
  const <- sqrt((lpos*lneg)/(lpos+lneg))
  
  # logX <- log(x)
  ipos <- rowSums(x[, contrast == 1, drop = FALSE])
  ineg <- rowSums(x[, contrast == -1, drop = FALSE])
  
  log(ipos / ineg)
}

#' Function to compute summed log-ratios from an SBP Matrix
#'
#' @param x A matrix with rows as samples (n) and columns as components (p).
#' @param y A serial binary partition matrix with rows as components (p) and
#'  columns as balances (p-1).

#' @return A transformation of samples for each contrast in the SBP matrix.
#' 
slr.fromSBP <- function(x, y){
  
  if(!identical(colnames(x), rownames(y))){
    
    stop("Component names for data matrix and balance matrix do not match.")
  }
  
  x <- as.matrix(x)
  
  if(any(x == 0)){
    
    message("Alert: Replacing 0s with next smallest value to calculate balances.")
    zeros <- x == 0
    x[zeros] <- min(x[!zeros])
  }
  
  res <- apply(y, 2, function(z) slr.fromContrast(x, z))
  rownames(res) <- as.character(1:nrow(res))
  return(res)
}


#' A wrapper function to fit Lasso at cross validation selected tuning parameters.
#' 
#' @param x The n by p predictors in the training set
#' @param y The response variable in the training set
#' @param xt The predictors in the test set
#' @param yt The reponse variable in the test set
#' 
run.glmnet <- function(x,y,xt,yt,lambda=NULL){
  cv_exact <- cv.glmnet(x=x,y=y,lambda=lambda)
  lambda <- log(cv_exact$lambda)
  lambda_new <- exp(seq(max(lambda),min(lambda)+2,length.out = 100))
  cv_exact <- cv.glmnet(x=x,y=y,lambda=lambda_new)
  refit_exact <- glmnet(x=x,y=y,family='gaussian',lambda=cv_exact$lambda.min)
  pred_exact <- predict(cv_exact, newx = xt, type = "response", s = 'lambda.min')
  
  return(list(beta.min=refit_exact$beta,
              beta=cv_exact$glmnet.fit$beta,
              lambda=cv_exact$lambda,
              mse.pred = mean((pred_exact-yt)^2)))
}

#' Function to evaluate the performance of selbal
#' 
#' @param x The n by p predictors in the training set
#' @param y The response variable in the training set
#' @param xt The predictors in the test set
#' @param yt The reponse variable in the test set
#' 
pred.from.selbal <- function(x,y,xt,yt){
  fit <- selbal::selbal(x,y,draw=F)
  design <- cbind(rep(1,nrow(x)), fit[[1]])
  bal_lm <- lm(y~design)
  
  contrast_selbal <- matrix(0,ncol(x),1)
  rownames(contrast_selbal) <- colnames(x)
  contrast_selbal[match(fit[[2]],rownames(contrast_selbal)),1] <- 1
  contrast_selbal[match(fit[[3]],rownames(contrast_selbal)),1] <- -1
  design[,2] <- balance.fromSBP(xt,contrast_selbal)
  mse.pred <- mean((yt - predict(bal_lm, newx=design))^2)
  
  return(list(mse.pred=mse.pred,pos=fit[[2]],neg=fit[[3]],coeff=bal_lm$coefficients))
} 


#' Function to evaluate the true positive rate of original variables.
#' 
#' @param beta The true coefficient vector.
#' @param beta_hat The estimated coefficient vector.
#' @param eps A threshold for declaring a non-zero discovery. Default is 1e-08.
#' @return The total number of discovery and true positive rate.
#' 
tpr.for.coef <- function(beta, beta_hat, eps = 1e-08){
  TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
  FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
  tpr <- TP/(TP + FN)
  S_hat <- sum((abs(beta_hat) > eps))
  out <- c(S_hat,tpr)
  names(out) <- c('S_hat','tpr')
  return(out)
}

#' Function to evaluate the true positive rate of the ilr transformed variables.
#' 
#' @param beta The true coefficient vector.
#' @param beta_hat The estimated coefficient vector.
#' @param sbp The sequential binary partition matrix. This matrix defines the balances.
#' @param eps A threshold for declaring a non-zero discovery. Default is 1e-08.
#' @return The total number of discovery and true positive rate.
tpr.for.coef.ilr <- function(beta,beta_hat,sbp,eps=1e-08){
  
  if (is.null(sbp)){
    stop('A sequential binary partition tree is needed for tpr evaluation!')
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
  tpr <- length(TP)/length(S0)
  out <- c(length(S_hat),tpr)
  names(out) <- c('S_hat','tpr')
  return(out)
}

