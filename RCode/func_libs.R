# library(gdata)
# library(glmnet)
library(Rcpp)
# library(dirmult)
# library(MGLM)
library(MASS)
library(compositions)
library(matrixStats)
# # source("COAT-master/coat.R")
# # solves t(beta)%*%xx%*%beta/2 - t(xy)%*%beta + fac||beta||_1 
# #or ||X%*%bet-y||_2^2/2 + fac||beta||_1 with xx=crossprod(X), xy=crossprod(x,y)
# # subject to t(cmat)%*%beta=0
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



log_var_correction <- function(W, method="ZR", ZRvalue=0.5, overdisp=NULL, ID=NULL){
  # a function to perform log transformation of error-in-variable covariates
  if (method=="ZR"){
    W[W==0] <- ZRvalue
    W <- log(W)
  }else if (method=="AH"){
    W <- W+0.5
    W <- log(W)
  }else if (method=="MOM"){
    ID_unique <- unique(ID)
    overdisp <- rep(1, nrow(W))
    for (i in ID_unique){
      overdisp[ID==i] <- est_overdisp_mom(W[ID==i,])
    }
    # overdisp <- est_overdisp_mom(W)
    W <- W+overdisp/2
    W <- log(W)
  }else if (method=="custom"){
    if (!is.null(overdisp)){
      W <- W+overdisp/2
      W <- log(W)
    }else{
      stop('Customized overdispersion is not provided!')
    }
  }else{
    stop('The input for method is not accepted!')
  }
  # W <- W-rowMeans(W)
  return(W)
}


est_overdisp_mom <- function(W){
  if (is.null(nrow(W))){
    W <- matrix(W,ncol=p)
  }
  # estimate the overdispersion parameter for samples from the same subject using moments method
  # W is n by p
  # overdisp is a vector of length n
  n <- nrow(W)
  p <- ncol(W)
  Ni <- rowSums(W)
  N <- sum(Ni)
  Pij <- W/Ni
  Pj <- colSums(W)/N
  Nc <- (N-sum(Ni^2)/N)/(n-1)
  Sj <- colMeans((scale(Pij, scale=F, center=Pj)^2)*Ni)
  Gj <- colSums(W*(1-Pij))/sum(Ni-1)
  lambdai <- sum(Sj+(Nc-1)*Gj)/sum(Sj-Gj) # the larger it is, the smaller variance is
  overdisp <- (lambdai+Ni)/(1+lambdai)
  print(lambdai)
  return(overdisp)
}

## ---- Part II ----

#' @description Random number generator for logistic normal distribution.
#' @param n sample size
#' @param mu mean vector
#' @param Sigma covariance matrix 
random_logistic_normal <- function(n,mu,Sigma){
  d <- MASS::mvrnorm(n=n,mu=mu,Sigma=Sigma)
  d <- matrix(d, nrow=n)
  d <- cbind(d,rep(0,n))
  d <- exp(d)
  compositions::acomp(d)
}

#' @param n sample size
#' @param p number of compositional variables - 1
#' @param q number of confounding factors
#' @param sig standard deviation of the noise added to response
data.sims <- function(n,p,q,sig,beta0,seed.id=1){
  
  set.seed(seed.id)
  # Generate confounding factors from normal distribution
  H <- matrix(rnorm(n*q),n,q)
  delta <- rnorm(q)
  
  # Generate coefficients Gamma and errors E from logistic normal distribution
  # Gamma is a matrix of dim q x (p+1), each row is a composition
  # E is a matrix of dim n x (p+1), each row is a composition
  Gamma <- random_logistic_normal(n=q,mu=rep(0,p),Sigma=diag(1,p))  
  E <- random_logistic_normal(n=n,mu=rep(0,p),Sigma=diag(1,p))  
  
  # Generate the predictor X
  # X is a matrix of dim n x (p+1), each row is a composition
  if (q>1){
    X <- acomp(t(sapply(1:n,function(i) acomp(colProds(power.acomp(Gamma,H[i,])))))) + E
  } else if (q==1){
    X <- acomp(t(sapply(1:n,function(i) acomp(power.acomp(Gamma,H[i,]))))) + E
  }
  x <- as(X,'matrix')
  
  ## -- Log contrast model --
  y <- log(x) %*% beta0 + H %*% delta + rnorm(n)*sig
  
  return(list(x=x,y=y))
}

#' @description Simulate data where the hidden confounder is independent of X
data.sims.noHX <- function(n,p,q,sig,beta0,seed.id=1){
  
  set.seed(seed.id)
  # Generate confounding factors from normal distribution
  H <- matrix(rnorm(n*q),n,q)
  delta <- rnorm(q)
  
  # Generate coefficients Gamma and errors E from logistic normal distribution
  # Gamma is a matrix of dim q x (p+1), each row is a composition
  # Gamma <- random_logistic_normal(n=q,mu=rep(0,p),Sigma=diag(1,p))  
  # Generate the predictor X
  X <- random_logistic_normal(n=n,mu=rep(0,p),Sigma=diag(1,p))  
  x <- as(X,'matrix')
  
  ## -- Log contrast model --
  y <- log(x) %*% beta0 + H %*% delta + rnorm(n)*sig
  
  return(list(x=x,y=y))
}

## Apply Trim transformation to Predictors and response
#' @param x Input data matrix of dimension nobs x nvars, compositional predictor
#' @param y A length n vector of response
# Select tuning parameter by minimizing the L1 estimation error of beta, given the true beta
fit_lc_confounding <- function(x,y,beta0){
  p <- ncol(x)
  n <- nrow(x)
  errorMat <- matrix(NA,nrow=2,ncol=3)
  rownames(errorMat) <- c('pre','post')
  colnames(errorMat) <- c('pred','l1estimation','l2estimation')
  lambda <- exp(seq(from=log(2), to=log(1e-4), length.out=50))
  
  # Run the constrained lasso without correcting the predictor
  C <- matrix(1,p,1)
  Z <- log(x)
  
  # Apply transformation
  Pc <- tcrossprod(C)/p
  Z_clr <- Z-Z%*%Pc
  z.svd <- svd(Z_clr)
  tau <- median(z.svd$d)
  Dt <- sapply(z.svd$d, function(a) ifelse(a>tau,tau,a))
  FM <- z.svd$u %*% diag(Dt/z.svd$d) %*% t(z.svd$u)
  Zt <- FM %*% Z_clr
  yt <- FM %*% y
  
  maxlam <- 2*max(abs(crossprod(Zt,yt)/n))
  lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=50))
  
  # res_post <- cv.func('ConstrLasso', x=Zt, y=yt, Cmat=C, nlam=50)
  cat('Fitting the corrected model ...\n')
  post <- ConstrLasso(yt,Zt,Cmat=C,lambda=lambda)
  tmp <- apply(post$bet, 2, function(a) sum(abs(a-beta0)))
  pred_classo_post <- Zt %*% post$bet[,which.min(tmp)]
  errorMat[2,1] <- mean((pred_classo_post - yt)^2)
  errorMat[2,2] <- min(tmp)
  errorMat[2,3] <- sum((post$bet[,which.min(tmp)] - beta0)^2)
  
  cat('Fitting the original model ...\n')
  pre <- ConstrLasso(y,x=Z,Cmat=C,lambda=lambda)
  tmp <- apply(pre$bet, 2, function(a) sum(abs(a-beta0)))
  pred_classo_pre <- Z %*% pre$bet[,which.min(tmp)]
  errorMat[1,1] <- mean((pred_classo_pre - y)^2)
  errorMat[1,2] <- min(tmp)
  errorMat[1,3] <- sum((pre$bet[,which.min(tmp)] - beta0)^2)
  
  return(errorMat)
}

## ---- Part III ----
gm.mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

clr <- function(x, base=exp(1)){
  x <- log((x / gm.mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}
#' Function to double center a given matrix
#' @param M a matrix of dimension n by n
#' @return a double centered matrix
GowerCentering <- function(M){
  d = dim(M); n = d[1]
  II = matrix(1,n,n)
  Mgc = (diag(n)-(1/n)*II)%*%M%*%(diag(n)-(1/n)*II)
}

# Function to compute the population correlation matrix
popGammajk = function(
  alpha1j, alpha1k, beta1, var_epsilon, var_epsilonj, var_epsilonk, U){
  varU = stats::var(U) #(1 / 12) * (0.5 - (-0.5))
  corrjk = ((alpha1j - alpha1k) * beta1 * varU) / 
    sqrt((beta1^2 * varU + var_epsilon) * 
           ((alpha1j - alpha1k)^2 * varU + var_epsilonj + var_epsilonk))
  return(abs(corrjk))
}
popGamma = function(
  alpha1, beta1, var_epsilon, var_epsilon2, U
){
  p = length(alpha1)
  if(length(var_epsilon2) == 1) var_epsilon2 = rep(var_epsilon2, p)
  
  rhoMat = matrix(0, p, p)
  for (j in 1:p){
    for (k in 1:p){
      if (k==j){next}
      else {
        rhoMat[j, k] = popGammajk(
          alpha1j = alpha1[j], alpha1k = alpha1[k], beta1 = beta1, 
          var_epsilon = var_epsilon, var_epsilonj = var_epsilon2[j], 
          var_epsilonk = var_epsilon2[k], U = U)
      }
    }
  }
  return(rhoMat)
}

# Function to compute the graph Laplacian
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

#' @param n_eig Number of clusters/eigenvectors to obtain
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

spectral.clustering.sign <- function(W, n_eig = 2, zeta=0, normalized=TRUE) {
  stopifnot(nrow(W) == ncol(W)) 
  
  n = nrow(W)    # number of vertices
  # We perturb the network by adding some links with low edge weights
  W <- W + zeta * mean(colSums(W))/n * tcrossprod(rep(1,n))
  g <- colSums(W) # degrees of vertices
  
  if(normalized){
    D_half = diag(1 / sqrt(g) )
    L = diag(n) - D_half %*% W %*% D_half
  } else {
    L = diag(g) - W 
  }
  ei = eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors and values of L
  # we will use sign of the eigenvectors corresponding to the 2nd smallest eigenvalue
  ei$vectors <- ei$vectors[,base::order(ei$values,decreasing=FALSE)]
  cl <- 2*(ei$vectors[, 2]>0)-1
  names(cl) <- rownames(W)
  # return the cluster membership
  return(cl) 
}

#' @param x A matrix with rows as samples and columns as components
#' @param method Method to evaluate the similairity between components
#' @param threshold Threshold value to determine whether two components are similar to each other. If threshold is passed, two components are deemed to 
#' be similar to each other. 
sbp.fromRSC <- function(A){
  target <- rep(0,ncol(A))
  names(target) <- colnames(A)
  
  # recursive function for binary partition with spectral clustering
  BinPartRSC = function(A,target){
    # cat('number of rows in A: ', rownames(A), '..\n')
    # ret$index <- rownames(A)
    bp <- rep(0,length(target))
    if (nrow(A) == 2) {
      bp[names(target) %in% rownames(A)] <- c(1,-1) # not necessarily have 2 elements
      return(bp)
    }
    if (nrow(A) == 1 || sum(A)==0) {# empty network
      bp[names(target) %in% rownames(A)] <- 1 # not necessarily have 2 elements
      return(bp)
    }
    
    # recursive partition
    obj <- spectral.clustering(A)
    bp[names(target) %in% rownames(A)] <- obj
    
    if (length(unique(obj))==2){
      ind1 <- which(obj==1)
      A1 <- matrix(A[ind1,ind1],nrow=length(ind1))
      rownames(A1) <- colnames(A1) <- rownames(A)[ind1]
      ind2 <- which(obj==-1)
      A2 <- matrix(A[ind2,ind2],nrow=length(ind2))
      rownames(A2) <- colnames(A2) <- rownames(A)[ind2]
      return(cbind(bp,BinPartRSC(A1,target),BinPartRSC(A2,target)))
    }
  }
  # Gather all binary partitions
  sbp <- BinPartRSC(A,target)
  # print(dim(sbp))
  if (min(colSums(abs(sbp)))==1){
    sbp <- sbp[,-which(colSums(abs(sbp))==1)]
  }
  # check if there is any column that has only one unique cluster. They will be removed for later analysis
  test_unique <- apply(sbp,2, function(a) c(1,-1) %in% unique(a))
  if (min(colSums(test_unique))==1){
    sbp <- sbp[,-which(colSums(test_unique)==1)]
    sbp <- matrix(sbp,nrow=nrow(A))
  }
  b.weight <- apply(sbp, 2, function(i) sum(abs(i)))
  b.order <- order(b.weight, decreasing = TRUE)
  sbp <- sbp[,b.order]
  sbp <- matrix(sbp,nrow=nrow(A))
  colnames(sbp) <- paste0('z',1:ncol(sbp))
  rownames(sbp) <- rownames(A)
  
  sbp
}

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

#' @param x An n by p data matrix of compositions
#' @return A set of clusterings with different distance measures
get.clusters <- function(x,K){
  
  cl <- NULL
  fit_coat <- coat(x)
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

#' @param x n by p compositional data matrix
#' @param dist_mat The p by p distance matrix
#' @param cl The cluster membership obtained a priori.
#' @return A sequential binary partition matrix and the summed log-ratio corresponding to the SBP. 
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

#' @description This function uses regularized spectral clustering instead of hierarchical clustering
#' with average linkage to get the SBP for each cluster. 
#' @param A p by p similarity graph defined a priori.
#' 
generate.sbp.for.clusters.rsc <- function(x, A, cl){
  sbp <- NULL
  K <- length(unique(cl))
  for (j in 1:K){
    if (sum(cl==j)==2){
      sbp[[j]] <- matrix(c(1,-1),ncol=1)
      rownames(sbp[[j]]) <- rownames(A[cl==j,cl==j])
    } else {
      sbp[[j]] <- sbp.fromRSC(A[cl==j,cl==j])
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
rgExpDecay <- function(
  p,  # number of variables
  rho, # the decaying correlation
  const = 0.01
){
  Sigma = matrix(rho, p, p)
  index = matrix(rep(1:p, p), p, p) - t(matrix(rep(1:p, p), p, p))
  Sigma = Sigma^(abs(index))
  Omega = solve(Sigma)
  Omega <- (Omega + t(Omega))/2
  Adj <- (abs(Omega)>1e-08) - diag(rep(1,p))
  # weights <- matrix(rho, p, p)
  # Omega <- Adj*weights
  # diag(Omega) <- ifelse(min(eigen(Omega)$val) > 0, 0, - min(eigen(Omega)$val)) + const
  # Sigma <- solve(Omega)
  # Sigma <- cov2cor(Sigma)
  # Omega <- solve(Sigma)
  return(list(Sigma = Sigma, Omega = Omega, Adj = Adj))
}
#' Generate random scale free graphs from the Barabasi-Albert model
#' @param n The number of vertices in the graph.
#' @param m Numeric constant, the number of edges to add in each time step. 
rg.scale.free.BA <- function(n, m=1, seed=1,const=0.01){
  set.seed(seed)
  g <- igraph::sample_pa(n,m=m,directed = FALSE)
  A <- as.matrix(as_adjacency_matrix(g,type = "both"))
  weights <- matrix(0, n, n)
  upperTriangle(weights, diag = F) <- runif((n*(n - 1))/2, 0.5, 1)*(2*rbinom((n*(n - 1))/2, 1, 0.5) - 1)
  weights <- weights + t(weights)
  Omega <- A*weights
  diag(Omega) <- ifelse(min(eigen(Omega)$val) > 0, 0, - min(eigen(Omega)$val)) + const
  Sigma <- solve(Omega)
  Sigma <- cov2cor(Sigma)
  Omega <- solve(Sigma)
  return(list(Omega = Omega, Sigma = Sigma, Adj=A))
}

#' Generate random graphs from the ER model
#' @param n The number of vertices in the graph.
#' @param m The number of edges in the graph	
rg.ER <- function(n, m, const=0.01){
  set.seed(1)
  g <- sample_gnm(n=n,m=m,directed = FALSE)
  A <- as.matrix(as_adjacency_matrix(g,type = "both"))
  weights <- matrix(0, n, n)
  upperTriangle(weights, diag = F) <- runif((n*(n - 1))/2, 0.5, 1)*(2*rbinom((n*(n - 1))/2, 1, 0.5) - 1)
  weights <- weights + t(weights)
  Omega <- A*weights
  diag(Omega) <- ifelse(min(eigen(Omega)$val) > 0, 0, - min(eigen(Omega)$val)) + const
  Sigma <- solve(Omega)
  Sigma <- cov2cor(Sigma)
  Omega <- solve(Sigma)
  return(list(Omega = Omega, Sigma = Sigma, Adj=A))
}

# rg.sbm <- function(n, pref.matrix, block.sizes, const=0.01){
#   g <- sample_sbm(n,pref.matrix, block.sizes, directed = FALSE)
#   A <- as.matrix(as_adjacency_matrix(g,type = "both"))
#   weights <- matrix(0, n, n)
#   upperTriangle(weights, diag = F) <- runif((n*(n - 1))/2, 0.5, 1)*(2*rbinom((n*(n - 1))/2, 1, 0.5) - 1)
#   weights <- weights + t(weights)
#   Omega <- A*weights
#   diag(Omega) <- ifelse(min(eigen(Omega)$val) > 0, 0, - min(eigen(Omega)$val)) + const
#   Sigma <- solve(Omega)
#   Sigma <- cov2cor(Sigma)
#   Omega <- solve(Sigma)
#   return(list(Omega = Omega, Sigma = Sigma, Adj=A))
# }

#' #' @description Simulation compositions from a log-normal basis distribution.
#' modelExpDecay <- function(
#'   p,  # number of variables
#'   rho, # the decaying correlation
#'   const = 0.01
#' ){
#'   Sigma = matrix(rho, p, p)
#'   index = matrix(rep(1:p, p), p, p) - t(matrix(rep(1:p, p), p, p))
#'   Sigma = Sigma^(abs(index))
#'   Omega = solve(Sigma)
#'   Adj <- (abs(Omega)>1e-08) - diag(rep(1,p))
#'   return(list(Sigma = Sigma, Omega = Omega, Adj = Adj))
#' }

#' #' @description This function calculates the jaccard similarity between rows of A
#' get_jaccard <- function(A){
#'   sim.jac <- matrix(0, nrow=nrow(A), ncol=nrow(A))
#'   rownames(sim.jac) <- rownames(A)
#'   colnames(sim.jac) <- rownames(A)
#'   
#'   #weighted jaccard
#'   pairs <- t(combn(1:nrow(A), 2))
#'   for (i in 1:nrow(pairs)){
#'     num <- sum(sapply(1:ncol(A), function(x)(min(A[pairs[i,1],x],A[pairs[i,2],x]))))
#'     den <- sum(sapply(1:ncol(A), function(x)(max(A[pairs[i,1],x],A[pairs[i,2],x]))))
#'     sim.jac[pairs[i,1],pairs[i,2]] <- num/den
#'     sim.jac[pairs[i,2],pairs[i,1]] <- num/den  
#'   }
#'   sim.jac[which(is.na(sim.jac))] <- 0
#'   diag(sim.jac) <- 1
#'   return(sim.jac)
#' }

#' Computate a supervised partition tree
#' @param X relative abundances
#' @param y covariate 
#' @param linkage
getSupervisedTree = function(y, X, linkage = "average", rho.type = "square"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # checks
  if(length(y) != n) stop("getSupervisedTree() error: dim(X)[1] != length(y)!")
  
  # calculate correlation of each pair of log-ratios with response y
  cormat = matrix(0, p, p) # diagonal == 1
  y_demeaned = y - mean(y)
  for (j in 1:(p - 1)){
    for (k in (j + 1):p){
      Zjk = log(X[, j]) - log(X[, k])
      Zjk_demeaned = Zjk - mean(Zjk)
      if(rho.type == "square" | rho.type == "squared" | rho.type == "s" | 
         rho.type == 2){
        # val = (cor(Zjk_demeaned, y_demeaned))^2
        val = (cor(Zjk, y))^2
      } else{
        # val = abs(cor(Zjk_demeaned, y_demeaned))
        val = abs(cor(Zjk, y))
      }
      # if(is.na(val)) stop("getSupervisedTree() : correlation = 0")
      # hopefully we never have to use this line below.
      # if(is.na(val)) val = 1 #######################################################################################
      cormat[j, k] = val
      cormat[k, j] = val
    }
  }
  # give the rows and columns the names of taxa in X, for sbp.fromHclust()
  rownames(cormat) = colnames(X)
  rownames(cormat) = colnames(X)
  
  # get dissimilarity matrix
  Gamma = 1 - cormat
  
  # get tree from hierarchical clustering
  btree_slr = hclust(as.dist(Gamma), method = linkage)
  
  return(btree_slr)
}
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

#' Compute Balances from an SBP Matrix
#'
#' @param x A matrix with rows as samples (N) and columns as components (D).
#' @param y A serial binary partition matrix with rows as components (D) and
#'  columns as balances (D-1).
#'
#' @return A transformation of samples for each balance in the SBP matrix.
#'
#' @author Thom Quinn
#'
#' @examples
#' library(balance)
#' data(iris)
#' x <- iris[,1:4]
#' sbp <- sbp.fromPBA(x)
#' balance.fromSBP(x, sbp)
#'
#' @export
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


require(glasso)
adjDGlasso <- function(
  X, #the n by p data matrix
  weights=1, #the weight for the penalty
  theta_star=NULL, #true precision matrix
  lambda = NULL, 
  FDR.type='BH', #FDR control procedure
  alpha = 0.1, #the significance level for FDR control
  quiet=TRUE
){
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, center = T, scale = F)
  empcov <- (1/n) * (t(X) %*% X) #empirical cov
  if (is.null(lambda)){
    lambda <- sqrt(log(p)/n)
  }
  
  if (!quiet){print('fit glasso')}
  Theta.hat.from.Glasso <- glasso(s=empcov, rho=lambda*weights, penalize.diagonal=FALSE)$wi
  if (!quiet){print('done')}
  
  if (!quiet){print('de-biasing glasso')}
  ## T.hat = Theta.hat - Theta.hat * (Sigma.hat - Theta.hat^{-1}) * Theta.hat
  ## Theta.hat and Sigma.hat are both symmetric.
  temp.mat <- empcov - chol2inv(chol(Theta.hat.from.Glasso))
  temp.vec <- as.vector(Theta.hat.from.Glasso %*% temp.mat %*% t(Theta.hat.from.Glasso))
  T.hat <- as.vector(Theta.hat.from.Glasso) - temp.vec
  T.hat <- matrix(T.hat,nrow = p)
  
  if (!quiet){print('done')}
  
  sigma.hat2 <- array(0,c(p,p))
  for (i in 1:p){
    for (j in 1:p){
      sigma.hat2[i,j] <- Theta.hat.from.Glasso[i,j]^2+Theta.hat.from.Glasso[i,i]*Theta.hat.from.Glasso[j,j]
    }
  }
  
  test.stat <- sqrt(n)*T.hat/sqrt(sigma.hat2)
  std.statistic <- NULL
  if (!is.null(theta_star)){
    std.statistic <- sqrt(n)*(T.hat - theta_star)/sqrt(sigma.hat2)
  }
  
  pvals <- 2*(pnorm(abs(test.stat), lower.tail=FALSE))
  pvals.vec <- lowerTriangle(pvals,diag=FALSE)
  adjpvals.vec <- p.adjust(pvals.vec, FDR.type)
  coeff <- diag(1,p) - cov2cor(T.hat)
  # coeff <- - pmax(pmin(T.hat, 1), -1)
  Qmat <- matrix(0, p, p)
  lowerTriangle(Qmat, diag=FALSE) <- adjpvals.vec
  Qmat <- Qmat + t(Qmat)
  diag(Qmat) <- rep(1, p)
  Qmat.fdr <- (Qmat <= alpha)
  
  # Qmat is the p by p matrix of qvalues;
  # Qmat.fdr is the thresholded matrix of qvalues based on alpha
  return(list(Theta=coeff, pvalue=pvals, qvalue = Qmat, qvalue.fdr=Qmat.fdr, statistic=std.statistic))
}


#' A wrapper function to fit Lasso at cross validation selected tuning parameters
#' @param x n by p training data
#' @param y response
#' @param xt test data
run.glmnet <- function(x,y,xt,yt,lambda=NULL){
  cv_exact <- cv.glmnet(x=x,y=y,lambda=lambda)
  # lambda <- log(cv_exact$lambda)
  # lambda_new <- exp(seq(max(lambda),min(lambda)+2,length.out = 100))
  # cv_exact <- cv.glmnet(x=x,y=y,lambda=lambda_new)
  refit_exact <- glmnet(x=x,y=y,family='gaussian',lambda=cv_exact$lambda.min)
  pred_exact <- predict(cv_exact, newx = xt, type = "response", s = 'lambda.min')
  
  return(list(beta.min=refit_exact$beta,
              beta=cv_exact$glmnet.fit$beta,
              lambda=cv_exact$lambda,
              mse.pred = mean((pred_exact-yt)^2)))
}

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

#' @param p Number of nodes
#' @param d0 Number of layers in the tree
#' @param deg Average node degree
sim.BTSBM.cov <- function(p,d0,deg,rho,sid=1){
  set.seed(sid)
  dt <- BTSBM(n=p,d=d0,a.seq=0.2^seq(0,d0),lambda=deg)
  A0 <- as.matrix(dt$A.list[[1]])
  rownames(A0) <- colnames(A0) <- paste0('s',1:p)
  Sigma <- rho*A0
  diag(Sigma) <- ifelse(min(eigen(Sigma)$val) > 0, 0, - min(eigen(Sigma)$val)) + 0.1
  Sigma <- cov2cor(Sigma)
  rownames(Sigma) <- colnames(Sigma) <- paste0('s',1:p)
  # sbp <- sbp.fromRSC(A0)
  HCD.result <- HCD(A0,method="SC", stopping="Fix", D=d0, n.min = round(deg/2),notree = FALSE)
  # phyloa <- ape::as.phylo(HCD.result$cluster.tree)
  # sbp <- philr::phylo2sbp(phyloa)
  # if (min(colSums(abs(sbp)))==2){
  #   sbp <- sbp[,-which(colSums(abs(sbp))==2)]
  # }
  # colnames(sbp) <- paste0('z',1:ncol(sbp))
  # rownames(sbp) <- paste0('s',rownames(sbp))
  # sbp <- sbp[match(paste0('s',1:nrow(A0)),rownames(sbp)),]
  
  return(list(S=Sigma, A=A0, HCD.result=HCD.result))
}

get.HCD.tree <- function(A,d0,deg){
  HCD.result <- HCD(A,method="SC",stopping="Fix",D = d0, n.min = deg,notree = FALSE)
  phyloa <- ape::as.phylo(HCD.result$cluster.tree)
  sbp <- philr::phylo2sbp(phyloa)
  if (min(colSums(abs(sbp)))==2){
    sbp <- sbp[,-which(colSums(abs(sbp))==2)]
  }
  colnames(sbp) <- paste0('z',1:ncol(sbp))
  rownames(sbp) <- paste0('s',rownames(sbp))
  sbp <- sbp[match(paste0('s',1:nrow(A)),rownames(sbp)),]
  
  return(list(sbp=sbp, HCD.result=HCD.result))
}
# image.plot(A0)
# balance.plot(WC,sbp)
# plot(HCD.result$cluster.tree)

#' This function evaluates the true positive rate of beta_hat
roc.for.coef <- function(beta_hat, beta, eps = 1e-08){
  TP = sum((abs(beta_hat) > eps) * (abs(beta) > eps))
  FN = sum((abs(beta_hat) <= eps) * (abs(beta) > eps))
  tpr <- TP/(TP + FN)
  S_hat <- sum((abs(beta_hat) > eps))
  out <- c(S_hat,tpr)
  names(out) <- c('S_hat','tpr')
  return(out)
}

#' @description An additional parameter sbp is needed in order to evaluate the variable selection
#' performance
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
  tpr <- length(TP)/length(S0)
  out <- c(length(S_hat),tpr)
  names(out) <- c('S_hat','tpr')
  return(out)
}

