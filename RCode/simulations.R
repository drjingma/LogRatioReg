rm(list=ls())

## Last update: 2022-05-13
today <- '20220513'

jid <- 2
set.seed(jid)
# devtools::install_github(repo = "malucalle/selbal")
# library(selbal)
# library(philr)
# library(propr)
# library(balance)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(ape)
library(igraph)
library(Matrix)
library(fields)
library(pheatmap)
library(ggplot2)

library(DiagrammeR)

CodeFolder <- "RCode/"
source(paste0(CodeFolder,"func_libs.R"))
source("RCode/SLR.R")

## Generate parameters
rho <- 0   # magnitude of the correlation in cov Sigma
p <- 50    # number of variables 
n <- 100   # number of samples
sigy <- 0.5   # noise level in regression
sigx <- 0.5   # noise level in covariates
r <- 3
s <- 3


# Generate x ----
U <- runif(2*n,-0.5,0.5)
alpha <- 10*c(rep(1,r)/r,-rep(1,s)/s,rep(0,p-1-r-s))
Sigma <- rgExpDecay(p-1,rho)$Sigma
rownames(Sigma) <- colnames(Sigma) <- paste0('s',1:(p-1))
xALR <- tcrossprod(as.matrix(U), as.matrix(alpha)) + mvtnorm::rmvnorm(2*n,mean=rep(0,p-1),sigma = Sigma) * sigx
x <- alrinv(xALR)
colnames(x) <- paste0('s',1:p)

## SVD on alr(x) and clr(x) ----
plot(svd(xALR)$d) # clear that there is one latent variable
plot(svd(xALR)$v[,1])
clrx <- t(apply(x,1,clr))
plot(svd(clrx)$d) # also clear that there is one latent variable
plot(svd(clrx)$v[,1])

## Generate y ----
theta <- 0.25
y <- theta * U + rnorm(2*n) * sigy
yb <- 1*(y>0)
plot(cor(clrx,as.matrix(y)))

## Clustering with the population similarity ----
popRho <- popGamma(c(alpha,0),theta,sigy^2,sigx^2,U) 
# between cluster values are not zero
image.plot(popRho)
image.plot(max(popRho)-popRho)

W = max(popRho)-popRho
W.eig <- eigen(W)
plot(W.eig$values)
L <- graph.laplacian(W,zeta=0.1)
L.eig <- eigen(L)
df <- data.frame(`1`=L.eig$vectors[,1],
                 `2`=L.eig$vectors[,2],
                 `3`=L.eig$vectors[,3]) %>%
  reshape2::melt() %>%
  ggplot(aes(x=rep(1:p,3),y=value)) + geom_point() +
  facet_wrap(~variable)
print(df) 

## Get the empirical similarity matrix ----
## psiMat is the regression coefficient matrix for continuous response and is equivalent to correlations due to scaling of y
## rhoMat is the correlation matrix for binary response
rhoMat <- matrix(0,p,p)
rownames(rhoMat) <- colnames(rhoMat) <- colnames(x)
psiMat <- rhoMat
AitchisonVar <- rhoMat
for (j in 1:p){
  for (k in 1:p){
    if (k==j){next}
    else {
      df <- data.frame(y=y,z=log(x[,j])-log(x[,k]))
      df$z <- scale(df$z)
      df$y <- df$y - mean(df$y)
      # df <- as.data.frame(scale(df))
      psiMat[j,k] <- abs(sum(df$y*df$z)/sum(df$z*df$z))
      rhoMat[j,k] <- abs(stats::cor(log(x[,j])-log(x[,k]),yb))
      AitchisonVar[j,k] <- var(log(x[,j])-log(x[,k]))
    }
  }
}


pheatmap(AitchisonVar,cluster_rows = F, cluster_cols = F,
         main='Heatmap of the Aitchison variation')

# Clustering with the empirical similarity ----
What.psi <- max(psiMat) - psiMat
pheatmap(What.psi,cluster_rows = F, cluster_cols = F,
         main='Heatmap of the regression coefficients')

Lhat.psi <- graph.laplacian(What.psi,zeta=0.1)
Lhat.psi.eig <- eigen(Lhat.psi)
plot(Lhat.psi.eig$values)

df.psi <- data.frame('1'=Lhat.psi.eig$vectors[,1],
                     '2'=Lhat.psi.eig$vectors[,2],
                     '3'=Lhat.psi.eig$vectors[,3]) %>%
  reshape2::melt() %>%
  ggplot(aes(x=rep(1:p,3),y=value)) + geom_point() +
  facet_wrap(~variable)
print(df.psi)
kmeans(Lhat.psi.eig$vectors[,1:3],centers=3,nstart=10)$cluster


What <- max(rhoMat) - rhoMat
pheatmap(What,show_colnames = T, show_rownames = T,
         cluster_rows = F, cluster_cols = F,
         main = 'Heatmap of the original correlation')

What.eig <- eigen(What)
Lhat <- graph.laplacian(What,zeta=0.1)
Lhat.eig <- eigen(Lhat)
plot(Lhat.eig$values)

dfhat <- data.frame('1'=Lhat.eig$vectors[,1],
                    '2'=Lhat.eig$vectors[,2],
                    '3'=Lhat.eig$vectors[,3]) %>%
  reshape2::melt() %>%
  ggplot(aes(x=rep(1:p,3),y=value)) + geom_point() +
  facet_wrap(~variable)
print(dfhat)

cl <- kmeans(Lhat.eig$vectors[,1:3],centers=3,nstart=10)$cluster
# kmeans(Lhat.eig$vectors[,1:2],centers=2,nstart=10)$cluster


## Build balance from these three subsets
SBP <- matrix(0,nrow=length(cl),ncol=3)
SBP[(cl==1),1] <- 1
SBP[(cl==2),1] <- -1
SBP[(cl==1),2] <- 1
SBP[(cl==3),2] <- -1
SBP[(cl==2),3] <- 1
SBP[(cl==3),3] <- -1
rownames(SBP) <- colnames(x)
z <- balance::balance.fromSBP(x,SBP)
cor(z,yb)
apply(abs(SBP),2,sum)


# object <- cv.slr(x,y,method='correlation',response.type = 'continuous', threshold = NULL,s0.perc = 0)
# plot(object$threshold)
# plot(object$cvm)
# object$threshold.1se
# object$threshold.min
# fit <- slr(x,y,method='correlation',response.type = 'continuous',threshold = object$threshold.min,s0.perc = 0)
# fit$sbp

object.wald <- cv.slr(x,y,method='wald',response.type = 'continuous', threshold = NULL,s0.perc = 0)
plot(object.wald$threshold)
plot(object.wald$cvm)
object.wald$threshold.1se
object.wald$threshold.min
fit.wald <- slr(x,y,method='wald',response.type = 'continuous',threshold = object.wald$threshold.1se,s0.perc = 0)
fit.wald$sbp

# object.wald <- cv.slr(x,yb,method='wald',response.type = 'binary', threshold = NULL,s0.perc = 0)
# fit.wald <- slr(x,yb,method='wald',response.type = 'binary', threshold = object.wald$threshold.1se,s0.perc = 0)

