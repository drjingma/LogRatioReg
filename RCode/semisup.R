rm(list=ls())

today <- '20220823'

jid <- 011
set.seed(jid)
# devtools::install_github(repo = "malucalle/selbal")
# library(selbal)
# library(propr)
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
p <- 30    # number of variables 
n <- 100   # number of labeled samples
m <- 10*n   # number of labeled and unlabeled samples
sigy <- 0.25   # noise level in regression
sigx <- 0.25   # noise level in covariates
r <- 3
s <- 3


# Generate x ----
U1 <- runif(m,-0.5,0.5)/2
alpha1 <- 1*c(rep(1,r)/r,-rep(1,s)/s,rep(0,p-1-r-s))
# alpha2 <- 10*c(rep(0,r+s),c(1,1)/2,-c(1,1)/2,rep(0,p-1-r-s-4))
Sigma <- rgExpDecay(p-1,rho)$Sigma
rownames(Sigma) <- colnames(Sigma) <- paste0('s',1:(p-1))
xALR <- tcrossprod(as.matrix(U1), as.matrix(alpha1)) + mvtnorm::rmvnorm(m,mean=rep(0,p-1),sigma = Sigma) * sigx
x <- alrinv(xALR)
colnames(x) <- paste0('s',1:p)


## Generate y ----
theta <- 10
y <- theta * U1 + rnorm(m) * sigy
yb <- 1*(y>0)
clrx <- t(apply(x,1,clr))
plot(cor(clrx,as.matrix(y)))


## Get the empirical similarity matrix ----
## psiMat is the regression coefficient matrix for continuous response and is equivalent to correlations due to scaling of y
## rhoMat is the correlation matrix for binary response
rhoMat <- matrix(0,p,p)
rownames(rhoMat) <- colnames(rhoMat) <- colnames(x)
psiMat <- rhoMat
AitchisonVar <- rhoMat
AitchisonVar.sub <- rhoMat
for (j in 1:p){
  for (k in 1:p){
    if (k==j){next}
    else {
      # df <- data.frame(y=y,z=log(x[,j])-log(x[,k]))
      # df$z <- scale(df$z)
      # df$y <- df$y - mean(df$y)
      # df <- as.data.frame(scale(df))
      # psiMat[j,k] <- abs(sum(df$y*df$z)/sum(df$z*df$z))
      # rhoMat[j,k] <- abs(stats::cor(log(x[,j])-log(x[,k]),yb))
      AitchisonVar[j,k] <- var(log(x[,j])-log(x[,k]))
      AitchisonVar.sub[j,k] <- var(log(x[1:n,j])-log(x[1:n,k]))
    }
  }
}

index <- 1:6
What.aitchison <- AitchisonVar[index,index]
pheatmap(What.aitchison,cluster_rows = F, cluster_cols = F,
         main='Heatmap of the Aitchison similarity')
h.aitchison <- hclust(as.dist(What.aitchison))
plot(h.aitchison)
cat('Full data',cutree(h.aitchison,k=2),'\n')

What.aitchison.sub <- AitchisonVar.sub[index,index]
h.aitchison.sub <- hclust(as.dist(What.aitchison.sub))
plot(h.aitchison.sub)
cat('Labeled data',cutree(h.aitchison.sub,k=2),'\n')

