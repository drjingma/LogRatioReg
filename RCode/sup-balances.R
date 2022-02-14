rm(list=ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(HCD)
# library(phyloseq)
library(philr)
# devtools::install_github(repo = "malucalle/selbal")
# library(selbal)
library(propr)
library(balance)
library(ape)
library(DiagrammeR)
library(igraph)
library(Matrix)
library(microbenchmark)

# CodeFolder <- "../RCode/"
# OutFolder <- "../Output/"
# source("../../Library/censoredGGM.R")
# source("../../Library/Grand_R_Functions.R")
source("RCode/func_libs.R")
source("Kristyn/Functions/HSClust.R")
# set.seed(12)
n <- 100
p <- 20
# G <- diag(p) - tcrossprod(rep(1,p))/p
# U <- matrix(c(0.5,0.5,-0.5,-0.5,
#               sqrt(1/2),-sqrt(1/2),0,0,
#               0,0,sqrt(1/2),-sqrt(1/2)),ncol=3)
# Uinv <- corpcor::pseudoinverse(t(U))

sig <- 0.1
rho <- 0.5
linkage <- 'average'

# Sigma0 <- rgExpDecay(p/4,rho)$Sigma
# Sigma <- as.matrix(Matrix::bdiag(list(Sigma0,Sigma0,Sigma0,Sigma0)))
Sigma <- rgExpDecay(p,rho)$Sigma
rownames(Sigma) <- colnames(Sigma) <- paste0('v',1:p)

## Generate latent variables from a log-normal distribution
## Observed compositions are defined from the latent variables after dividing by total.
mu <- c(rep(log(p),5),rep(0,p-5))
logW <- mvrnorm(n=n*2, mu=mu, Sigma=Sigma)
W <- exp(logW) # basis
colnames(W) <- colnames(Sigma)
rownames(W) <- paste0('s',1:(n*2))
WC <- sweep(W,1,rowSums(W), FUN='/')
WCLR <- t(apply(WC,1,clr))

normalize <- function(contrast){
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")
  
  lpos <- sum(contrast == 1)
  lneg <- sum(contrast == -1)
  const <- sqrt((lpos*lneg)/(lpos+lneg))
  contrast[contrast==1] = 1/lpos
  contrast[contrast==-1] = -1/lneg
  
  const * contrast
}
## Generate coefficient
## Oracle will be fitting an OLS model, rather than a Lasso regression
# Construct the oracle binary tree
h_oracle <- hclust(as.dist(1-abs(Sigma)), method=linkage)
# sbp_oracle <- sbp.fromHclust(h_oracle)
# U <- apply(sbp_oracle,2,normalize) # sbp is a partition matrix 
# ba_oracle <- balance.fromSBP(WC, sbp_oracle)
# theta_oracle <- rep(0,p-1)
# theta_oracle[c(3,4)] <- 2
# beta_oracle <- U %*% theta_oracle
# Construct a regression coefficient so that all variables within positively correlated blocks have the same sign
beta_oracle <- 10*c(rep(1,p/4),
                 rep(-1,p/4), rep(0, p/2))
# beta_oracle <- c(rep(1,5),
#                  rep(-1,5), 
#                  rep(1,2),rep(-1,2),rep(0,3),
#                  rep(0,3))
names(beta_oracle) <- colnames(W)
# beta_lc <- c(1,0,-0.8,0,0.4,0,-1.5,0,1.2,0,-0.3,rep(0,p-11))
# beta_lc <- c(1,0.4,1.2,-1.5,-0.8,-0.3,rep(0,p-6))
# beta_lc <- c(1,-0.8,0.4,0,-0.6,0,-1.5,0,1.2,0,0.3,rep(0,p-11))
beta_support <- which(abs(beta_oracle)>0)

## Generate response
yAll <-  log(WC) %*% beta_oracle + rnorm(2*n) * sig
y <- yAll[1:n,]
y_test <- yAll[-(1:n),]

## ---- Set the error matrix ----
err <- matrix(NA,nrow=5,ncol=1)
rownames(err) <- c('supervised','propr','pba','classo','oracle')
colnames(err) <- c('pred')

## Benchmark the computational time of various methods.
# mbm <- microbenchmark(
#   "supervised" = {
    S <- matrix(0,p,p)
    rownames(S) <- colnames(S) <- colnames(W)
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        newx <- log(WC[1:n,j]) - log(WC[1:n,k])
        newx <- newx - mean(newx)
        S[j,k] <- S[k,j] <- cor(newx,y)^2
      }
    }
    h_supervised <- hclust(as.dist(1-S),method = linkage)
    plot(hclust(as.dist(1-S),method = 'single'))
    # ggtree(h_supervised) + geom_point(aes(shape=isTip, color=isTip), size=3)
    sbp_supervised <- sbp.fromHclust(h_supervised)
    ba_supervised <- balance.fromSBP(WC,sbp_supervised)
    
    a <- HSClust(S)
    sbp_spec <- sbp.fromHSClust(a$allLevels)
#     fit_supervised <- run.glmnet(x=ba_supervised[1:n,],y,xt=ba_supervised[-(1:n),],y_test)
#     
#     pheatmap(1-S,cluster_cols = F,cluster_rows = F)
#     pheatmap(1-S)
#     h1 <- spectral.clustering(S, n_eig = 2)
#     group1 <- which(h1==1)
#     group2 <- which(h1==-1)
#     h2 <- spectral.clustering(S[group1,group1], n_eig = 2)
#     group21 <- which(h2==1)
#     group22 <- which(h2==-1)
#     spectral.clustering(S[group21,group21], n_eig = 2)
#     spectral.clustering(S[group22,group22], n_eig = 2)
#     #   },
# #   "classo" = {
# #     z <- log(WC[1:n,])
# #     z_test <- log(WC[-(1:n),])
# #     fit_classo <- cv.func('ConstrLasso',y=y,x=z,C=matrix(1,p,1), nlam = 100, nfolds=10)
# #   },
# #   "propr" = {
# pr <- propr::propr(WC[1:n,], metric = "phs")
# h_propr <- hclust(as.dist(pr@matrix),method = linkage)
# plot(h_propr)
# sbp_propr <- sbp.fromHclust(h_propr)
# ba_propr <- balance.fromSBP(WC,sbp_propr)
# fit_propr <- run.glmnet(x=ba_propr[1:n,],y,xt=ba_propr[-(1:n),],y_test)
# #   },
# #   'pba' = {
# #     sbp_pba <- pba(WC)
# #     fit_pba <- run.glmnet(x=sbp_pba@pba[1:n,],y,xt=sbp_pba@pba[-(1:n),],y_test)
# #     # err[4,1] <- fit_pba$mse.pred
# #     # roc_pba <- apply(fit_pba$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_pba@sbp))
# #     #
# #   },times=10,unit = 's',
# #   check = NULL)
# # print(autoplot(mbm))
# err[1,1] <- fit_supervised$mse.pred
# roc_supervised <- apply(fit_supervised$beta, 2, function(a) roc.for.coef.LR(a, beta_oracle, sbp_supervised))
# err[2,1] <- fit_propr$mse.pred
# roc_propr <- apply(fit_propr$beta, 2, function(a) roc.for.coef.LR(a, beta_oracle, sbp_propr))
# # fit_oracle <- run.glmnet(x=ba_oracle[1:n,],y,xt=ba_oracle[-(1:n),],y_test)
# # err[5,1] <- fit_oracle$mse.pred
# # pred_classo_1 <- z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
# # err[5,1] <- mean((pred_classo_1 - y_test)^2)
# # roc_classo <- apply(fit_classo$bet, 2, function(a) roc.for.coef(a,beta_lc))
# # fit_selbal <- pred.from.selbal(WC[1:n,],y,WC[-(1:n),],y_test)
# # err[6,1] <- fit_selbal$mse.pred
# plot(roc_supervised[1,],roc_supervised[2,],type='l')
# lines(roc_propr[1,], roc_propr[2,],col='red')
# print(err)
# 
