# rm(list=ls())

library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# library(phyloseq)
library(philr)
library(selbal)
library(propr)
library(balance)
library(ape)
library(DiagrammeR)
library(igraph)
library(Matrix)
library(microbenchmark)

# CodeFolder <- "../RCode/"
# OutFolder <- "../Output/"
source("RCode/func_libs.R")
source("COAT-master/coat.R")

set.seed(12)
n <- 50
p <- 20
# G <- diag(p) - tcrossprod(rep(1,p))/p
# U <- matrix(c(0.5,0.5,-0.5,-0.5,
#               sqrt(1/2),-sqrt(1/2),0,0,
#               0,0,sqrt(1/2),-sqrt(1/2)),ncol=3)
# Uinv <- corpcor::pseudoinverse(t(U))

sig <- 0.5
rho <- 0.5
linkage <- 'average'

Sigma <- rgExpDecay(p,rho)$Sigma

## Generate latent variables from a log-normal distribution
## Observed compositions are defined from the latent variables after dividing by total.
mu <- c(rep(log(p),5),rep(0,p-5))
logW <- mvrnorm(n=n*2, mu=mu, Sigma=Sigma)
W <- exp(logW) # basis
colnames(W) <- paste0('s',1:p)
WC <- sweep(W,1,rowSums(W), FUN='/')
WCLR <- t(apply(WC,1,clr))

## Generate coefficient
# beta_lc <- c(1,0,-0.8,0,0.4,0,-1.5,0,1.2,0,-0.3,rep(0,p-11))
beta_lc <- c(1,0.4,1.2,-1.5,-0.8,-0.3,rep(0,p-6))
# beta_lc <- c(1,-0.8,0.4,0,-0.6,0,-1.5,0,1.2,0,0.3,rep(0,p-11))
beta_support <- which(abs(beta_lc)>0)
names(beta_lc) <- colnames(W)

## Generate response
yAll <-  log(WC) %*% beta_lc + rnorm(n) * sig
y <- yAll[1:n,]
y_test <- yAll[-(1:n),]

## ---- Set the error matrix ----
err <- matrix(NA,nrow=6,ncol=1)
rownames(err) <- c('supervised','propr','coat','pba','classo','selbal')
colnames(err) <- c('pred')


## Benchmark the computational time of various methods.
mbm <- microbenchmark(
  "supervisedKristyn" = {
    btree = getSupervisedTree(WC[1:n, ], y, linkage)
    Xb = computeBalances(btree, WC)
    fit_supervisedKristyn <- run.glmnet(x=Xb[1:n,],y,xt=Xb[-(1:n),],y_test)
  },
  "supervised" = { 
    S <- matrix(0,p,p)
    rownames(S) <- colnames(S) <- colnames(W)
    for (j in 1:(p-1)){
      for (k in (j+1):p){
        newx <- log(WC[1:n,j]) - log(WC[1:n,k])
        newx <- newx - mean(newx)
        newy <- y - mean(y)
        S[j,k] <- S[k,j] <- abs(cor(newx,y))
        # S[j,k] <- S[k,j] <- abs(crossprod(newx,newy)/(sqrt(crossprod(newx)) * sqrt(crossprod(newy))))
      }
    }
    h_supervised <- hclust(as.dist(1-S),method = linkage)
    # plot(h_supervised)
    # ggtree(h_supervised) + geom_point(aes(shape=isTip, color=isTip), size=3)
    sbp_supervised <- sbp.fromHclust(h_supervised)
    ba_supervised <- balance.fromSBP(WC,sbp_supervised)
    fit_supervised <- run.glmnet(x=ba_supervised[1:n,],y,xt=ba_supervised[-(1:n),],y_test)
  },
  # "classo" = {
    # z <- log(WC[1:n,])
    # z_test <- log(WC[-(1:n),])
  #   fit_classo <- cv.func('ConstrLasso',y=y,x=z,C=matrix(1,p,1), nlam = 100, nfolds=10)
  # },
  'coat' = {
    d_coat <- 1-coat(WC[1:n,])$corr
    rownames(d_coat) <- colnames(d_coat) <- colnames(W)
    sbp_coat <- sbp.fromHclust(hclust(as.dist(d_coat),method = linkage))
    ba_coat <- balance.fromSBP(WC,sbp_coat)
    fit_coat <- run.glmnet(x=ba_coat[1:n,],y,xt=ba_coat[-(1:n),],y_test)
  },
  # "propr" = {
  #   pr <- propr::propr(WC[1:n,], metric = "phs")
  #   sbp_propr <- sbp.fromHclust(hclust(as.dist(pr@matrix),method = linkage))
  #   ba_propr <- balance.fromSBP(WC,sbp_propr)
  #   fit_propr <- run.glmnet(x=ba_propr[1:n,],y,xt=ba_propr[-(1:n),],y_test)
  # },  
  'pba' = {
    sbp_pba <- pba(WC)
    fit_pba <- run.glmnet(x=sbp_pba@pba[1:n,],y,xt=sbp_pba@pba[-(1:n),],y_test)
    # err[4,1] <- fit_pba$mse.pred
    # roc_pba <- apply(fit_pba$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_pba@sbp))
    # 
  },times=10,unit = 's',
  check = NULL)
# print(autoplot(mbm))

# err[1,1] <- fit_supervised$mse.pred
# roc_supervised <- apply(fit_supervised$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_supervised))

# err[2,1] <- fit_propr$mse.pred
# roc_propr <- apply(fit_propr$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_propr))

# err[3,1] <- fit_coat$mse.pred
# roc_coat <- apply(fit_coat$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_coat))

# pred_classo_1 <- z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
# err[5,1] <- mean((pred_classo_1 - y_test)^2)
# roc_classo <- apply(fit_classo$bet, 2, function(a) roc.for.coef(a,beta_lc))
#
# fit_selbal <- pred.from.selbal(WC[1:n,],y,WC[-(1:n),],y_test)
# # err[6,1] <- fit_selbal$mse.pred
#

compare.mse = c(fit_supervised$mse.pred, fit_supervisedKristyn$mse.pred, fit_pba$mse.pred)
names(compare.mse) = c("slr", "slr-kristyn", "pba")
print(compare.mse)



