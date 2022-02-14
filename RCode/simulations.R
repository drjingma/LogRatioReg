## Generate parameters
rm(list=ls())
jid <- 12
set.seed(jid)
# devtools::install_github(repo = "malucalle/selbal")
# Aim 3: balances regression
library(philr)
library(propr)
library(balance)
library(tidyverse)
library(cowplot)
library(ape)
library(DiagrammeR)
# library(selbal)
library(igraph)
library(Matrix)
library(fields)
library(pheatmap)

CodeFolder <- "RCode/"
source(paste0(CodeFolder,"func_libs.R"))
source("RCode/slr-main.R")

DEBUG <- 'aim3'
today <- '20220201'
linkage <- 'average'   # controls the hclust
rho <- 0.5 # magnitude of the correlation in cov Sigma
p <- 60   # number of variables 
n <- 100   # number of samples
sig <- 0.5 # noise level

today <- paste0(DEBUG,'_',today,"_p",p,"_n",n,'_sig',sig*10,'_rho',rho*10)

Sigma <- rgExpDecay(p,rho)$Sigma
# sds <- rlnorm(p,meanlog = 0.16, sdlog = 0.45)
# Sigma <- diag(sds) %*% Sigma %*% diag(sds)
rownames(Sigma) <- colnames(Sigma) <- paste0('s',1:p)
h_oracle <- hclust(as.dist(1-cov2cor(Sigma)),method = linkage)
sbp_oracle <- sbp.fromHclust(h_oracle)

## Generate compositions from normalized log-normal
mu <- c(rep(log(p),5),rep(0,p-5))
logW <- mvrnorm(n=n*2, mu=mu, Sigma=diag(p)) 
W <- exp(logW) # basis
colnames(W) <- paste0('s',1:p)
WC <- sweep(W,1,rowSums(W), FUN='/')
WCLR <- t(apply(WC,1,clr))

## Generate response ----
theta <- 0.5
yAll <-  (log(WC[,1]) + log(WC[,2]) + log(WC[,3]) - log(WC[,4]) - log(WC[,5]) - log(WC[,6])) * theta + rnorm(n) * sig
y_train <- yAll[1:n]
y_test <- yAll[-(1:n)]
X_train <- WC[1:n,]
X_test <- WC[-c(1:n),]

out <- slr(WC,yAll)
print(out$index)

# asinTransform <- function(p) { asin(sqrt(p)) }
# plot(svd(rhoMat)$d)
# pheatmap(rhoMat,show_colnames = T, show_rownames = T, 
#          cluster_rows = F, cluster_cols = F,
#          main = 'Heatmap of the original correlation')
# pheatmap(asinTransform(rhoMat),show_colnames = T, show_rownames = T, 
#          cluster_rows = F, cluster_cols = F,
#          main = 'Heatmap of the arcsin correlation')
# spectral.clustering(asinTransform(rhoMat))
# spectral.clustering(rhoMat[1:6,1:6])

# head(sort(apply(rhoMat,1,sum), decreasing = T))
# setdiff(names(head(sort(apply(rhoMat,1,max), decreasing = T))),names(index))
# setdiff(names(head(sort(apply(rhoMat,1,sum), decreasing = T))),names(index))
# a <- apply(rhoMat,1,max)
# b <- apply(rhoMat,1,sum)
# a <- a/max(a)
# b <- b/max(b)
# alpha <- 1
# head(sort(alpha*a + (1-alpha)*b,decreasing = T))

# rhoMat.svd <- svd(rhoMat)
# rhoMat_approx_1 <-  tcrossprod(rhoMat.svd$u[,1], rhoMat.svd$v[,1]) * rhoMat.svd$d[1]
# rownames(rhoMat_approx_1) <- colnames(rhoMat_approx_1) <- rownames(rhoMat)
# pheatmap(rhoMat_approx_1,show_colnames = T, show_rownames = T, cluster_rows = F, cluster_cols = F,
#          main = 'Heatmap of first SVD approximated correlation')
# cl <- spectral.clustering(rhoMat_approx_1)
# index <- which(cl==1)
# spectral.clustering(rhoMat[index,index])

# setdiff(names(head(sort(apply(rhoMat_approx_1,1,max), decreasing = T),n=sum(abs(sbp_oracle[,jid])))),names(index))
# setdiff(names(head(sort(apply(rhoMat_approx_1,1,sum), decreasing = T),n=sum(abs(sbp_oracle[,jid])))),names(index))
# 
# rhoMat_approx_2 <-  rhoMat.svd$u[,1:2] %*% diag(rhoMat.svd$d[1:2]) %*% t(rhoMat.svd$v[,1:2])
# rownames(rhoMat_approx_2) <- colnames(rhoMat_approx_2) <- rownames(rhoMat)
# pheatmap(rhoMat_approx_2,show_colnames = T, show_rownames = T, cluster_rows = F, cluster_cols = F)
# spectral.clustering(rhoMat_approx_2)
# head(sort(apply(rhoMat_approx_2,1,max), decreasing = T))
# head(sort(apply(rhoMat_approx_2,1,sum), decreasing = T))

# rho4 <- rhoMat_approx_1
# rho4[which(rho4<0.55)] <- 0
# pheatmap(rho4,show_colnames = T, show_rownames = T, cluster_rows = F, cluster_cols = F)
# rho4 <- rhoMat
# rho4[which(rho4<0.27)] <- 0
# plot(hclust(as.dist(1-rhoMat)))

# spectral.clustering(rhoMat)

# library(pROC)
# library(ggplot2)
# source("RCode/Selbal_Functions.R")
# out <- selbal(WC,yAll)

## ---- Set the error matrix ----
err <- matrix(NA,nrow=6,ncol=2)
rownames(err) <- c('slr','propr','coat','pba','classo','selbal')
colnames(err) <- c('training', 'test')
x = X_train; y=y_train; lambda=seq(0.1,0.7,0.02)
# fit <- cv.slr(x = X_train, y=y_train, lambda=seq(0.3,0.7,0.02))
# refit <- slr(X_train, y_train, lambda=0.75)
# err[1,1] <- mean((y_train - refit$fitted.values)^2)
# err[1,2] <- mean((y_test - predict(refit,X_test))^2)
## Tuning parameters

# ## Construct SBP from proportionality matrix
# pr <- propr::propr(WC[1:n,], metric = "phs")
# sbp_propr <- sbp.fromHclust(hclust(as.dist(pr@matrix),method = linkage))
# ba_propr <- balance.fromSBP(WC,sbp_propr)
# fit_propr <- run.glmnet(x=ba_propr[1:n,],y,xt=ba_propr[-(1:n),],y_test,lambda = lambda)
# err[2,1] <- fit_propr$mse.pred
# roc_propr <- apply(fit_propr$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_propr))
# 
# # d_coat <- 1-coat(WC[1:n,])$corr
# # rownames(d_coat) <- colnames(d_coat) <- colnames(W)
# # sbp_coat <- sbp.fromHclust(hclust(as.dist(d_coat),method = linkage))
# # ba_coat <- balance.fromSBP(WC,sbp_coat)
# # fit_coat <- run.glmnet(x=ba_coat[1:n,],y,xt=ba_coat[-(1:n),],y_test)
# # err[3,1] <- fit_coat$mse.pred
# # roc_coat <- apply(fit_coat$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_coat))
# 
# ## Construct principal balances
# sbp_pba <- pba(WC)
# fit_pba <- run.glmnet(x=sbp_pba@pba[1:n,],y,xt=sbp_pba@pba[-(1:n),],y_test,lambda = lambda)
# err[4,1] <- fit_pba$mse.pred
# roc_pba <- apply(fit_pba$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_pba@sbp))
# 
## Fit compositional Lasso using the same lambda
# nlam <- 100
# maxlam <- 2*max(abs(crossprod(log(X_train),y_train)/n))
# lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))
# fit <- ConstrLasso(y=y_train,x=log(X_train),Cmat=matrix(1,p,1))
# fit_classo <- cv.func('ConstrLasso',y=y_train,x=log(X_train),Cmat=matrix(1,p,1), lambda = lambda, nfolds=10)
# pred_classo <- fit_classo$int[which.min(fit_classo$cvm)] + log(X_test) %*% fit_classo$bet[,which.min(fit_classo$cvm)]
# err[5,2] <- mean((pred_classo - y_test)^2)
# err[5,1] <- mean((fit_classo$int[which.min(fit_classo$cvm)] + log(X_train) %*% fit_classo$bet[,which.min(fit_classo$cvm)] - y_train)^2)
# roc_classo <- apply(fit_classo$bet, 2, function(a) roc.for.coef(a,beta_lc))
# 
# print(err)
# plot(roc_supervised[1,],roc_supervised[2,],type='b')
# lines(roc_classo[1,],roc_classo[2,],col='red')
# 
# res <- list(sup=roc_supervised,
#             pba=roc_pba,
#             propr=roc_propr,
#             classo=roc_classo,
#             err=err, beta=beta_lc, mu=mu)
# 
# if (!is.null(res)) {save(res, file=paste0(OutFolder, today, '_',jid,'.rda'))}
# 
