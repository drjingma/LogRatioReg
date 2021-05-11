rm(list=ls())

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

# set.seed(12)
n <- 100#50
p <- 200#20
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
seed = sample(1:1000, 1)
print(paste0("### ---------- ", seed, " ----------------------------------###"))
set.seed(seed)
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

################################################################################
# the original sup-balances code, for the most part, at least

# ## ---- Set the error matrix ----
# err <- matrix(NA,nrow=6,ncol=1)
# rownames(err) <- c('supervised','propr','coat','pba','classo','selbal')
# colnames(err) <- c('pred')
# 
# 
# ## Benchmark the computational time of various methods.
# mbm <- microbenchmark(
#   "supervised" = { 
#     S <- matrix(0,p,p)
#     rownames(S) <- colnames(S) <- colnames(W)
#     for (j in 1:(p-1)){
#       for (k in (j+1):p){
#         newx <- log(WC[1:n,j]) - log(WC[1:n,k])
#         newx <- newx - mean(newx)
#         newy <- y - mean(y)
#         S[j,k] <- S[k,j] <- abs(cor(newx,y))
#         # S[j,k] <- S[k,j] <- abs(crossprod(newx,newy)/(sqrt(crossprod(newx)) * sqrt(crossprod(newy))))
#       }
#     }
#     h_supervised <- hclust(as.dist(1-S),method = linkage)
#     # plot(h_supervised)
#     # ggtree(h_supervised) + geom_point(aes(shape=isTip, color=isTip), size=3)
#     sbp_supervised <- sbp.fromHclust(h_supervised)
#     ba_supervised <- balance.fromSBP(WC,sbp_supervised)
#     fit_supervised <- run.glmnet(x=ba_supervised[1:n,],y,xt=ba_supervised[-(1:n),],y_test)
#   },
#   "classo" = {
#     z <- log(WC[1:n,])
#     z_test <- log(WC[-(1:n),])
#     fit_classo <- cv.func('ConstrLasso',y=y,x=z,C=matrix(1,p,1), nlam = 100, nfolds=10)
#   },
#   'coat' = {
#     d_coat <- 1-coat(WC[1:n,])$corr
#     rownames(d_coat) <- colnames(d_coat) <- colnames(W)
#     sbp_coat <- sbp.fromHclust(hclust(as.dist(d_coat),method = linkage))
#     ba_coat <- balance.fromSBP(WC,sbp_coat)
#     fit_coat <- run.glmnet(x=ba_coat[1:n,],y,xt=ba_coat[-(1:n),],y_test)
#   },
#   "propr" = {
#     pr <- propr::propr(WC[1:n,], metric = "phs")
#     sbp_propr <- sbp.fromHclust(hclust(as.dist(pr@matrix),method = linkage))
#     ba_propr <- balance.fromSBP(WC,sbp_propr)
#     fit_propr <- run.glmnet(x=ba_propr[1:n,],y,xt=ba_propr[-(1:n),],y_test)
#   },
#   'pba' = {
#     sbp_pba <- pba(WC)
#     fit_pba <- run.glmnet(x=sbp_pba@pba[1:n,],y,xt=sbp_pba@pba[-(1:n),],y_test)
#     # err[4,1] <- fit_pba$mse.pred
#     # roc_pba <- apply(fit_pba$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_pba@sbp))
#     # 
#   },times = 1, #times=10,
#   unit = 's',
#   check = NULL)
# # print(autoplot(mbm))
# 
# err[1,1] <- fit_supervised$mse.pred
# # roc_supervised <- apply(fit_supervised$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_supervised))
# roc_supervised <- apply(fit_supervised$beta, 2, 
#                         function(a) tpr.for.coef.ilr(beta_lc, a, sbp_supervised))
# 
# 
# err[2,1] <- fit_propr$mse.pred
# # roc_propr <- apply(fit_propr$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_propr))
# roc_propr <- apply(fit_propr$beta, 2, function(a) tpr.for.coef.ilr(beta_lc, a, sbp_propr))
# 
# err[3,1] <- fit_coat$mse.pred
# # roc_coat <- apply(fit_coat$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_coat))
# roc_coat <- apply(fit_coat$beta, 2, function(a) tpr.for.coef.ilr(beta_lc, a, sbp_coat))
# 
# pred_classo_1 <- z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
# err[5,1] <- mean((pred_classo_1 - y_test)^2)
# # roc_classo <- apply(fit_classo$bet, 2, function(a) roc.for.coef(a,beta_lc))
# roc_classo <- apply(fit_classo$bet, 2, function(a) tpr.for.coef(beta_lc, a))
# 
# pred_classo_2 <- fit_classo$int[which.min(fit_classo$cvm)] + 
#   z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
# err[5,1] <- mean((pred_classo_2 - y_test)^2)
# # 
# # # fit_selbal <- pred.from.selbal(WC[1:n,],y,WC[-(1:n),],y_test)
# # # err[6,1] <- fit_selbal$mse.pred
# # 
# # print(err)
# 
# err[4,1] <- fit_pba$mse.pred
# # roc_pba <- apply(fit_pba$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_pba@sbp))
# roc_pba <- apply(fit_pba$beta, 2, function(a) tpr.for.coef.ilr(beta_lc, a, sbp_pba@sbp))
# 
# roc_dat = data.frame(
#   S.hat = c(roc_classo[1, ], 
#             roc_supervised[1, ], 
#             roc_propr[1, ],
#             roc_coat[1, ], 
#             roc_pba[1, ]
#             ), 
#   TPR = c(roc_classo[2, ], 
#           roc_supervised[2, ], 
#           roc_propr[2, ],
#           roc_coat[2, ],
#           roc_pba[2, ]
#           ), 
#   type = c(rep("classo", ncol(roc_classo)), 
#            rep("slr", ncol(roc_supervised)),
#            rep("propr", ncol(roc_propr)),
#            rep("coat", ncol(roc_coat)), 
#            rep("pba", ncol(roc_pba))
#            )
# )
# roc_dat$type = factor(roc_dat$type, 
#                       levels = c("slr", "classo", "coat", "propr", "pba"))
# ggplot(roc_dat, aes(x = S.hat, y = TPR, color = type, linetype = type, shape = type)) + 
#   geom_path() + 
#   geom_point(size = 3) +
#   # xlim(0, 40) +
#   theme(text = element_text(size = 20))
# print(paste0("### ---------- ", seed, " ----------------------------------###"))
# err
















################################################################################
# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

numSims = 50

evals = foreach(
  b = 1:numSims, 
  .combine = cbind, 
  .noexport = c("ConstrLassoC0")
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  source("RCode/func_libs.R")
  
  # data
  ## Generate latent variables from a log-normal distribution
  ## Observed compositions are defined from the latent variables after dividing by total.
  logW <- mvrnorm(n=n*2, mu=mu, Sigma=Sigma)
  W <- exp(logW) # basis
  colnames(W) <- paste0('s',1:p)
  WC <- sweep(W,1,rowSums(W), FUN='/')
  WCLR <- t(apply(WC,1,clr))
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
    "classo" = {
      z <- log(WC[1:n,])
      z_test <- log(WC[-(1:n),])
      fit_classo <- cv.func('ConstrLasso',y=y,x=z,C=matrix(1,p,1), nlam = 100, nfolds=10)
    },
    'coat' = {
      d_coat <- 1-coat(WC[1:n,])$corr
      rownames(d_coat) <- colnames(d_coat) <- colnames(W)
      sbp_coat <- sbp.fromHclust(hclust(as.dist(d_coat),method = linkage))
      ba_coat <- balance.fromSBP(WC,sbp_coat)
      fit_coat <- run.glmnet(x=ba_coat[1:n,],y,xt=ba_coat[-(1:n),],y_test)
    },
    "propr" = {
      pr <- propr::propr(WC[1:n,], metric = "phs")
      sbp_propr <- sbp.fromHclust(hclust(as.dist(pr@matrix),method = linkage))
      ba_propr <- balance.fromSBP(WC,sbp_propr)
      fit_propr <- run.glmnet(x=ba_propr[1:n,],y,xt=ba_propr[-(1:n),],y_test)
    },
    'pba' = {
      sbp_pba <- pba(WC)
      fit_pba <- run.glmnet(x=sbp_pba@pba[1:n,],y,xt=sbp_pba@pba[-(1:n),],y_test)
      # err[4,1] <- fit_pba$mse.pred
      # roc_pba <- apply(fit_pba$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_pba@sbp))
      # 
    },times = 1, #times=10,
    unit = 's',
    check = NULL)
  # print(autoplot(mbm))
  
  err[1,1] <- fit_supervised$mse.pred
  # roc_supervised <- apply(fit_supervised$beta, 2, 
  #                         function(a) tpr.for.coef.ilr(beta_lc, a, sbp_supervised))
  
  
  err[2,1] <- fit_propr$mse.pred
  # roc_propr <- apply(fit_propr$beta, 2, function(a) tpr.for.coef.ilr(beta_lc, a, sbp_propr))
  
  err[3,1] <- fit_coat$mse.pred
  # roc_coat <- apply(fit_coat$beta, 2, function(a) tpr.for.coef.ilr(beta_lc, a, sbp_coat))
  
  pred_classo_1 <- z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
  err[5,1] <- mean((pred_classo_1 - y_test)^2)
  # roc_classo <- apply(fit_classo$bet, 2, function(a) tpr.for.coef(beta_lc, a))
  
  # pred_classo_2 <- fit_classo$int[which.min(fit_classo$cvm)] +
  #   z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
  # err[5,1] <- mean((pred_classo_2 - y_test)^2) 
  
  err[4,1] <- fit_pba$mse.pred
  # roc_pba <- apply(fit_pba$beta, 2, function(a) tpr.for.coef.ilr(beta_lc, a, sbp_pba@sbp))
  
  err
}


dim(evals)
eval.means = apply(evals, 1, mean)
evals.df = data.frame(t(evals))
data.gg = melt(evals.df)
data.gg$variable = factor(data.gg$variable, 
                          levels = c("supervised", "classo", "coat", "propr", "pba"))
plt = ggplot(data.gg, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot() + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank())

ggsave("MSE.pdf", 
       plot = plt, 
       width = 6, 
       height = 4, 
       units = "in")
