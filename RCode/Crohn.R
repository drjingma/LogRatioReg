rm(list=ls())
library(selbal)
CodeFolder <- "RCode/"
source(paste0(CodeFolder,"func_libs.R"))
source("RCode/SLR.R")

library(pheatmap)
library(dplyr)
library(ggplot2)
# library("RColorBrewer")
# display.brewer.all()

x3 <- Crohn[,1:48]
y3 <- Crohn[,49]

# zeros in the covariate 
qplot(apply(x3==0,1,mean))

# check number of samples where each pair of variables are nonzero


# Run selbal.cv function
# BAL.Crohn <- selbal.cv(x = x3, y = y3, n.fold = 5, n.iter = 10,
#                        covar = NULL, logit.acc = "AUC")
# subset_selbal <- ifelse(BAL.Crohn$global.balance$Group=='DEN',1,-1)
# names(subset_selbal) <- BAL.Crohn$global.balance$Taxa
# zS <- balance.fromSBP(x3.adj[,match(names(subset_selbal),colnames(x3.adj))],as.matrix(subset_selbal))
# cor(zS,y)

# Run SLR with 2-level spectral clustering
y <- 2*(as.numeric(y3)-1.5)
x3.adj <- cmultRepl2(x3, zero.rep = 'bayes') # relative abundances
x3.clr <- t(apply(x3.adj,1,alr))
cm <- cor(x3.clr)
p <- ncol(x3.adj)
qplot(abs(cor(x3.adj,y)))
object <- slr(x3.adj,y,method='correlation',response.type='binary',threshold = 0.19,s0.perc = 0)
predictor <- predict.slr(object,newdata=x3.adj,response.type = 'binary')
misClasificError <- mean(predictor != y)
print(paste('Accuracy',1-misClasificError))

# object <- cv.slr(x3.adj,y,method='wald',threshold = seq(0.5,0.99,0.01),s0.perc = 0)
object <- cv.slr(x3.adj,y,method='correlation',response.type='binary',threshold = seq(0.1,0.25,0.01),s0.perc = 0)
plot(object$cvm)
object$threshold.min


rhoMat <- matrix(0,p,p)
rownames(rhoMat) <- colnames(rhoMat) <- colnames(x3.adj)
for (j in 1:p){
  for (k in 1:p){
    if (k==j){next}
    else {
      # testmat[j,k] <- cor.test(log(x[,j])-log(x[,k]),y)$p.value
      rhoMat[j,k] <- abs(stats::cor(log(x3.adj[,j])-log(x3.adj[,k]),y))
    }
  }
}
W <- max(rhoMat)-rhoMat
pheatmap(W,show_colnames = T, show_rownames = T,
         cluster_rows = T, cluster_cols = T,
         main = 'Heatmap of the original correlation')

# qplot(rowSums(max(rhoMat)-rhoMat)) # need to regularize the high degree. 

L <- graph.laplacian(W,normalized = TRUE,zeta=0.01)
L.eig <- eigen(L)
range(L.eig$values)
plot(L.eig$values)
df <- data.frame('PC1'=L.eig$vectors[,1],
                 'PC2'=L.eig$vectors[,2],'PC3'=L.eig$vectors[,3]) %>%
  reshape2::melt() %>%
  ggplot(aes(x=rep(1:p,3),y=value)) + geom_point() +
  facet_wrap(~variable)
print(df) 

cl <- kmeans(L.eig$vectors[,1:3],centers=3,nstart = 100)$cluster
## Build balance from these three subsets
SBP <- matrix(0,nrow=length(cl),ncol=3)
rownames(SBP) <- colnames(x3.adj)
SBP[(cl==1),1] <- 1
SBP[(cl==2),1] <- -1
SBP[(cl==1),2] <- 1
SBP[(cl==3),2] <- -1
SBP[(cl==2),3] <- 1
SBP[(cl==3),3] <- -1
apply(abs(SBP),2,sum)
z <- balance::balance.fromSBP(x3.adj,SBP)
cor(z,y)
index <- which.max(abs(cor(z,y)))
plot(z[,index],y)

# df <- data.frame()
summary(glm(y/2~z[,index],family = 'binomial'))


subset <- spectral.clustering(max(rhoMat)-rhoMat,2)
subset2 <- spectral.clustering.sign(1-rhoMat)

zz <- balance.fromSBP(x3.adj[,match(names(subset),colnames(x3.adj))],as.matrix(subset))
cor(zz,y)

## view heatmap with cluster labels
metadata_col <- data.frame(cluster=as.factor(subset),
                           row.names=colnames(rhoMat))
metadata_row <- data.frame(cluster=as.factor(subset),
                           row.names=rownames(rhoMat))
# ann_colors = list(
#   cluster = c('1'='#4DAF4A',"2"='#984EA3',"3"='#999999')
# )
ann_colors = list(
  cluster = c('1'='#4DAF4A',"-1"='#984EA3')
)
W <- 1-rhoMat

plt.rho.joint <- pheatmap(1-rhoMat,show_rownames = T,show_colnames = T, 
                          # treeheight_row = 0, treeheight_col = 0,
                          # color = viridis(256),
                          # fontsize = 26,
                          annotation_colors=ann_colors,
                          annotation_row = metadata_row,annotation_col = metadata_col)

## Perform spectral clustering on each subset of variables using the original correlation matrix
index <- which(subset==1)
rhoMat_A <- rhoMat[index,index]
rhoMat_B <- rhoMat[-index,-index]
pheatmap(1-rhoMat_A,show_colnames = T, show_rownames = T,
         cluster_rows = T, cluster_cols = T,
         main = 'Heatmap of the original correlation')
pheatmap(1-rhoMat_B,show_colnames = T, show_rownames = T,
         cluster_rows = T, cluster_cols = T,
         main = 'Heatmap of the original correlation')
# intersect(BAL.Crohn$global.balance$Taxa, rownames(rhoMat_A))
# intersect(BAL.Crohn$global.balance$Taxa, rownames(rhoMat_B))

# 1st subset
# rhoMat_A <- rhoMat_A
# if (APPROXIMATION){
#   rhoMat_A <- low.rank(rhoMat_A,1)
# }
subset_A <- spectral.clustering.sign(1-rhoMat_A)
table(subset_A)
zA <- balance.fromSBP(x3.adj[,match(names(subset_A),colnames(x3.adj))],as.matrix(subset_A))
cor(zA,y)

# table(spectral.clustering(rhoMat[index,index]))
# Cluster the first subset again
rhoMat_AA <- rhoMat_A[which(subset_A==1),which(subset_A==1)]
rhoMat_AB <- rhoMat_A[-which(subset_A==1),-which(subset_A==1)]
# if (APPROXIMATION){
#   # rhoMat_AA <- low.rank(rhoMat_AA,1)
#   rhoMat_AB <- low.rank(rhoMat_AB,1)
# }
pheatmap(rhoMat_AB-diag(diag(rhoMat_AB)),show_colnames = T, show_rownames = T,
         cluster_rows = T, cluster_cols = T,
         main = 'Heatmap of the original correlation')

subset_AA <- spectral.clustering.sign(1-rhoMat_AA)
table(subset_AA)
subset_AB <- spectral.clustering.sign(1-rhoMat_AB)
table(subset_AB)
z_AA <- balance.fromSBP(x3.adj[,match(names(subset_AA),colnames(x3.adj))],as.matrix(subset_AA))
z_AB <- balance.fromSBP(x3.adj[,match(names(subset_AB),colnames(x3.adj))],as.matrix(subset_AB))
cor(z_AA,y)
cor(z_AB,y)

# 2nd subset
# rhoMat_B <- rhoMat_B
# if (APPROXIMATION){
#   rhoMat_B <- low.rank(rhoMat_B,1)
# } 
subset_B <- spectral.clustering.sign(1-rhoMat_B)
table(subset_B)
zB <- balance.fromSBP(x3.adj[,match(names(subset_B),colnames(x3.adj))],as.matrix(subset_B))
cor(zB,y)

# Cluster the second subset again
rhoMat_BA <- rhoMat_B[which(subset_B==1),which(subset_B==1)]
rhoMat_BB <- rhoMat_B[-which(subset_B==1),-which(subset_B==1)]
# if (APPROXIMATION){
#   rhoMat_BA <- low.rank(rhoMat_BA,1)
#   rhoMat_BB <- low.rank(rhoMat_BB,1)
# }
subset_BA <- spectral.clustering.sign(1-rhoMat_BA)
table(subset_BA)
subset_BB <- spectral.clustering.sign(1-rhoMat_BB)
table(subset_BB)
z_BA <- balance.fromSBP(x3.adj[,match(names(subset_BA),colnames(x3.adj))],as.matrix(subset_BA))
z_BB <- balance.fromSBP(x3.adj[,match(names(subset_BB),colnames(x3.adj))],as.matrix(subset_BB))
cor(z_BA,y)
cor(z_BB,y)

# pheatmap(rhoMat_AA,show_colnames = T, show_rownames = T,
#          cluster_rows = T, cluster_cols = T,
#          main = 'Heatmap of the original correlation')
# pheatmap(rhoMat_AB,show_colnames = T, show_rownames = T,
#          cluster_rows = T, cluster_cols = T,
#          main = 'Heatmap of the original correlation')
# intersect(BAL.Crohn$global.balance$Taxa, rownames(rhoMat_AB))
# intersect(BAL.Crohn$global.balance$Taxa, rownames(rhoMat_AA))

# out <- slr(x,y)
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
# rhoMat <-  tcrossprod(rhoMat.svd$u[,1], rhoMat.svd$v[,1]) * rhoMat.svd$d[1]
# rownames(rhoMat) <- colnames(rhoMat) <- rownames(rhoMat)
# pheatmap(rhoMat,show_colnames = T, show_rownames = T, cluster_rows = F, cluster_cols = F,
#          main = 'Heatmap of first SVD approximated correlation')
# cl <- spectral.clustering(rhoMat)
# index <- which(cl==1)
# spectral.clustering(rhoMat[index,index])

# setdiff(names(head(sort(apply(rhoMat,1,max), decreasing = T),n=sum(abs(sbp_oracle[,jid])))),names(index))
# setdiff(names(head(sort(apply(rhoMat,1,sum), decreasing = T),n=sum(abs(sbp_oracle[,jid])))),names(index))
# 
# rhoMat_2 <-  rhoMat.svd$u[,1:2] %*% diag(rhoMat.svd$d[1:2]) %*% t(rhoMat.svd$v[,1:2])
# rownames(rhoMat_2) <- colnames(rhoMat_2) <- rownames(rhoMat)
# pheatmap(rhoMat_2,show_colnames = T, show_rownames = T, cluster_rows = F, cluster_cols = F)
# spectral.clustering(rhoMat_2)
# head(sort(apply(rhoMat_2,1,max), decreasing = T))
# head(sort(apply(rhoMat_2,1,sum), decreasing = T))

# rho4 <- rhoMat
# rho4[which(rho4<0.55)] <- 0
# pheatmap(rho4,show_colnames = T, show_rownames = T, cluster_rows = F, cluster_cols = F)
# rho4 <- rhoMat
# rho4[which(rho4<0.27)] <- 0
# plot(hclust(as.dist(1-rhoMat)))

# spectral.clustering(rhoMat)

# library(pROC)
# library(ggplot2)
# source("RCode/Selbal_Functions.R")
# out <- selbal(x,yAll)

## ---- Set the error matrix ----
# err <- matrix(NA,nrow=6,ncol=2)
# rownames(err) <- c('slr','propr','coat','pba','classo','selbal')
# colnames(err) <- c('training', 'test')
# x = X_train; y=y_train; lambda=seq(0.1,0.7,0.02)
# fit <- cv.slr(x = X_train, y=y_train, lambda=seq(0.3,0.7,0.02))
# refit <- slr(X_train, y_train, lambda=0.75)
# err[1,1] <- mean((y_train - refit$fitted.values)^2)
# err[1,2] <- mean((y_test - predict(refit,X_test))^2)
## Tuning parameters

# ## Construct SBP from proportionality matrix
# pr <- propr::propr(x[1:n,], metric = "phs")
# sbp_propr <- sbp.fromHclust(hclust(as.dist(pr@matrix),method = linkage))
# ba_propr <- balance.fromSBP(x,sbp_propr)
# fit_propr <- run.glmnet(x=ba_propr[1:n,],y,xt=ba_propr[-(1:n),],y_test,lambda = lambda)
# err[2,1] <- fit_propr$mse.pred
# roc_propr <- apply(fit_propr$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_propr))
# 
# # d_coat <- 1-coat(x[1:n,])$corr
# # rownames(d_coat) <- colnames(d_coat) <- colnames(W)
# # sbp_coat <- sbp.fromHclust(hclust(as.dist(d_coat),method = linkage))
# # ba_coat <- balance.fromSBP(x,sbp_coat)
# # fit_coat <- run.glmnet(x=ba_coat[1:n,],y,xt=ba_coat[-(1:n),],y_test)
# # err[3,1] <- fit_coat$mse.pred
# # roc_coat <- apply(fit_coat$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_coat))
# 
# ## Construct principal balances
# sbp_pba <- pba(x)
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
