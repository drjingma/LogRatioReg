rm(list=ls())
path = "Results/s2/"
library(ggplot2)
library(tidyverse)
filenames <- list.files(path)

res <- lapply(paste0(path,filenames),function(a) {load(a);return(res)})
err <- sapply(res,function(a) a$err[c(1,2,4,5),])
df15 <- data.frame(t(err))
colnames(df15) <- rownames(res[[1]]$err)[c(1,2,4,5)]
# print(summary(df1))

roc_sup = Reduce('+', lapply(res,function(a) a$sup) )/length(res)
# roc_coat = Reduce('+', lapply(res,function(a) a$coat) )/length(res)
roc_pba = Reduce('+', lapply(res,function(a) a$pba) )/length(res)
roc_propr = Reduce('+', lapply(res,function(a) a$propr) )/length(res)
roc_classo = Reduce('+', lapply(res,function(a) a$classo) )/length(res)

# df1 <- data.frame(rbind(df15,df12),'rho'=c(rep('5',nrow(df15)),rep('2',nrow(df12))))
df1 <- reshape2::melt(df15) %>%
  dplyr::rename(method=variable,MSE=value)%>%
  mutate('method'=factor(method,c('supervised','classo','coat','propr','pba')))
#
# my.labs <- list(bquote(rho==.(0.2)),bquote(rho==.(0.5)))

plt_err <- ggplot(df1, aes(x=method, y=MSE, fill=method)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4)  +
  scale_y_continuous(trans='log2') +
  theme(legend.position = 'none',axis.title.x = element_blank())
  # scale_x_discrete(labels=my.labs)

df1 %>%
  group_by(method) %>%
  summarise(mean=mean(MSE))

df <- rbind(
  mutate(as.data.frame(t(roc_classo)),'method'='classo'),
  mutate(as.data.frame(t(roc_propr)), 'method'='propr'),
  mutate(as.data.frame(t(roc_pba)),   'method'='pba'),
  mutate(as.data.frame(t(roc_sup)),   'method'='supervised')
) %>%
  mutate('method'=factor(method,c('supervised','classo','propr','pba')))

plt_roc <-  ggplot(df,aes(S_hat,tpr,color=method)) +
  geom_line(aes(linetype=method)) + geom_point(aes(pch=method)) +
  scale_linetype_manual(values=c('solid',"twodash","dotted",'dashed','dotdash')) +
  scale_shape_manual(values=c(16,17,18,8,9))+
  # scale_colour_manual(values=my_cols[c(3,1,2,4,5)])+
  xlim(0,40) +
  # ylim(0.4,1) +
  # coord_fixed() +
  # facet_grid(.~network) +
  xlab("Total number of discovery") + ylab("True positive rate") +
  theme(legend.position = c(0.65, 0.3))

# plt <- plot_grid(plt_err,plt_roc)

# + facet_wrap(pair~network,labeller = label_wrap_gen(multi_line=FALSE))
# save_plot(plt,filename = '../Output/aim3_expdecay2.png', base_asp = 2.5, base_height = 3.5)





# apply supervised log-ratios, using CV to select lambda
slr = cvSLR(y = y, X = WC[1:n,], intercept = T, nfolds = 10, lambda = lambda,
            rho.type = 'squared', linkage = linkage, standardize = T)
btree = slr$btree
# plot(btree)

# choose lambda
lam.min.idx = which.min(slr$cvm)
lam.min = slr$lambda[lam.min.idx]
a0 = slr$int[lam.min.idx]
thetahat = slr$bet[, lam.min.idx]
thetahat

Uhat = getU(btree = btree)
betahat = getBeta(thetahat, U = Uhat)

# evaluate model #

# 1. prediction error #
# 1a. on training set #
# get prediction error on training set
Yhat.train = a0 + computeBalances(WC[1:n,], btree) %*% thetahat
PE.train = as.vector(crossprod(y - Yhat.train) / n)
# 1b. on test set #
# get prediction error on test set
Yhat.test = a0 + computeBalances(WC[-c(1:n),], btree) %*% thetahat

PE.test = as.vector(crossprod(y_test - Yhat.test) / n)
# 2. estimation accuracy #
# 2a. estimation of beta #
EA1 = sum(abs(betahat - beta_lc))
EA2 = as.vector(sqrt(crossprod(betahat - beta_lc)))
EAInfty = max(abs(betahat - beta_lc))

SBP = sbp.fromHclust(btree)
non0.thetahat = (thetahat != 0)
sel.cols.SBP = SBP[, non0.thetahat]
non0.betahat = apply(sel.cols.SBP, 1, function(row) any(row != 0))
is0.betahat = !non0.betahat
# beta
non0.beta = abs(beta_lc) > 10e-8
is0.beta = abs(beta_lc) <= 10e-8
# FP
FP = sum(is0.beta & non0.betahat)
# FN
FN = sum((non0.beta != non0.betahat) & non0.beta)
# TPR
TPR = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
# beta sparsity
bspars = sum(non0.beta)