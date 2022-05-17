set.seed(11)
library(philr)
library(propr)
library(balance)
library(tidyverse)
library(cowplot)
library(ape)
library(DiagrammeR)
library(igraph)
library(Matrix)

library(phyloseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
theme_set(theme_cowplot())
load("../../Microbiome/Data/TwinsUK/twinsUK.rda")
source("../../Library/Grand_R_Functions.R")
# source("libs_hurdle.R")

twins2Physeq1
summary(colSums(otu_table(twins2Physeq1)))
qplot(log10(colSums(otu_table(twins2Physeq1))),binwidth=0.1) +
  xlab("Logged counts-per-sample")

# Remove one sample that has low depth
twins2Physeq1 <- prune_samples(colSums(otu_table(twins2Physeq1)) > 1000, twins2Physeq1)
summary(colSums(otu_table(twins2Physeq1)))
depths <- colSums(otu_table(twins2Physeq1))
qplot(log10(depths),binwidth=0.1) + xlab("Logged counts-per-sample")

# Remove family dependency
twins2Physeq2 <- prune_samples(!duplicated(sample_data(twins2Physeq1)$family_id), twins2Physeq1)
depths <- colSums(otu_table(twins2Physeq2))
qplot(log10(depths),binwidth=0.1) + xlab("Logged counts-per-sample")

# Add taxonomy and total read counts to this data.frame
raw.counts <- t(as(otu_table(twins2Physeq2), "matrix")) # n by p
prevdf = data.frame(prevalence = colMeans(raw.counts>0),
                    TotalAbundance = taxa_sums(twins2Physeq2),
                    as(tax_table(twins2Physeq2), "matrix"))

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf)[(prevdf$prevalence >= 0.1)]
twins2Physeq3 = prune_taxa(keepTaxa, twins2Physeq2)

X <- 1+t(as(otu_table(twins2Physeq3), "matrix")) # n by p
p <- ncol(X)
n <- nrow(X)
y <- log(as(sample_data(twins2Physeq3), "data.frame")$host.age)
ind.train <- sample(1:nrow(X), ceiling(0.7*nrow(X)))
ind.test <- setdiff(1:nrow(X), ind.train) 
X.train <- X[ind.train,]
X.test <- X[ind.test,]
y.train <- y[ind.train]
y.test <- y[ind.test]

CodeFolder <- "RCode/"
source(paste0(CodeFolder,"func_libs.R"))

today <- '20210701'
linkage <- 'average'   # controls the hclust

## ---- Set the error matrix ----
err <- matrix(NA,nrow=6,ncol=1)
rownames(err) <- c('supervised','propr','coat','pba','classo','spetral')
colnames(err) <- c('pred')

## Tuning parameters
nlam <- 100
maxlam <- 2*max(abs(crossprod(log(X[ind.train,]),y.train)/length(y.train)))
lambda <- exp(seq(from=log(maxlam), to=log(1e-4), length.out=nlam))

## Evaluate the correlation between y and each of the log ratio
## Construct SBP from supervised log-ratio matrix

S <- matrix(0,p,p)
rownames(S) <- colnames(S) <- colnames(X.train)
for (j in 1:(p-1)){
  for (k in (j+1):p){
    newx <- log(X.train[,j]) - log(X.train[,k])
    S[j,k] <- S[k,j] <- abs(cor(newx,y.train))^2
  }
}
h_supervised <- hclust(as.dist(1-S),method = linkage)
sbp_supervised <- sbp.fromHclust(h_supervised)
normalize <- function(contrast){
  if(any(!contrast %in% c(-1, 0, 1))) stop("Contrast must contain [-1, 0, 1] only.")

  lpos <- sum(contrast == 1)
  lneg <- sum(contrast == -1)
  const <- sqrt((lpos*lneg)/(lpos+lneg))
  contrast[contrast==1] = 1/lpos
  contrast[contrast==-1] = -1/lneg

  const * contrast
}
U <- apply(sbp_supervised,2,normalize)

ba_supervised <- balance.fromSBP(X,sbp_supervised)
fit_supervised <- run.glmnet(x=ba_supervised[ind.train,],y.train,xt=ba_supervised[ind.test,],y.test,lambda = lambda)
err[1,1] <- fit_supervised$mse.pred

a <- HSClust(S)
sbp_spec <- sbp.fromHSClust(a$allLevels)
rownames(sbp_spec) <- colnames(X)
ba_spec <- balance.fromSBP(X,sbp_spec)
fit_spec <- run.glmnet(x=ba_spec[ind.train,],y.train,xt=ba_spec[ind.test,],y.test,lambda = lambda)
err[6,1] <- fit_spec$mse.pred
# roc_supervised <- apply(fit_supervised$beta, 2, function(a) roc.for.coef.LR(a, beta_lc, sbp_supervised))

## Fit compositional Lasso using the same lambda
z <- log(X.train)
z_test <- log(X.test)
fit_classo <- cv.func('ConstrLasso',y=y.train,x=z,C=matrix(1,p,1), lambda = lambda, nfolds=10)
pred_classo <- fit_classo$int[which.min(fit_classo$cvm)] + z_test %*% fit_classo$bet[,which.min(fit_classo$cvm)]
err[5,1] <- mean((pred_classo - y.test)^2)
# roc_classo <- apply(fit_classo$bet, 2, function(a) roc.for.coef(a,beta_lc))

print(err)
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
