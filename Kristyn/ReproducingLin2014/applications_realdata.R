workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
setwd(workdir)

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "classic_lasso.R"))
source(paste0(functions_path, "compositional_lasso.R"))
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "coat.R"))
source(paste0(functions_path, "principlebalances.R"))
source(paste0(functions_path, "selbal.R"))

# Dr. Ma sources
source("RCode/func_libs.R")
source("COAT-master/coat.R")

# libraries
library(mvtnorm)
library(balance)
# devtools::install_github(repo = "UVic-omics/selbal")
library(selbal)
library(microbenchmark)
library(ggplot2)
library(logratiolasso) # bates & tibshirani 2019
image_path = "/home/kristyn/Pictures"

# data
DataFolder <- "/Data/"
load(paste0(workdir, DataFolder, "BMI.rda"))
dim(raw_data) # 98 x 89
dim(X) # 98 x 87
dim(X.prop) # 98 x 87

# preprocessing: take out rare taxa/OTU, to avoid 0-vector columns when splitting
#   into training and test sets
prop_zero_per_sample = as.vector(apply(X, 1, function(a) sum(a == 0.5) / ncol(X)))
prop_zero_per_otu = as.vector(apply(X, 2, function(a) sum(a == 0.5) / nrow(X)))
X = X[, prop_zero_per_otu < 0.95] # take out OTUs with less than 5% nonzero
X.prop = X.prop[, prop_zero_per_otu < 0.97]

# separate into training and test sets
# Split the data into K folds
n = dim(X)[1]
set.seed(2)
shuffle = sample(1:n)
id_traintest = (shuffle %% 2) + 1
n_traintest = as.vector(table(id_traintest))
X.tr = X[id_traintest == 1, ]
X.prop.tr = sweep(X.tr, MARGIN = 1, STATS = rowSums(X.tr), FUN = "/")
X.te = X[id_traintest == 2, ]
X.prop.te = sweep(X.te, MARGIN = 1, STATS = rowSums(X.te), FUN = "/")
y.tr = y[id_traintest == 1]
y.te = y[id_traintest == 2]

################################################################################
#  compositional lasso
################################################################################

# y = log(X) beta + epsilon
# subject to constraint sum_{i = 1}^p beta_i == 1

# betahat from compositional Lasso
test_cvCompLASSO = cvCompositionalLASSO(log(X.prop.tr) ,y.tr, lambda_seq = NULL, n_lambda = 30, k = 10)
# plot(test_cvCompLASSO$cvm ~ test_cvCompLASSO$lambda_seq, type = "l")
betahatcompositional = test_cvCompLASSO$beta_mat[, test_cvCompLASSO$cvm_idx]
betahat0compositional = test_cvCompLASSO$beta0_vec[test_cvCompLASSO$cvm_idx]

# check constraint
sum(betahatcompositional)

compositionalLassofit = function(x){
  betahat0compositional + log(x) %*% betahatcompositional
}

compositionalLassoyhat = apply(X.prop.te, 1, compositionalLassofit)
compositionalLassoMSE = mean(compositionalLassoyhat - y.te)^2

################################################################################
#  subcompositional lasso
################################################################################

# y = log(X) beta + epsilon
# subject to constraint sum_{i = 1}^p beta_i == 1

# betahat from compositional Lasso
test_cvSubCompLASSO =  cv.func(
  method="ConstrLasso", y.tr, log(X.prop.tr), Cmat=matrix(1, dim(X.prop.tr)[2], 1), 
  lambda=NULL, nlam=50, intercept=TRUE, scaling=TRUE, nfolds=10, maxiter=1000, 
  tol=1e-4, seed=0)

betahatsubcomp = test_cvSubCompLASSO$bet[, which.min(test_cvSubCompLASSO$cvm)]
betahat0subcomp = test_cvSubCompLASSO$int[which.min(test_cvSubCompLASSO$cvm)]

# check constraint
sum(betahatsubcomp) # very small

subcompLassofit = function(x){
  betahat0subcomp + log(x) %*% betahatsubcomp
}

subcompLassoyhat = apply(X.prop.te, 1, subcompLassofit)
subcompLassoMSE = mean(subcompLassoyhat - y.te)^2

# compare to my compositional lasso
print(cbind(betahatsubcomp, betahatcompositional))
sum((betahatsubcomp - betahatcompositional)^2)
which(betahatsubcomp != 0)
which(betahatcompositional != 0)

# what happens when they use the same lambda sequence?
lambda_seq = test_cvSubCompLASSO$lambda
test_cvCompLASSO2 = cvCompositionalLASSO(log(X.prop.tr) ,y.tr, lambda_seq = lambda_seq, k = 10)
betahatcomp2 = test_cvCompLASSO2$beta_mat[, test_cvCompLASSO2$cvm_idx]
betahat0comp2 = test_cvCompLASSO2$beta0_vec[test_cvCompLASSO2$cvm_idx]
# check constraint
sum(betahatcomp2)
compLassofit2 = function(x){
  betahat0comp2 + log(x) %*% betahatcomp2
}
compLassoyhat2 = apply(X.prop.te, 1, compLassofit2)
compLassoMSE2 = mean(compLassoyhat2 - y.te)^2
# compare to my compositional lasso
print(cbind(betahatsubcomp, betahatcomp2))
sum((betahatsubcomp - betahatcomp2)^2)
which(betahatsubcomp != 0)
which(betahatcomp2 != 0)
lambda_min_subcomp = test_cvSubCompLASSO$lambda[which.min(test_cvSubCompLASSO$cvm)]
lambda_min_comp2 = test_cvCompLASSO2$lambda_seq[test_cvCompLASSO2$cvm_idx]

# what happens when they use the same lambda value?
lambda_val = lambda_min_subcomp
# my version
test_cvCompLASSO.lam = fitCompositionalLASSO(log(X.prop.tr) ,y.tr, lambda_val)
betahatcomp.lam = test_cvCompLASSO.lam$beta_mat[, 1]
betahat0complam = test_cvCompLASSO.lam$beta0_vec
# check constraint
sum(betahatcomp.lam)
# implemented version
test_cvSubCompLASSO.lam =  ConstrLasso(
  y.tr, log(X.prop.tr), Cmat = matrix(1, dim(log(X.prop.tr))[2], 1), 
  lambda = lambda_val, nlam = 1, intercept=TRUE, scaling=TRUE, maxiter=1000, 
  tol=1e-4)
betahatsubcomp.lam = test_cvSubCompLASSO.lam$bet
betahat0subcomp.lam = test_cvSubCompLASSO.lam$int
# compare
# check constraint
sum(betahatcomp.lam)
sum(betahatsubcomp.lam)
# check lambdas are actually same
test_cvCompLASSO.lam$lambda_seq
test_cvSubCompLASSO.lam$lambda
# compare betas
print(cbind(betahatsubcomp.lam, betahatcomp.lam))
sum((betahatsubcomp.lam - betahatcomp.lam)^2)
which(betahatsubcomp.lam != 0)
which(betahatcomp.lam != 0)
lambda_min_subcomp = test_cvSubCompLASSO$lambda[which.min(test_cvSubCompLASSO$cvm)]
lambda_min_comp2 = test_cvCompLASSO2$lambda_seq[test_cvCompLASSO2$cvm_idx]

################################################################################
# classic lasso
################################################################################

# betahat from classic Lasso
test_cvClassicLASSO = cvLASSO(log(X.prop.tr) ,y.tr, lambda_seq = NULL, n_lambda = 30, k = 10)
# plot(test_cvClassicLASSO$cvm ~ test_cvClassicLASSO$lambda_seq, type = "l")
betahatclassic = test_cvClassicLASSO$beta_mat[, test_cvClassicLASSO$cvm_idx]
betahat0classic = test_cvClassicLASSO$beta0_vec[test_cvClassicLASSO$cvm_idx]

# consider constraint
sum(betahatclassic)

classicLassofit = function(x){
  betahat0classic + log(x) %*% betahatclassic
}

classicLassoyhat = apply(X.prop.te, 1, classicLassofit)
classicLassoMSE = mean(classicLassoyhat - y.te)^2

# compare classic lasso and compositional lasso
sum(betahatcompositional)
sum(betahatclassic)

################################################################################
# supervised log-ratios lasso
################################################################################

test_slrLASSO = fitSLRLasso(X.prop.tr, y.tr, linkage = "complete") # works now!
# this is why we need to take out the rare OTUs -- so that there is lower chance
# of splitting train and test set and ending up with 0 columns in train set

# fit to test data
slrLassofit = function(x){
  xb = computeBalances(test_slrLASSO$btree, x)
  predict(test_slrLASSO$glmnet, newx = xb, type = "response")
}
slrLassoyhat = slrLassofit(X.prop.te)
slrLassoMSE = mean(slrLassoyhat - y.te)^2

# plot dendrogram
slr_btree = test_slrLASSO$btree
slr_btree$labels = 1:length(slr_btree$labels)
# plot(slr_btree)

################################################################################
#  coat
################################################################################

test_coatLASSO = fitCOATLasso(X.prop.tr, y.tr, linkage = "complete")
coatLassofit = function(x){
  xb = computeBalances(test_coatLASSO$btree, x)
  predict(test_coatLASSO$glmnet, newx = xb, type = "response")
}
coatLassoyhat = coatLassofit(X.prop.te)
coatLassoMSE = mean(coatLassoyhat - y.te)^2

# plot dendrogram
coat_btree = test_coatLASSO$btree
coat_btree$labels = 1:length(coat_btree$labels)
# plot(coat_btree)

################################################################################
#  selbal
################################################################################

selbalfit = fitselbal(X.prop.tr, y.tr)
selbalyhat = selbalfit(X.prop.te)
selbalMSE = mean(selbalyhat - y.te)^2

################################################################################
#  principle balances
################################################################################

# trying principle balances
test_pbLASSO = fitPBLasso(X.prop.tr, y.tr, lambda = NULL)
pbLassofit = function(x){
  sbp_pba = pba(x)
  xb = sbp_pba@pba
  predict(test_pbLASSO$glmnet, newx = xb, type = "response")
}
pbLassoyhat = pbLassofit(X.prop.te)
pbLassoMSE = mean(pbLassoyhat - y.te)^2

################################################################################
#  comparing methods' results
################################################################################


# seeing proportion of 0 beta elts
comp.0beta = sum(betahatcompositional == 0) / length(betahatcompositional)
classic.0beta = sum(betahatclassic == 0) / length(betahatclassic)
slr.0beta = sum(test_slrLASSO$betahat == 0) / length(test_slrLASSO$betahat)
coat.0beta = sum(test_coatLASSO$betahat == 0) / length(test_coatLASSO$betahat)
pb.0beta = sum(test_pbLASSO$betahat == 0) / length(test_pbLASSO$betahat)
propr.0beta = sum(test_proprLASSO$betahat == 0) / length(test_proprLASSO$betahat)
# selbal?
compare.0beta = data.frame(compositionalLasso = comp.0beta, 
                           classicLasso = classic.0beta, 
                           coat = coat.0beta,
                           pb = pb.0beta, 
                           propr = propr.0beta, 
                           slr = slr.0beta)
round(compare.0beta, 3)

# compare MSEs of three methods
compare.mse = data.frame(compositionalLasso = compositionalLassoMSE, 
                         classicLasso = classicLassoMSE, 
                         coat = coatLassoMSE,
                         principle = pbLassoMSE, 
                         selbal = selbalMSE, 
                         propr = proprLassoMSE, 
                         supervised = slrLassoMSE)
round(compare.mse, 3)

# 5-fold CV to compute MSE
#
#
#
#
#
#
#
#
#


### plot proportion of zeros in OTUs and samples
ggdata = rbind(data.frame(var = prop_zero_per_otu, type = "otu"), 
               data.frame(var = prop_zero_per_sample, type = "sample"))
ggbreaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)
ggplot(ggdata, aes(y = var, x = type, fill = type)) + geom_boxplot() + theme_bw() + 
  scale_y_continuous("proportion of zeros", breaks = ggbreaks, labels = ggbreaks) + 
  theme(panel.grid.minor = element_blank())
# ggsave("propzeros102620.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 2,
#        height = 4,
#        units = c("in")
# )               

# plot predicted response against true response
ggdata2 = data.frame(
  predicted = as.vector(c(
    compositionalLassoyhat, 
    classicLassoyhat, 
    coatLassoyhat, 
    pbLassoyhat, 
    selbalyhat, 
    proprLassoyhat, 
    slrLassoyhat)), 
  trueresponse = rep(y.te, 7), 
  type =  rep(c(
    "CompLasso", "Lasso", "COAT", "Principle","selbal", "propr", "Supervised"), 
    c(rep(length(y), 7)))
)
ggplot(ggdata2, aes(x = trueresponse, y = predicted, color = type)) + 
  facet_wrap(vars(type)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, alpha = 0.25) + 
  theme_bw()
# ggsave("predvtruey110220.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 6,
#        height = 4,
#        units = c("in")
# )