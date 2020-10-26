workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
setwd(workdir)

# Kristyn sources
source("Kristyn/Functions/classic_lasso.R")
source("Kristyn/Functions/compositional_lasso.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/coat.R")

# Dr. Ma sources
source("RCode/func_libs.R")
source("COAT-master/coat.R")

# libraries
library(mvtnorm)
library(balance)
library(microbenchmark)
library(ggplot2)
image_path = "/home/kristyn/Pictures"

# data
DataFolder <- "/Data/"
load(paste0(workdir, DataFolder, "BMI.rda"))
dim(raw_data) # 98 x 89
dim(X) # 98 x 87
dim(X.prop) # 98 x 87

# separate into training and test sets
# Split the data into K folds
n = dim(X)[1]
set.seed(1995)
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

test_slrLASSO = fitSLRLasso(X.prop.tr, y.tr, linkage = "average")
# warning msg is about computing separate statistics for each fold

# possible explanation for needing noise due to equal columns:
#   splitting the data on already sparse bacteria
set.seed(1) # just in case
test_slrLASSO = fitSLRLasso(X.prop.tr, y.tr, linkage = "average")
# there is no error when we don't split.

# possible solution: take out taxa that have too many 0's
prop_zero_per_sample = as.vector(apply(X, 1, function(a) sum(a == 0.5) / ncol(X)))
prop_zero_per_otu = as.vector(apply(X, 2, function(a) sum(a == 0.5) / nrow(X)))
X2 = X[, prop_zero_per_otu < 0.97] # take out OTUs with less than 5% nonzero
X.tr2 = X2[id_traintest == 1, ]
X.prop.tr2 = sweep(X.tr2, MARGIN = 1, STATS = rowSums(X.tr2), FUN = "/") # may not be necessary
X.te2 = X2[id_traintest == 2, ]
X.prop.te2 = sweep(X.te2, MARGIN = 1, STATS = rowSums(X.te2), FUN = "/") #  may not be necessary
y.tr2 = y.tr # y's don't change
y.te2 = y.te
# retrying with new data
test_slrLASSO = fitSLRLasso(X.prop.tr2, y.tr2, linkage = "average")

# fit to test data
slrLassofit = function(x){
  # names(x) = test_slrLASSO$btree$labels
  xb = computeBalances(test_slrLASSO$btree, x)
  xb %*% test_slrLASSO$betahat
}
slrLassoyhat = slrLassofit(X.prop.te2)
slrLassoMSE = mean(slrLassoyhat - y.te2)^2

# plot dendrogram
slr_btree = test_slrLASSO$btree
slr_btree$labels = 1:length(slr_btree$labels)
plot(slr_btree)

# retrying with linkage = 'complete'
test_slrLASSO2 = fitSLRLasso(X.prop.tr2, y.tr2, linkage = "complete")
slrLassofit2 = function(x){
  # names(x) = test_slrLASSO$btree$labels
  xb = computeBalances(test_slrLASSO2$btree, x)
  xb %*% test_slrLASSO2$betahat
}
slrLassoyhat2 = slrLassofit2(X.prop.te2)
slrLassoMSE2 = mean(slrLassoyhat2 - y.te2)^2

# plot dendrogram
slr_btree2 = test_slrLASSO2$btree
slr_btree2$labels = 1:length(slr_btree2$labels)
plot(slr_btree2)

### rerun classic and compositional lassos
# compositional lasso #
# betahat from compositional Lasso
test_cvCompLASSO = cvCompositionalLASSO(log(X.prop.tr2) ,y.tr2, lambda_seq = NULL, n_lambda = 30, k = 10)
# plot(test_cvCompLASSO$cvm ~ test_cvCompLASSO$lambda_seq, type = "l")
betahatcompositional = test_cvCompLASSO$beta_mat[, test_cvCompLASSO$cvm_idx]
betahat0compositional = test_cvCompLASSO$beta0_vec[test_cvCompLASSO$cvm_idx]
# check constraint
sum(betahatcompositional)
compositionalLassofit = function(x){
  betahat0compositional + log(x) %*% betahatcompositional
}
compositionalLassoyhat = apply(X.prop.te2, 1, compositionalLassofit)
compositionalLassoMSE = mean(compositionalLassoyhat - y.te2)^2
# classic lasso #
# betahat from classic Lasso
test_cvClassicLASSO = cvLASSO(log(X.prop.tr2) ,y.tr2, lambda_seq = NULL, n_lambda = 30, k = 10)
# plot(test_cvClassicLASSO$cvm ~ test_cvClassicLASSO$lambda_seq, type = "l")
betahatclassic = test_cvClassicLASSO$beta_mat[, test_cvClassicLASSO$cvm_idx]
betahat0classic = test_cvClassicLASSO$beta0_vec[test_cvClassicLASSO$cvm_idx]
# consider constraint
sum(betahatclassic)
classicLassofit = function(x){
  betahat0classic + log(x) %*% betahatclassic
}
classicLassoyhat = apply(X.prop.te2, 1, classicLassofit)
classicLassoMSE = mean(classicLassoyhat - y.te2)^2
# compare classic lasso and compositional lasso
sum(betahatcompositional)
sum(betahatclassic)

### trying coat
test_coatLASSO = fitCOATLasso(X.prop.tr2, y.tr2, linkage = "average")
coatLassofit = function(x){
  # names(x) = test_slrLASSO$btree$labels
  xb = computeBalances(test_coatLASSO$btree, x)
  xb %*% test_coatLASSO$betahat
}
coatLassoyhat = coatLassofit(X.prop.te2)
coatLassoMSE = mean(coatLassoyhat - y.te2)^2

# plot dendrogram
coat_btree = test_coatLASSO$btree
coat_btree$labels = 1:length(coat_btree$labels)
plot(coat_btree)

# trying principle balances
test_pbLASSO = fitPBLasso(X.prop.tr2, y.tr2, lambda = NULL)




# seeing proportion of 0 beta elts
comp.0beta = sum(betahatcompositional == 0) / length(betahatcompositional)
classic.0beta = sum(betahatclassic == 0) / length(betahatclassic)
slr.avg.0beta = sum(test_slrLASSO$betahat == 0) / length(test_slrLASSO$betahat)
slr.comp.0beta = sum(test_slrLASSO2$betahat == 0) / length(test_slrLASSO2$betahat)
coat.0beta = sum(test_coatLASSO$betahat == 0) / length(test_coatLASSO$betahat)
compare.0beta = data.frame(compositionalLasso = comp.0beta, 
                           classicLasso = classic.0beta, 
                           coatLasso = coat.0beta,
                           slrLassoAvg = slr.avg.0beta, 
                           slrLassoComp = slr.comp.0beta)
compare.0beta

# compare MSEs of three methods
compare.mse = data.frame(compositionalLasso = compositionalLassoMSE, 
                         classicLasso = classicLassoMSE, 
                         coatLasso = coatLassoMSE,
                         slrLassoAvg = slrLassoMSE, 
                         slrLassoComp = slrLassoMSE2)
compare.mse









ggdata = rbind(data.frame(var = prop_zero_per_otu, type = "otu"), 
               data.frame(var = prop_zero_per_sample, type = "sample"))
ggbreaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)
ggplot(ggdata, aes(y = var, x = type, fill = type)) + geom_boxplot() + theme_bw() + 
  scale_y_continuous("proportion of zeros", breaks = ggbreaks, labels = ggbreaks) + 
  theme(panel.grid.minor = element_blank())
ggsave("propzeros102620.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       scale = 1,
       width = 2,
       height = 4,
       units = c("in")
)               
z