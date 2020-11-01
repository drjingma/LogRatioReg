workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
setwd(workdir)

# Kristyn sources
source("Kristyn/Functions/classic_lasso.R")
source("Kristyn/Functions/compositional_lasso.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/coat.R")
source("Kristyn/Functions/principlebalances.R")

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

### trying coat
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

# trying principle balances
test_pbLASSO = fitPBLasso(X.prop.tr, y.tr, lambda = NULL)
pbLassofit = function(x){
  sbp_pba = pba(x)
  xb = sbp_pba@pba
  predict(test_pbLASSO$glmnet, newx = xb, type = "response")
}
pbLassoyhat = pbLassofit(X.prop.te)
pbLassoMSE = mean(pbLassoyhat - y.te)^2




# seeing proportion of 0 beta elts
comp.0beta = sum(betahatcompositional == 0) / length(betahatcompositional)
classic.0beta = sum(betahatclassic == 0) / length(betahatclassic)
slr.0beta = sum(test_slrLASSO$betahat == 0) / length(test_slrLASSO$betahat)
coat.0beta = sum(test_coatLASSO$betahat == 0) / length(test_coatLASSO$betahat)
pb.0beta = sum(test_pbLASSO$betahat == 0) / length(test_pbLASSO$betahat)
compare.0beta = data.frame(compositionalLasso = comp.0beta, 
                           classicLasso = classic.0beta, 
                           coatLasso = coat.0beta,
                           pbLasso = pb.0beta, 
                           slrLasso = slr.0beta)
compare.0beta

# compare MSEs of three methods
compare.mse = data.frame(compositionalLasso = compositionalLassoMSE, 
                         classicLasso = classicLassoMSE, 
                         coatLasso = coatLassoMSE,
                         PBLasso = pbLassoMSE, 
                         slrLasso = slrLassoMSE)
compare.mse


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
  predicted = as.vector(c(compositionalLassoyhat, classicLassoyhat, coatLassoyhat, slrLassoyhat)), 
  trueresponse = rep(y.te, 4), 
  type =  rep(c("compositionalLasso", "classicLasso", "coatLasso", "slrLasso"), 
              c(length(compositionalLassoyhat), length(classicLassoyhat), length(coatLassoyhat), length(slrLassoyhat)))
)
ggplot(ggdata2, aes(x = trueresponse, y = predicted, color = type)) + 
  facet_wrap(vars(type)) +
  geom_point() + 
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







