# Purpose: compare slr to other methods on data sets
# Date: 8/15/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs"

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(selbal)

source("RCode/func_libs.R")
source("slr_analyses/Functions/slrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

# tuning parameter settings
hparam = "1se"
K = 10
scaling = TRUE

filter.perc = 0.8 # 0.8, 1
split.perc = 0.7 # 0.7, 0.8

file.end = paste0(
  "/HIV",
  "_split", split.perc, 
  "_filter", filter.perc, 
  "_hparam", hparam, 
  "_gbm")

##############################################################################
# HIV: one of the HIV data sets in selbal package -- has binary response
#   n = 155 samples, 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 covariate (MSM), 
#   1 response (HIV_Status - binary)
W = selbal::HIV[, 1:60]
W.origin <- W
W <- W[,apply(W==0,2,mean)<filter.perc]
covar = data.frame(MSM = selbal::HIV[, 61])
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::HIV[, 62]
# levels(Y) # (control, case)
Y2 = ifelse(Y == "Pos", 1, 0)

################################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

################################################################################
# fit methods
################################################################################
p = ncol(X)

# classo #######################################################################
# if(hparam == "min"){
#   classo = codalasso(
#     X_gbm, Y2, numFolds = K, gamma = 0, type.measure = "AUC", stratify = FALSE)
# } else if(hparam == "1se"){
#   classo = codalasso(
#     X_gbm, Y2, numFolds = K, gamma = 1, type.measure = "AUC", stratify = FALSE)
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   classo,
#   paste0(
#     output_dir, file.end,
#     "_classo",
#     ".rds"))

cl = readRDS(
  paste0(
    output_dir, file.end,
    "_classo",
    ".rds"))

# slr - spectral ###############################################################
# slrspeccv = cv.slr(
#   x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "spectral",
#   response.type = "binary", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "auc",
#   scale = scaling, trace.it = FALSE)
# if(hparam == "min"){
#   slrspec = slr(
#     x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "spectral",
#     response.type = "binary", s0.perc = 0, zeta = 0,
#     threshold = slrspeccv$threshold[slrspeccv$index["min",]], 
#     positive.slope = TRUE)
# } else if(hparam == "1se"){
#   slrspec = slr(
#     x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "spectral",
#     response.type = "binary", s0.perc = 0, zeta = 0,
#     threshold = slrspeccv$threshold[slrspeccv$index["1se",]], 
#     positive.slope = TRUE)
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   slrspeccv,
#   paste0(
#     output_dir, file.end,
#     "_slrcv_spectral",
#     ".rds"))
# saveRDS(
#   slrspec,
#   paste0(
#     output_dir, file.end,
#     "_slr_spectral",
#     ".rds"))

slrspeccv = readRDS(
  paste0(
    output_dir, file.end,
    "_slrcv_spectral",
    ".rds"))
slrspec = readRDS(
  paste0(
    output_dir, file.end,
    "_slr_spectral",
    ".rds"))

# slr - hierarchical ###########################################################
# slrhiercv = cv.slr(
#   x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "binary", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "auc",
#   scale = scaling, trace.it = FALSE)
# if(hparam == "min"){
#   slrhier = slr(
#     x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "hierarchical",
#     response.type = "binary", s0.perc = 0, zeta = 0,
#     threshold = slrhiercv$threshold[slrhiercv$index["min",]], 
#     positive.slope = TRUE)
# } else if(hparam == "1se"){
#   slrhier = slr(
#     x = X_gbm, y = Y2, screen.method = "wald", cluster.method = "hierarchical",
#     response.type = "binary", s0.perc = 0, zeta = 0,
#     threshold = slrhiercv$threshold[slrhiercv$index["1se",]], 
#     positive.slope = TRUE)
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   slrhiercv,
#   paste0(
#     output_dir, file.end,
#     "_slrcv_hierarchical",
#     ".rds"))
# saveRDS(
#   slrhier,
#   paste0(
#     output_dir, file.end,
#     "_slr_hierarchical",
#     ".rds"))

slrhiercv = readRDS(
  paste0(
    output_dir, file.end,
    "_slrcv_hierarchical",
    ".rds"))
slrhier = readRDS(
  paste0(
    output_dir, file.end,
    "_slr_hierarchical",
    ".rds"))

# selbal #######################################################################
# if(hparam == "min"){
#   slbl = selbal::selbal.cv(x = X_gbm, y = Y, n.fold = K, opt.cri = "min")
# } else if(hparam == "1se"){
#   slbl = selbal::selbal.cv(x = X_gbm, y = Y, n.fold = K, opt.cri = "1se")
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   slbl,
#   paste0(
#     output_dir, file.end,
#     "_selbal",
#     ".rds"))

slbl = readRDS(
  paste0(
    output_dir, file.end,
    "_selbal",
    ".rds"))

# selbal - covariate ###########################################################
# if(hparam == "min"){
#   slbl2 = selbal::selbal.cv(
#     x = X_gbm, y = Y, n.fold = K, covar = covar, opt.cri = "min")
# } else if(hparam == "1se"){
#   slbl2 = selbal::selbal.cv(
#     x = X_gbm, y = Y, n.fold = K, covar = covar, opt.cri = "1se")
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   slbl2,
#   paste0(
#     output_dir, file.end,
#     "_selbal_covar",
#     ".rds"))

# slblc = readRDS(
#   paste0(
#     output_dir, file.end,
#     "_selbal_covar",
#     ".rds"))

# codacore #####################################################################
# library(codacore)
# if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
#   reticulate::use_condaenv("anaconda3")
# }
# if(hparam == "min"){
#   codacore0 = codacore::codacore(
#     x = X_gbm, y = Y2, logRatioType = "ILR",
#     objective = "binary classification", cvParams = list(numFolds = K), 
#     lambda = 0) 
# } else if(hparam == "1se"){
#   codacore0 = codacore::codacore(
#     x = X_gbm, y = Y2, logRatioType = "ILR",
#     objective = "binary classification", cvParams = list(numFolds = K), 
#     lambda = 1) 
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   codacore0,
#   paste0(
#     output_dir, file.end,
#     "_codacore",
#     ".rds"))

# cdcr = readRDS(
#   paste0(
#     output_dir, file.end,
#     "_codacore",
#     ".rds"))

# codacore - 1 balance #########################################################
# library(codacore)
# if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
#   reticulate::use_condaenv("anaconda3")
# }
# if(hparam == "min"){
#   codacore1 = codacore::codacore(
#     x = X_gbm, y = Y2, logRatioType = "ILR",
#     objective = "binary classification", cvParams = list(numFolds = K),
#     maxBaseLearners = 1,
#     lambda = 0)
# } else if(hparam == "1se"){
#   codacore1 = codacore::codacore(
#     x = X_gbm, y = Y2, logRatioType = "ILR",
#     objective = "binary classification", cvParams = list(numFolds = K),
#     maxBaseLearners = 1,
#     lambda = 1)
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   codacore1,
#   paste0(
#     output_dir, file.end,
#     "_codacore1",
#     ".rds"))

cdcr1 = readRDS(
  paste0(
    output_dir, file.end,
    "_codacore1",
    ".rds"))

# log-ratio lasso ############################################################
# library(logratiolasso)
# source("slr_analyses/Functions/logratiolasso.R")
# W_gbm.c = scale(log(X_gbm), center = TRUE, scale = FALSE)
# if(hparam == "min"){
#   lrl <- cv_two_stage(
#     z = W_gbm.c, y = Y2, n_folds = K, family="binomial", gamma = 0)
#   # lrl.betahat = lrl$beta_min
# } else if(hparam == "1se"){
#   lrl <- cv_two_stage(
#     z = W_gbm.c, y = Y2, n_folds = K, family="binomial", gamma = 1)
#   # lrl.betahat = lrl$beta_gammase
# } else{
#   stop("invalid hparam setting (method for selecting hyperparameter(s)).")
# }
# saveRDS(
#   lrl,
#   paste0(
#     output_dir, file.end,
#     "_lrlasso",
#     ".rds"))

lrl = readRDS(
  paste0(
    output_dir, file.end,
    "_lrlasso",
    ".rds"))

################################################################################
# get active sets and selected balances (if applicable)
################################################################################

# classo #######################################################################
# selected variables
cl.betahat = cl$cll$betas[-1]
# positive/negative effect on response
colnames(X)[cl.betahat > 0 & abs(cl.betahat) > 1e-8] # positive effect
colnames(X)[cl.betahat < 0 & abs(cl.betahat) > 1e-8] # negative effect
sum(abs(cl.betahat) > 1e-8)

# for overleaf
cat(paste(str_replace_all(
  colnames(X)[cl.betahat > 0 & abs(cl.betahat) > 1e-8],
  fixed("_"), "\\_"), collapse = ", "))
cat(paste(str_replace_all(
  colnames(X)[cl.betahat < 0 & abs(cl.betahat) > 1e-8],
  fixed("_"), "\\_"), collapse = ", "))

# slr - spectral ###############################################################
# SBP
slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrspec.fullSBP) = colnames(X)
slrspec.fullSBP[match(
  names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp
# thetahat
slrspec.coefs = getCoefsBM(
  coefs = coefficients(slrspec$fit), sbp = slrspec.fullSBP)
# numerator (I+) / denominator (I-) of selected balance
names(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs > 0]
names(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs < 0]
sum(slrspec.fullSBP[, 1] != 0)

# for overleaf
cat(paste(str_replace_all(
  names(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs > 0],
  fixed("_"), "\\_"), collapse = ", "))
cat(paste(str_replace_all(
  names(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs < 0],
  fixed("_"), "\\_"), collapse = ", "))

# # slr - hierarchical ###########################################################
# # SBP
# slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
# rownames(slrhier.fullSBP) = colnames(X)
# slrhier.fullSBP[match(
#   names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
# # thetahat
# slrhier.coefs = getCoefsBM(
#   coefs = coefficients(slrhier$fit), sbp = slrhier.fullSBP)
# # numerator (I+) / denominator (I-) of selected balance
# rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs > 0]
# rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs < 0]
# sum(slrhier.fullSBP[, 1] != 0)
# 
# # for overleaf
# cat(paste(str_replace_all(
#   rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs > 0],
#   fixed("_"), "\\_"), collapse = ", "))
# cat(paste(str_replace_all(
#   rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs < 0],
#   fixed("_"), "\\_"), collapse = ", "))

# selbal #######################################################################
# numerator (I+) / denominator (I-) of selected balance
slbl$global.balance[slbl$global.balance$Group == "NUM", "Taxa"] # 1
slbl$global.balance[slbl$global.balance$Group == "DEN", "Taxa"] # 1
length(slbl$global.balance$Group)

# for overleaf
cat(paste(str_replace_all(
  slbl$global.balance[slbl$global.balance$Group == "NUM", "Taxa"],
  fixed("_"), "\\_"), collapse = ", "))
cat(paste(str_replace_all(
  slbl$global.balance[slbl$global.balance$Group == "DEN", "Taxa"] ,
  fixed("_"), "\\_"), collapse = ", "))

# # selbal - covariate ###########################################################
# # numerator (I+) / denominator (I-) of selected balance
# slblc$global.balance[slblc$global.balance$Group == "NUM", "Taxa"] # 1
# slblc$global.balance[slblc$global.balance$Group == "DEN", "Taxa"] # 1

# # codacore #####################################################################
# length(cdcr$ensemble) # one balance selected
# # numerator (I+) / denominator (I-) of 1st selected balance
# cdcr$ensemble[[1]]$slope # positive (if negative, num -> den & vice versa)
# colnames(X)[cdcr$ensemble[[1]]$hard$numerator]
# colnames(X)[cdcr$ensemble[[1]]$hard$denominator]
# sum(c(cdcr$ensemble[[1]]$hard$numerator, cdcr$ensemble[[1]]$hard$denominator))
# 
# # for overleaf
# cat(paste(str_replace_all(
#   colnames(X)[cdcr$ensemble[[1]]$hard$numerator],
#   fixed("_"), "\\_"), collapse = ", "))
# cat(paste(str_replace_all(
#   colnames(X)[cdcr$ensemble[[1]]$hard$denominator],
#   fixed("_"), "\\_"), collapse = ", "))

# codacore1 ####################################################################
length(cdcr1$ensemble) # one balance selected
# numerator (I+) / denominator (I-) of selected balance
cdcr1$ensemble[[1]]$slope # positive (if negative, num -> den & vice versa)
colnames(X)[cdcr1$ensemble[[1]]$hard$numerator]
colnames(X)[cdcr1$ensemble[[1]]$hard$denominator]
sum(c(cdcr1$ensemble[[1]]$hard$numerator, cdcr1$ensemble[[1]]$hard$denominator))

# for overleaf
cat(paste(str_replace_all(
  colnames(X)[cdcr1$ensemble[[1]]$hard$numerator],
  fixed("_"), "\\_"), collapse = ", "))
cat(paste(str_replace_all(
  colnames(X)[cdcr1$ensemble[[1]]$hard$denominator],
  fixed("_"), "\\_"), collapse = ", "))

# log-ratio lasso ############################################################
# positive/negative effect on response
colnames(X)[lrl$beta_min > 0 & abs(lrl$beta_min) > 1e-8] # positive effect
colnames(X)[lrl$beta_min < 0 & abs(lrl$beta_min) > 1e-8] # negative effect
sum(abs(lrl$beta_min) > 1e-8)

# for overleaf
cat(paste(str_replace_all(
  colnames(X)[lrl$beta_min > 0 & abs(lrl$beta_min) > 1e-8],
  fixed("_"), "\\_"), collapse = ", "))
cat(paste(str_replace_all(
  colnames(X)[lrl$beta_min < 0 & abs(lrl$beta_min) > 1e-8] ,
  fixed("_"), "\\_"), collapse = ", "))














