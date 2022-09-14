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
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4

################################################################################
# sCD14: another HIV data set in selbal package
#   n = 151 samples (a subset from sCD14 data set), 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
W = selbal::sCD14[, 1:60]
X = sweep(W, 1, rowSums(W), FUN='/')
Y = selbal::sCD14[, 61]

################################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

################################################################################
# fit methods
################################################################################
p = ncol(X)

# classo #######################################################################
# classo = cv.func(
#   method="ConstrLasso", y = Y, x = log(X_gbm), Cmat = matrix(1, p, 1),
#   nlam = nlam, nfolds = K, tol = tol, intercept = intercept,
#   scaling = scaling)
# saveRDS(
#   classo,
#   paste0(
#     output_dir, "/sCD14",
#     "_classo",
#     "_gbm",
#     ".rds"))

cl = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_classo",
    "_gbm",
    ".rds"))

# slr - spectral ###############################################################
# slrspeccv = cv.slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "mse",
#   scale = scaling, trace.it = FALSE)
# slrspec = slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "spectral",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   threshold = slrspeccv$threshold[slrspeccv$index["1se",]], 
#   positive.slope = TRUE)
# saveRDS(
#   slrspeccv,
#   paste0(
#     output_dir, "/sCD14",
#     "_slrcv_spectral",
#     "_gbm",
#     ".rds"))
# saveRDS(
#   slrspec,
#   paste0(
#     output_dir, "/sCD14",
#     "_slr_spectral",
#     "_gbm",
#     ".rds"))

slrspeccv = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slrcv_spectral",
    "_gbm",
    ".rds"))
slrspec = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slr_spectral",
    "_gbm",
    ".rds"))

# slr - hierarchical ###########################################################
# slrhiercv = cv.slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   nfolds = K, type.measure = "mse",
#   scale = scaling, trace.it = FALSE)
# slrhier = slr(
#   x = X_gbm, y = Y, screen.method = "wald", cluster.method = "hierarchical",
#   response.type = "continuous", s0.perc = 0, zeta = 0,
#   threshold = slrhiercv$threshold[slrhiercv$index["1se",]], 
#   positive.slope = TRUE)
# saveRDS(
#   slrhiercv,
#   paste0(
#     output_dir, "/sCD14",
#     "_slrcv_hierarchical",
#     "_gbm",
#     ".rds"))
# saveRDS(
#   slrhier,
#   paste0(
#     output_dir, "/sCD14",
#     "_slr_hierarchical",
#     "_gbm",
#     ".rds"))

slrhiercv = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slrcv_hierarchical",
    "_gbm",
    ".rds"))
slrhier = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_slr_hierarchical",
    "_gbm",
    ".rds"))

# selbal #######################################################################
# slbl = selbal::selbal.cv(x = X_gbm, y = Y, n.fold = K)
# saveRDS(
#   slbl,
#   paste0(
#     output_dir, "/sCD14",
#     "_selbal",
#     "_gbm",
#     ".rds"))

slbl = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_selbal",
    "_gbm",
    ".rds"))

# codacore #####################################################################
# library(codacore)
# if(getwd() == "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"){
#   reticulate::use_condaenv("anaconda3")
# }
# codacore0 = codacore::codacore(
#   x = X_gbm, y = Y, logRatioType = "ILR",
#   objective = "regression", cvParams = list(numFolds = K))
# saveRDS(
#   codacore0,
#   paste0(
#     output_dir, "/sCD14",
#     "_codacore",
#     "_gbm",
#     ".rds"))

cdcr = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_codacore",
    "_gbm",
    ".rds"))

# log-ratio lasso ##############################################################
# library(logratiolasso)
# source("slr_analyses/Functions/logratiolasso.R")
# W_gbm.c = scale(log(X_gbm), center = TRUE, scale = FALSE)
# lrl <- cv_two_stage(
#   z = W_gbm.c, y = Y, n_folds = K, family="gaussian")
# saveRDS(
#   lrl,
#   paste0(
#     output_dir, "/sCD14",
#     "_lrlasso",
#     "_gbm",
#     ".rds"))

lrl = readRDS(
  paste0(
    output_dir, "/sCD14",
    "_lrlasso",
    "_gbm",
    ".rds"))



################################################################################
# get active sets and selected balances (if applicable)
################################################################################

# classo #######################################################################
# selected variables
oneSErule = min(cl$cvm) + cl$cvsd[which.min(cl$cvm)] * 1
cl.lam.idx = which(cl$cvm <= oneSErule)[1]
cl.a0 = cl$int[cl.lam.idx]
cl.betahat = cl$bet[, cl.lam.idx]
# positive/negative effect on response
colnames(X)[cl.betahat > 0 & abs(cl.betahat) > 1e-8] # positive effect
colnames(X)[cl.betahat < 0 & abs(cl.betahat) > 1e-8] # negative effect
sum(abs(cl.betahat) > 1e-8)

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
rownames(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs > 0]
rownames(slrspec.coefs$llc.coefs)[slrspec.coefs$llc.coefs < 0]
sum(slrspec.fullSBP[, 1] != 0)

# slr - hierarchical ###########################################################
# SBP
slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrhier.fullSBP) = colnames(X)
slrhier.fullSBP[match(
  names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
# thetahat 
slrhier.coefs = getCoefsBM(
  coefs = coefficients(slrhier$fit), sbp = slrhier.fullSBP)
# numerator (I+) / denominator (I-) of selected balance
rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs > 0]
rownames(slrhier.coefs$llc.coefs)[slrhier.coefs$llc.coefs < 0]
sum(slrhier.fullSBP[, 1] != 0)

# selbal #######################################################################
# numerator (I+) / denominator (I-) of selected balance
slbl$global.balance[slbl$global.balance$Group == "NUM", "Taxa"] # 3
slbl$global.balance[slbl$global.balance$Group == "DEN", "Taxa"] # 2

# codacore #####################################################################
length(cdcr$ensemble) # one balance selected
# numerator (I+) / denominator (I-) of selected balance
cdcr$ensemble[[1]]$slope # positive (if negative, num -> den & vice versa)
colnames(X)[cdcr$ensemble[[1]]$hard$numerator]
colnames(X)[cdcr$ensemble[[1]]$hard$denominator]

# log-ratio lasso ############################################################
# positive/negative effect on response
colnames(X)[lrl$beta_min > 0 & abs(lrl$beta_min) > 1e-8] # positive effect
colnames(X)[lrl$beta_min < 0 & abs(lrl$beta_min) > 1e-8] # negative effect
sum(abs(lrl$beta_min) > 1e-8)



################################################################################
# make some plots for supervised log-ratios methods
################################################################################
library(tidyverse)

# slr - spectral ###############################################################

# type.measure across different threshold values
cv_data.spec = data.frame(
  Threshold = slrspeccv$threshold, 
  CVMetric = slrspeccv$cvm,
  Lower = slrspeccv$cvm - slrspeccv$cvsd, 
  Upper = slrspeccv$cvm + slrspeccv$cvsd
)
ggplot(cv_data.spec, aes(x = Threshold, y = CVMetric)) + 
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) + #, width = 0.2, col = "gray") +
  geom_vline(xintercept = slrspeccv$threshold.min, linetype = "dashed",
             col = "blue") +
  geom_vline(xintercept = slrspeccv$threshold.1se, linetype = "dotted",
           col = "lightblue") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-spec: threshold vs CV metric")

# plot the balance regression model with the data
bal_data.spec = data.frame(
  ilrX = slr.fromContrast(X_gbm, as.vector(slrspec.fullSBP)),
  y = Y
)
bal_data.spec
slrspec.coefs

lmfit = lm(y ~ ilrX, data = bal_data.spec)
summary(lmfit)


ggplot(bal_data.spec, aes(x = ilrX, y = Y)) + 
  geom_point() + 
  geom_abline(
    aes(
      intercept = slrspec.coefs$a0, 
      slope = slrspec.coefs$bm.coefs, 
      color = "slr-spec")) + 
  geom_abline(
    aes(
      intercept = cl.a0, 
      slope = 0, 
      color = "classo")) + 
  scale_color_manual(values = c(
    "slr-spec" = "black", "classo" = "blue")
    ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-spec: balance regression model")

# slr - hierarchical ###########################################################

# type.measure across different threshold values
cv_data.hier = data.frame(
  Threshold = slrhiercv$threshold, 
  CVMetric = slrhiercv$cvm,
  Lower = slrhiercv$cvm - slrhiercv$cvsd, 
  Upper = slrhiercv$cvm + slrhiercv$cvsd
)
ggplot(cv_data.hier, aes(x = Threshold, y = CVMetric)) + 
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) + #, width = 0.2, col = "gray") +
  geom_vline(xintercept = slrhiercv$threshold.min, linetype = "dashed",
             col = "blue") +
  geom_vline(xintercept = slrhiercv$threshold.1se, linetype = "dotted",
             col = "lightblue") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-hier: threshold vs CV metric")

# plot the balance regression model with the data
bal_data.hier = data.frame(
  ilrX = slr.fromContrast(X_gbm, as.vector(slrhier.fullSBP)),
  y = Y
)
bal_data.hier
slrhier.coefs
ggplot(bal_data.hier, aes(x = ilrX, y = Y)) + 
  geom_point() + 
  geom_abline(
    aes(
      intercept = slrhier.coefs$a0, 
      slope = slrhier.coefs$bm.coefs, 
      color = "slr-hier")) + 
  geom_abline(
    aes(
      intercept = cl.a0, 
      slope = 0, 
      color = "classo")) + 
  scale_color_manual(values = c(
    "slr-hier" = "black", "classo" = "blue")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-hier: balance regression model")

################################################################################
# plots for manuscript?
################################################################################
library(tidyverse)


# slr - spectral ###############################################################

# type.measure across different threshold values
cv_data.spec = data.frame(
  Threshold = slrspeccv$threshold, 
  CVMetric = slrspeccv$cvm,
  Lower = slrspeccv$cvm - slrspeccv$cvsd, 
  Upper = slrspeccv$cvm + slrspeccv$cvsd
)
ggplot(cv_data.spec, aes(x = Threshold, y = CVMetric)) + 
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) + #, width = 0.2, col = "gray") +
  geom_vline(xintercept = slrspeccv$threshold.min, linetype = "dashed",
             col = "blue") +
  geom_vline(xintercept = slrspeccv$threshold.1se, linetype = "dotted",
             col = "lightblue") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-spec: threshold vs CV metric")

# plot the balance regression model with the data
bal_data.spec = data.frame(
  ilrX = slr.fromContrast(X_gbm, as.vector(slrspec.fullSBP)),
  y = Y
)
bal_data.spec
slrspec.coefs
ggplot(bal_data.spec, aes(x = ilrX, y = Y)) + 
  geom_point() + 
  geom_abline(
    aes(
      intercept = slrspec.coefs$a0, 
      slope = slrspec.coefs$bm.coefs, 
      color = "slr-spec")) + 
  geom_abline(
    aes(
      intercept = cl.a0, 
      slope = 0, 
      color = "classo")) + 
  scale_color_manual(values = c(
    "slr-spec" = "black", "classo" = "blue")
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-spec: balance regression model")

# slr - hierarchical ###########################################################

# type.measure across different threshold values
cv_data.hier = data.frame(
  Threshold = slrhiercv$threshold, 
  CVMetric = slrhiercv$cvm,
  Lower = slrhiercv$cvm - slrhiercv$cvsd, 
  Upper = slrhiercv$cvm + slrhiercv$cvsd
)
ggplot(cv_data.hier, aes(x = Threshold, y = CVMetric)) + 
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) + #, width = 0.2, col = "gray") +
  geom_vline(xintercept = slrhiercv$threshold.min, linetype = "dashed",
             col = "blue") +
  geom_vline(xintercept = slrhiercv$threshold.1se, linetype = "dotted",
             col = "lightblue") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-hier: threshold vs CV metric")

# plot the balance regression model with the data
bal_data.hier = data.frame(
  ilrX = slr.fromContrast(X_gbm, as.vector(slrhier.fullSBP)),
  y = Y
)
bal_data.hier
slrhier.coefs
ggplot(bal_data.hier, aes(x = ilrX, y = Y)) + 
  geom_point() + 
  geom_abline(
    aes(
      intercept = slrhier.coefs$a0, 
      slope = slrhier.coefs$bm.coefs)) + 
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  ggtitle("slr-hier: balance regression model")

# classo #######################################################################

cl.data = data.frame(
  "index" = 1:length(cl$lambda),
  "value" = cl$lambda,
  "CVmean" = cl$cvm,
  "CVse" = cl$cvsd
)
cl.data.mlt = pivot_longer(
  cl.data, 
  cols = c("index", "value"),
  names_to = "lambda", 
  values_to = "value"
)
cl.idx.pts = data.frame(
  "x" = c(which.min(cl$cvm), which(cl$cvm <= oneSErule)[1])
)
cl.idx.pts$y = cl.data$CVmean[cl.idx.pts$x]
cl.idx.pts$lambda = "index"
cl.val.pts = data.frame(
  "x" = cl$lambda[cl.idx.pts$x], 
  "y" = cl.idx.pts$y, 
  "lambda" = "value"
)
cl.pts = rbind(cl.idx.pts, cl.val.pts)

cl_plt = ggplot(cl.data.mlt, aes(x = value)) + 
  geom_ribbon(aes(ymin = CVmean - CVse, ymax = CVmean + CVse), fill = "grey70") +
  geom_line(aes(y = CVmean)) + 
  geom_point(data = cl.pts, mapping = aes(x = x, y = y)) + 
  facet_wrap(vars(lambda), scales = "free")
ggsave(
  filename = paste0(
    "20220906",
    "_", "sCD14", "_classo.png"),
  plot = cl_plt,
  width = 6, height = 3, units = c("in")
)


