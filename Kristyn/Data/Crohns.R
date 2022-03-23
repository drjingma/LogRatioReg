# Date: 3/13/2022
rm(list=ls())

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)
library(propr)

library(pROC)

source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/supervisedlogratioseta.R")
source("Kristyn/Functions/HSClust.R")
source("Kristyn/Functions/slrnew.R")
source("Kristyn/Functions/codalasso.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/helper_functions.R")

library(selbal)
library(zCompositions)
library(gplots)
library(pheatmap)

# Crohn: a data set in selbal package
#   n = 975 samples, 
#   p = 48 taxa (counts for microbial taxa at genus level), 
#   1 response (y - binary)
W = selbal::Crohn[, 1:48]
X = sweep(W, 1, rowSums(W), FUN='/')
y = ifelse(selbal::Crohn[, 49] == "CD", 1, 0)

# getSlrMatrix(y = y, X = X, type = "similarity")
# need to deal with 0's in X, first -- log(Xij) = Inf if Xij = 0.

################################################################################
# Strategy 1: Replace 0's with 0.5 (used in Lin et al. 2014 [classo])
W_0.5 = W
W_0.5[W_0.5 == 0] = 0.5
# any(as.vector(W_0.5) == 0)
X_0.5 = sweep(W_0.5, 1, rowSums(W_0.5), FUN='/')

slrcor_0.5 = getSlrMatrix(y = y, X = X_0.5, type = "similarity")
fields::image.plot(slrcor_0.5)

# heatmap(slrcor_0.5, col = rainbow(256),scale = "none", symm = TRUE)
heatmap.2(
  slrcor_0.5, 
  # col = rainbow(256), 
  scale = "none", symm = TRUE, trace = "none", 
  # key = FALSE, 
  cexRow = 0.5, cexCol = 0.5, offsetRow = -0.25, offsetCol = -0.25, 
  margins = c(0.25, 0.25))

pheatmap(slrcor_0.5)

slrfit = slr(x = X_0.5, y = y, classification = TRUE, rank1approx = TRUE)

################################################################################
# Strategy 2: GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

slrcor_gbm = getSlrMatrix(y = y, X = X_gbm, type = "similarity")
fields::image.plot(slrcor_gbm)

# heatmap(slrcor_gbm, col = rainbow(256),scale = "none", symm = TRUE)
heatmap.2(
  slrcor_gbm, col = rainbow(256), scale = "none", symm = TRUE, trace = "none", 
  key = FALSE, cexRow = 0.5, cexCol = 0.5, offsetRow = -0.25, offsetCol = -0.25, 
  margins = c(0.25, 0.25))
