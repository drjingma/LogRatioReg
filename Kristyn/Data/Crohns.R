# Date: 3/7/2022
rm(list=ls())

library(zCompositions)

source("Kristyn/Functions/supervisedlogratios.R")

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

################################################################################
# Strategy 2: GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = cmultRepl2(W, zero.rep = "bayes")

slrcor_gbm = getSlrMatrix(y = y, X = X_gbm, type = "similarity")
fields::image.plot(slrcor_gbm)

