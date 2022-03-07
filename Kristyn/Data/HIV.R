# Date: 3/7/2022
rm(list=ls())

library(zCompositions)

source("Kristyn/Functions/supervisedlogratios.R")

################################################################################
# HIV: a data set in selbal package
#   n = 155 samples, 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 covariate (MSM), 
#   1 response (HIV_Status - binary)
################################################################################
W1 = selbal::HIV[, 1:60]
X1 = sweep(W1, 1, rowSums(W1), FUN='/')
y1 = ifelse(as.character(selbal::HIV[, 62]) == "Pos", 1, 0)

################################################################################
# Strategy 1: Replace 0's with 0.5 (used in Lin et al. 2014 [classo])
W1_0.5 = W1
W1_0.5[W1_0.5 == 0] = 0.5
# any(as.vector(W1_0.5) == 0)
X1_0.5 = sweep(W1_0.5, 1, rowSums(W1_0.5), FUN='/')

slrcor1_0.5 = getSlrMatrix(y = y1, X = X1_0.5, type = "similarity")
fields::image.plot(slrcor1_0.5)

################################################################################
# Strategy 2: GBM (used in Rivera-Pinto et al. 2018 [selbal])
X1_gbm = cmultRepl2(W1, zero.rep = "bayes")

slrcor1_gbm = getSlrMatrix(y = y1, X = X1_gbm, type = "similarity")
fields::image.plot(slrcor1_gbm)




################################################################################
################################################################################
################################################################################




################################################################################
# sCD14: a data set in selbal package
#   n = 151 samples (a subset from HIV data set), 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
################################################################################
W2 = selbal::sCD14[, 1:60]
X2 = sweep(W2, 1, rowSums(W2), FUN='/')
y2 = selbal::sCD14[, 61]

################################################################################
# Strategy 1: Replace 0's with 0.5 (used in Lin et al. 2014 [classo])
W2_0.5 = W2
W2_0.5[W2_0.5 == 0] = 0.5
# any(as.vector(W2_0.5) == 0)
X2_0.5 = sweep(W2_0.5, 1, rowSums(W2_0.5), FUN='/')

slrcor2_0.5 = getSlrMatrix(y = y2, X = X2_0.5, type = "similarity")
fields::image.plot(slrcor2_0.5)

################################################################################
# Strategy 2: GBM (used in Rivera-Pinto et al. 2018 [selbal])
X2_gbm = cmultRepl2(W2, zero.rep = "bayes")

slrcor2_gbm = getSlrMatrix(y = y2, X = X2_gbm, type = "similarity")
fields::image.plot(slrcor2_gbm)




