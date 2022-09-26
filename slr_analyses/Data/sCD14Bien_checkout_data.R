# Purpose: compare slr to other methods on data sets
# Date: 9/25/2022
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
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4

################################################################################
# sCD14: another HIV data set, this time given by Bien 2020 (trac paper)
#   n = 152 samples (a subset from sCD14 data set), 
#   p = 539 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
sCD14 = readRDS("Data/sCD14.RDS")
W = sCD14$x
X = sweep(W, 1, rowSums(W), FUN='/')
Y = sCD14$y

# check whether some variables are too rare
#   remove those that are in less than 25% of patients

n = nrow(W)
p = ncol(W)
dim(W)

W.nonzero = W != 0
num.samples = apply(W.nonzero, 2, sum)
prop.samples = num.samples / n

prop.df = data.frame(
  Taxa = names(prop.samples),
  Proportion = prop.samples
)
rownames(prop.df) = NULL
prop.plt = ggplot(data=prop.df, aes(x=Proportion, y=reorder(Taxa, Proportion))) +
  geom_bar(stat="identity")
prop.plt
##############################################################################
# 0-Handling -- GBM (used in Rivera-Pinto et al. 2018 [selbal])
X_gbm = selbal::cmultRepl2(W, zero.rep = "bayes")



