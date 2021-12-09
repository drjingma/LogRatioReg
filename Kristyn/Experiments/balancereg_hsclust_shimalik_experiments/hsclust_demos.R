# Purpose: demonstrate hierarchical spectral clustering 
# Date: 11/23/2021

################################################################################
# libraries and settings

library(limSolve)
library(mvtnorm)
library(Matrix)
library(glmnet)
library(compositions)
library(stats)

library(balance) # for sbp.fromHclust()
library(selbal)
library(propr)
library(sClust)
source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/HSClust.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/simulatedata.R")

# Sigma with 2 blocks ##########################################################

# Settings to toggle with
sigma.settings = "2blockSigma"
rho.type = "square" # 1 = "absolute value", 2 = "square"
theta.settings = "dense" # block, dense
# "block2" => choose j corresp. to block of correlated variables in block 2
# "dense" => j = 1
values.theta = 1
linkage = "average"
tol = 1e-4
nlam = 100
intercept = TRUE
K = 10
n = 100 #100
p = 200 #200
cor_ij = 0.2 # 0.2, 0.5
scaling = TRUE

# Population parameters
sigma_eps = 0.1 # 0.1, 0.5
num.blocks = 2
SigmaWblock = matrix(cor_ij, p / num.blocks, p / num.blocks)
for(i in 1:nrow(SigmaWblock)) SigmaWblock[i, i] = 1
SigmaW = as.matrix(bdiag(SigmaWblock, SigmaWblock))

# using hierarchical clustering
SigmaWtree = hclust(as.dist(1 - SigmaW), method = linkage)
sbp.SigmaW.hclust = sbp.fromHclust(SigmaWtree)
# rownames(sbp.SigmaW.hclust) = colnames(sbp.SigmaW.hclust) = NULL

# using hierarchical spectral clustering
hsclust.SigmaW = HSClust(
  W = getSimilarityMatrix(unnormalized_similarity_matrix = SigmaW), 
  levelMax = p - 1)
sbp.SigmaW.hsclust = sbp.fromHSClust(levels_matrix = hsclust.SigmaW$allLevels)

# inspect
all.equal(sbp.sort(sbp.SigmaW.hclust), sbp.sort(sbp.SigmaW.hsclust))
sbp.SigmaW.hclust # seems to leave out one variable at a time in binary partition
sbp.SigmaW.hsclust # partition produces a more symmetric tree

U = getU(sbp = sbp.SigmaW.hsclust) # transformation matrix

# theta settings
indices.theta = 1
theta = rep(0, p - 1)
values.theta = 1
values.theta = rep(values.theta, length(indices.theta))
theta[indices.theta] = values.theta
theta = as.matrix(theta)

# about beta
beta = as.vector(getBeta(theta, U = U))
names(beta) <- paste0('s', 1:p)
non0.beta = (beta != 0)
slr.non0.beta = abs(beta) > 10e-8
is0.beta = abs(beta) <= 10e-8
bspars = sum(non0.beta)

# Population parameters, continued
muW = c(rep(log(p), 5), rep(0, p - 5))
names(muW) = names(beta)

##############################################################################
# simulate data
set.seed(123) ############################################################################### take this part out!!!
fake.data = simulateBalanceReg(
  mu = muW, Sigma = SigmaW, U = U, n = 2 * n, theta = theta, 
  sigma.noise = sigma_eps)
colnames(fake.data$X) = names(beta)

# subset out training and test sets
X = fake.data$X[1:n, ]
X.test = fake.data$X[-(1:n), ]
Y <- fake.data$y[1:n, , drop = TRUE]
Y.test <- fake.data$y[-(1:n), , drop = TRUE]

##############################################################################
# supervised log-ratios (a balance regression method)
##############################################################################

# slr, using hierarchical clustering
slr.hclust = fitSLR(
  y = Y, X = X, intercept = intercept, 
  rho.type = rho.type, linkage = linkage, standardize = scaling)
slr.btree.hclust = slr.hclust$btree
plot(slr.btree.hclust)
slr.SBP.hclust = sbp.fromHclust(slr.btree.hclust)
# rownames(slr.SBP.hclust) = colnames(slr.SBP.hclust) = NULL

# slr, using hierarchical spectral clustering
slrMat = getSupervisedDistanceMatrix(y = Y, X = X, rho.type = rho.type)
slr.hsclust = HSClust(
  W = getSimilarityMatrix(unnormalized_similarity_matrix = slrMat), 
  levelMax = p - 1)
slr.SBP.hsclust = sbp.fromHSClust(levels_matrix = slr.hsclust$allLevels, row_names = names(beta))

#inspect
slr.SBP.hclust
dim(slr.SBP.hclust)
slr.SBP.hsclust
dim(slr.SBP.hsclust)

library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidyverse)
require(tidygraph)

# edges.hclust = getEdgesFromSBP(slr.SBP.hclust)
plotSBP(sbp = slr.SBP.hclust)
plotSBP(sbp = slr.SBP.hsclust)


