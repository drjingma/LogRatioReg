# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 12/19/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_hsclust_experiments/outputs"

library(mvtnorm)

library(Matrix)
library(glmnet)

source("RCode/func_libs.R")
source("Kristyn/Functions/supervisedlogratios.R")
source("Kristyn/Functions/HSClust.R")

# helper functions
source("Kristyn/Functions/metrics.R")
source("Kristyn/Functions/simulatedata.R")

# for plots
library(ggraph) # make dendrogram
library(igraph) # transform dataframe to graph object: graph_from_data_frame()
library(tidygraph)

# diagonal Sigma - different sigma_eps values
diagSigma_sigmaeps = readRDS(paste0(
  output_dir, 
  "/diagSigma_threshold_maxsigma",
  ".rds"
))
diagSigma_sigmaeps

expdecaySigma_rho = readRDS(paste0(
  output_dir, 
  "/expdecaySigma_threshold_maxrho",
  ".rds"
))

