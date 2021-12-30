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

# diagonal Sigma - different sigma_eps values ##################################
diagSigma_sigmaeps = readRDS(paste0(
  output_dir, 
  "/diagSigma_threshold_maxsigma",
  ".rds"
))
sigma_eps_seq = seq(0, 4.5, by = 0.1)

# preview
length(diagSigma_sigmaeps)
idx = 10

# save plots
for(i in 1:length(diagSigma_sigmaeps)){
  # plt1.tmp = diagSigma_sigmaeps[[i]]$true_tree
  # ggsave(
  #   filename = paste0(
  #     "20211220_",
  #     "diagSigma_threshold_maxsigma", 
  #     "_noise", sigma_eps_seq[i],
  #     "_SigmaWtree.pdf"),
  #   plot = plt1.tmp,
  #   width = 8, height = 5, units = c("in")
  # )
  plt2.tmp = diagSigma_sigmaeps[[i]]$slr_tree
  # ggsave(
  #   filename = paste0(
  #     "20211220_",
  #     "diagSigma_threshold_maxsigma", 
  #     "_noise", sigma_eps_seq[i],
  #     "_slrtree.pdf"),
  #   plot = plt2.tmp,
  #   width = 8, height = 5, units = c("in")
  # )
}


# exponential decay Sigma - different rho ######################################
expdecaySigma_rho = readRDS(paste0(
  output_dir, 
  "/expdecaySigma_threshold_maxrho",
  ".rds"
))
sigma_eps = 0 ###################
rho_seq = seq(0, 1, by = 0.1)
rho_seq = rho_seq[-length(rho_seq)]

# preview
length(expdecaySigma_rho)
idx = 10
expdecaySigma_rho[[idx]]$true_tree
expdecaySigma_rho[[idx]]$slr_tree

# save plots
for(i in 1:length(expdecaySigma_rho)){
  plt1.tmp = expdecaySigma_rho[[i]]$true_tree
  # ggsave(
  #   filename = paste0(
  #     "20211220_",
  #     "expdecaySigma_threshold_maxrho", 
  #     "_noise", sigma_eps,
  #     "_rho", rho_seq[i],
  #     "_SigmaWtree.pdf"),
  #   plot = plt1.tmp,
  #   width = 8, height = 5, units = c("in")
  # )
  plt2.tmp = expdecaySigma_rho[[i]]$slr_tree
  # ggsave(
  #   filename = paste0(
  #     "20211220_",
  #     "expdecaySigma_threshold_maxrho", 
  #     "_noise", sigma_eps,
  #     "_rho", rho_seq[i],
  #     "_slrtree.pdf"),
  #   plot = plt2.tmp,
  #   width = 8, height = 5, units = c("in")
  # )
}

# 10 block Sigma - different sigma_eps values ##################################
block10Sigma_sigmaeps = readRDS(paste0(
  output_dir, 
  "/10blockSigma_threshold_maxsigma",
  ".rds"
))
sigma_eps_seq = seq(0, 4.5, by = 0.1)

# preview
length(block10Sigma_sigmaeps)
idx = 10

# save plots
for(i in 1:length(block10Sigma_sigmaeps)){
  plt1.tmp = block10Sigma_sigmaeps[[i]]$true_tree
  ggsave(
    filename = paste0(
      "20211220_",
      "10block_threshold_maxsigma",
      "_noise", sigma_eps_seq[i],
      "_SigmaWtree.pdf"),
    plot = plt1.tmp,
    width = 8, height = 5, units = c("in")
  )
  plt2.tmp = block10Sigma_sigmaeps[[i]]$slr_tree
  ggsave(
    filename = paste0(
      "20211220_",
      "10block_threshold_maxsigma", 
      "_noise", sigma_eps_seq[i],
      "_slrtree.pdf"),
    plot = plt2.tmp,
    width = 8, height = 5, units = c("in")
  )
}
