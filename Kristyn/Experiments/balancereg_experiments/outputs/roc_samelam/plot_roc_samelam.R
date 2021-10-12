rm(list=ls())
# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Date: 10/11/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_experiments/outputs/roc_samelam/roccurves_samelam"

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)
library(dplyr)

numSims = 100
rng.seed = 123

# Settings to toggle with
sigma.settings = "lin14Sigma" # 2blockSigma, 4blockSigma, 10blockSigma, lin14Sigma
rho.type = "square" # 1 = "absolute value", 2 = "square"
theta.settings = "multsparse" # "dense", "sparse", "both", "multsparse"
# if "2blockSigma" then "dense"
# if "4blockSigma", then "2blocks"
# if "10blockSigma", then "pairperblock" or "1blockpair4halves"
# if "lin14Sigma" then "dense" or "multsparse"
mu.settings = "" # matchbeta
linkage = "average"
tol = 1e-4
nlam = 200
intercept = TRUE
K = 10
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
cor_ij = 0.2 # 0.2, 0.5
scaling = TRUE
sigma_eps = 0.1  # 0.1, 0.5

if(sigma.settings == "lin14Sigma"){
  if(mu.settings == "matchbeta"){
    file.end = paste0( # for old simulations
      "_dim", n, "x", p,
      "_", sigma.settings,
      "_", theta.settings,
      "_", mu.settings,
      "_noise", sigma_eps,
      "_rho", rho,
      "_int", intercept,
      "_scale", scaling,
      "_K", K,
      "_seed", rng.seed,
      ".rds")
  } else{
    file.end = paste0( # for old simulations
      "_dim", n, "x", p,
      "_", sigma.settings,
      "_", theta.settings,
      "_noise", sigma_eps,
      "_rho", rho,
      "_int", intercept,
      "_scale", scaling,
      "_K", K,
      "_seed", rng.seed,
      ".rds")
  }
} else{ # for block-diagonal Sigma, either "2blockSigma" or "4blockSigma"
  file.end = paste0(
    "_dim", n, "x", p, 
    "_", sigma.settings,
    "_", theta.settings, 
    "_noise", sigma_eps,
    "_cor", cor_ij, 
    "_int", intercept,
    "_scale", scaling,
    "_K", K,
    "_seed", rng.seed,
    ".rds")
}

has.selbal = FALSE
has.coat = FALSE
has.oracle = TRUE
has.propr = TRUE

################################################################################
# average roc curves over lambda

cl.rocs = data.frame()
slr.rocs = data.frame()
or.rocs = data.frame()
pr.rocs = data.frame()
rocs = data.frame()
for(i in 1:numSims){
  cl.roc.tmp = readRDS(paste0(output_dir, "/classo_roc", i, file.end))
  cl.roc.tmp$sim = i
  cl.rocs = rbind(cl.rocs, cl.roc.tmp)
  slr.roc.tmp = readRDS(paste0(output_dir, "/slr_roc", i, file.end))
  slr.roc.tmp$sim = i
  slr.rocs = rbind(slr.rocs, slr.roc.tmp)
  or.roc.tmp = readRDS(paste0(output_dir, "/oracle_roc", i, file.end))
  or.roc.tmp$sim = i
  or.rocs = rbind(or.rocs, or.roc.tmp)
  pr.roc.tmp = readRDS(paste0(output_dir, "/propr_roc", i, file.end))
  pr.roc.tmp$sim = i
  pr.rocs = rbind(pr.rocs, pr.roc.tmp)
  
  cl.roc.tmp$Method = "classo"
  slr.roc.tmp$Method = "slr"
  or.roc.tmp$Method = "oracle"
  pr.roc.tmp$Method = "propr"
  rocs = rbind(rocs, cl.roc.tmp, slr.roc.tmp, or.roc.tmp, pr.roc.tmp)
}

rocs = rocs %>% 
  group_by(Method, lambda) %>% 
  summarize(
    S_hat = mean(S_hat),
    tpr = mean(tpr),
    TP = mean(TP)
  ) %>% 
  ungroup()

################################################################################
# plot roc curves

# plot
rocs.gg = rocs
levels.gg = c("classo", "slr", "oracle", "propr")
rocs.gg$Method = factor(rocs.gg$Method, levels = levels.gg)
tp_roc = ggplot(
  rocs.gg, aes(x = S_hat, y = TP, color = Method)) + 
  geom_line(alpha = 0.5, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
tpr_roc = ggplot(
  rocs.gg, aes(x = S_hat, y = tpr, color = Method)) + 
  geom_line(alpha = 0.5, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
ggarrange(tp_roc, tpr_roc)
if(sigma.settings == "lin14Sigma" & mu.settings == "matchbeta"){
  ggsave(
    filename = paste0(
      "20211011_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_", mu.settings, "_rocs_samelam.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
} else{
  ggsave(
    filename = paste0(
      "20211011_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_rocs.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
}


