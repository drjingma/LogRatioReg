rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 4/13/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_slrs"

source("Kristyn/Functions/util.R")

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100

sigma.settings = "latentVarModel_corX"
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_eps1 = 0.1
sigma_eps2 = 0.1
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.25 # 1, 0.5, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0
rho_alrXj = 0.2

file.end0 = paste0(
  "_", sigma.settings,
  "_", paste0(
    paste(which(SBP.true == 1), collapse = ""), "v", 
    paste(which(SBP.true == -1), collapse = "")),
  "_dim", n, "x", p, 
  "_noisey", sigma_eps1, 
  "_noisex", sigma_eps2,
  "_b0", b0, 
  "_b1", b1, 
  "_a0", a0, 
  "_theta", theta.value, 
  "_rho", rho_alrXj)

################################################################################
# plot metrics

# import metrics
classo_sims_list = list()
slr_sims_list = list()
slr_hdr_original_sims_list = list()
slr_hdr_mean_sims_list = list()
slr_hdr_median_sims_list = list()
for(i in 1:numSims){
  print(i)
  
  # classo
  cl_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cl_sim_tmp) = NULL
  classo_sims_list[[i]] = data.table(cl_sim_tmp)
  
  # slr
  slr_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_sim_tmp) = NULL
  slr_sims_list[[i]] = data.table(slr_sim_tmp)
  
  # slr - hi-degree reg original
  slr_hdr_original_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hdr_original_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hdr_original_sim_tmp) = NULL
  slr_hdr_original_sims_list[[i]] = data.table(slr_hdr_original_sim_tmp)
  
  # slr - hi-degree reg mean
  slr_hdr_mean_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hdr_mean_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hdr_mean_sim_tmp) = NULL
  slr_hdr_mean_sims_list[[i]] = data.table(slr_hdr_mean_sim_tmp)
  
  # slr - hi-degree reg median
  slr_hdr_median_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hdr_median_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hdr_median_sim_tmp) = NULL
  slr_hdr_median_sims_list[[i]] = data.table(slr_hdr_median_sim_tmp)
  
}
# metrics boxplots
classo_sims.gg = reshape2::melt(as.data.frame(rbindlist(classo_sims_list)))
classo_sims.gg$Method = "classo"
#
slr_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_sims_list)))
slr_sims.gg$Method = "slr"
#
slr_hdr_original_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_hdr_original_sims_list)))
slr_hdr_original_sims.gg$Method = "slr-hdr-orig"
#
slr_hdr_mean_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_hdr_mean_sims_list)))
slr_hdr_mean_sims.gg$Method = "slr-hdr-mean"
#
slr_hdr_median_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_hdr_median_sims_list)))
slr_hdr_median_sims.gg$Method = "slr-hdr-med"

data.gg = rbind(
  classo_sims.gg,
  slr_sims.gg,
  slr_hdr_original_sims.gg, 
  slr_hdr_mean_sims.gg, 
  slr_hdr_median_sims.gg)

data.gg_main = data.gg %>% 
  dplyr::filter(
    variable %in% c(
      "PEtr", "PEte", 
      "EA1", "EA2", "EAInfty",
      "FP", "FN", "TPR", "precision", 
      "Fscore", "time"
    )
  )
plt_main = ggplot(
  data.gg_main, 
  aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(variable), scales = "free_y") +
  geom_boxplot() +
  stat_summary(
    fun = mean, geom = "point", shape = 4, size = 1.5,
    color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank())
plt_main
ggsave(
  filename = paste0(
    "20220413_slr_hdr",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 8, height = 6, units = c("in")
)
