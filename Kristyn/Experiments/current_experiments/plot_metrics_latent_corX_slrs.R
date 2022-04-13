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
# SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
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
slr_am_sims_list = list()
slr_ap_sims_list = list()
slr_am_ap_sims_list = list()
slr_hdr_sims_list = list()
slr_hdr_ap_sims_list = list()
slr_eig12_sims_list = list()
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
  
  # slr - amini
  slr_am_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_amini_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_am_sim_tmp) = NULL
  slr_am_sims_list[[i]] = data.table(slr_am_sim_tmp)
  
  # slr - approx
  slr_ap_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_ap_sim_tmp) = NULL
  slr_ap_sims_list[[i]] = data.table(slr_ap_sim_tmp)
  
  # slr - amini + approx
  slr_am_ap_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_amini_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_am_ap_sim_tmp) = NULL
  slr_am_ap_sims_list[[i]] = data.table(slr_am_ap_sim_tmp)
  
  # slr - high-degree reg
  slr_hdr_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hdr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hdr_sim_tmp) = NULL
  slr_hdr_sims_list[[i]] = data.table(slr_hdr_sim_tmp)
  
  # slr - high-degree reg + approx
  slr_hdr_ap_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hdr_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hdr_ap_sim_tmp) = NULL
  slr_hdr_ap_sims_list[[i]] = data.table(slr_hdr_ap_sim_tmp)
  
  # slr - 2 leading eigenvectors
  slr_eig12_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_eig12_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_eig12_sim_tmp) = NULL
  slr_eig12_sims_list[[i]] = data.table(slr_eig12_sim_tmp)
  
}
# metrics boxplots
classo_sims.gg = reshape2::melt(as.data.frame(rbindlist(classo_sims_list)))
classo_sims.gg$Method = "classo"
#
slr_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_sims_list)))
slr_sims.gg$Method = "slr"
#
slr_am_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_am_sims_list)))
slr_am_sims.gg$Method = "slr-am"
#
slr_ap_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_ap_sims_list)))
slr_ap_sims.gg$Method = "slr-ap"
#
slr_am_ap_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_am_ap_sims_list)))
slr_am_ap_sims.gg$Method = "slr-am-ap"
#
slr_hdr_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_hdr_sims_list)))
slr_hdr_sims.gg$Method = "slr-hdr"
#
slr_hdr_ap_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_hdr_ap_sims_list)))
slr_hdr_ap_sims.gg$Method = "slr-hdr-ap"
#
slr_eig12_sims.gg = reshape2::melt(as.data.frame(rbindlist(slr_eig12_sims_list)))
slr_eig12_sims.gg$Method = "slr-eig12"

data.gg = rbind(
  classo_sims.gg,
  slr_sims.gg,
  slr_am_sims.gg, 
  slr_ap_sims.gg, 
  slr_am_ap_sims.gg,
  slr_hdr_sims.gg,
  slr_hdr_ap_sims.gg,
  slr_eig12_sims.gg)

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
    "20220413",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 8, height = 6, units = c("in")
)

# data.gg_main2 = data.gg_main %>%
#   dplyr::filter(
#     !(Method %in% c("mslr-cv", "mslr-cv-appr") & variable == "time")
#   )
# plt_main2 = ggplot(
#   data.gg_main2, 
#   aes(x = Method, y = value, color = Method)) +
#   facet_wrap(vars(variable), scales = "free_y") +
#   geom_boxplot() +
#   stat_summary(
#     fun = mean, geom = "point", shape = 4, size = 1.5,
#     color = "red") +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(), 
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#     axis.title.y = element_blank())
# plt_main2
# ggsave(
#   filename = paste0(
#     "20220406",
#     file.end0,
#     "_", "metrics_exclude", ".pdf"),
#   plot = plt_main2,
#   width = 8, height = 6, units = c("in")
# )
# 
# data.gg_pos = data.gg %>% dplyr::filter(
#   variable %in% c(
#     "FP+", "FN+", "TPR+"
#   )
# )
# plt_pos = ggplot(
#   data.gg_pos, 
#   aes(x = Method, y = value, color = Method)) +
#   facet_wrap(vars(variable), scales = "free_y") +
#   geom_boxplot() +
#   stat_summary(
#     fun = mean, geom = "point", shape = 4, size = 1.5,
#     color = "red") +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(), 
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#     axis.title.y = element_blank())
# data.gg_neg = data.gg %>% dplyr::filter(
#   variable %in% c(
#     "FP-", "FN-", "TPR-"
#   )
# )
# plt_neg = ggplot(
#   data.gg_neg, 
#   aes(x = Method, y = value, color = Method)) +
#   facet_wrap(vars(variable), scales = "free_y") +
#   geom_boxplot() +
#   stat_summary(
#     fun = mean, geom = "point", shape = 4, size = 1.5,
#     color = "red") +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(), 
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#     axis.title.y = element_blank())
# ggarrange(plt_pos, plt_neg, nrow = 2)
# ggsave(
#   filename = paste0(
#     "20220315",
#     file.end0, 
#     "_", "metrics_posneg", ".pdf"),
#   plot = last_plot(),
#   width = 8, height = 5, units = c("in")
# )
