rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 4/27/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_slrs"

source("Kristyn/Functions/util.R")

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)

numSims = 100

sigma.settings = "latentVarModel"
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
b1 = 1 # 1, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0

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
  "_theta", theta.value)

################################################################################
# plot metrics

# import metrics
classo_sims_list = list()
slrk_sims_list = list()
slrk_am_sims_list = list()
slrk_hdr_sims_list = list()
slrk_mg_sims_list = list()
slrk_mg_am_sims_list = list()
slrk_mg_hdr_sims_list = list()
slrc_sims_list = list()
slrc_am_sims_list = list()
slrc_hdr_sims_list = list()
slrc_mg_sims_list = list()
slrc_mg_am_sims_list = list()
slrc_mg_hdr_sims_list = list()
slr1sc_sims_list = list()
slr1sc_mg_sims_list = list()
for(i in 1:numSims){
  print(i)
  
  # compositional lasso
  cl_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cl_sim_tmp) = NULL
  classo_sims_list[[i]] = data.table::data.table(cl_sim_tmp)
  
  ###
  
  # slr - 1minusGamma - kmeans
  slrk_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_kmeans_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrk_sim_tmp) = NULL
  slrk_sims_list[[i]] = data.table::data.table(slrk_sim_tmp)
  
  # slr - 1minusGamma - kmeans - amini
  slrk_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_kmeans_amini_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrk_sim_tmp) = NULL
  slrk_am_sims_list[[i]] = data.table::data.table(slrk_sim_tmp)
  
  # slr - 1minusGamma - kmeans - hdr
  slrk_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_kmeans_hdr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrk_sim_tmp) = NULL
  slrk_hdr_sims_list[[i]] = data.table::data.table(slrk_sim_tmp)
  
  #
  
  # slr - maxGammaMinusGamma - kmeans
  slrk_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_maxGamma_kmeans_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrk_sim_tmp) = NULL
  slrk_mg_sims_list[[i]] = data.table::data.table(slrk_sim_tmp)
  
  # slr - maxGammaMinusGamma - kmeans - amini
  slrk_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_maxGamma_kmeans_amini_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrk_sim_tmp) = NULL
  slrk_mg_am_sims_list[[i]] = data.table::data.table(slrk_sim_tmp)
  
  # slr - maxGammaMinusGamma - kmeans - hdr
  slrk_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_maxGamma_kmeans_hdr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrk_sim_tmp) = NULL
  slrk_mg_hdr_sims_list[[i]] = data.table::data.table(slrk_sim_tmp)
  
  ###
  
  # slr - 1minusGamma - cut
  slrc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_cut_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrc_sim_tmp) = NULL
  slrc_sims_list[[i]] = data.table::data.table(slrc_sim_tmp)
  
  # slr - 1minusGamma - cut - amini
  slrc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_cut_amini_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrc_sim_tmp) = NULL
  slrc_am_sims_list[[i]] = data.table::data.table(slrc_sim_tmp)
  
  # slr - 1minusGamma - cut - hdr
  slrc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_cut_hdr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrc_sim_tmp) = NULL
  slrc_hdr_sims_list[[i]] = data.table::data.table(slrc_sim_tmp)
  
  #
  
  # slr - maxGammaMinusGamma - cut
  slrc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_maxGamma_cut_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrc_sim_tmp) = NULL
  slrc_mg_sims_list[[i]] = data.table::data.table(slrc_sim_tmp)
  
  # slr - maxGammaMinusGamma - cut - amini
  slrc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_maxGamma_cut_amini_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrc_sim_tmp) = NULL
  slrc_mg_am_sims_list[[i]] = data.table::data.table(slrc_sim_tmp)
  
  # slr - maxGammaMinusGamma - cut - hdr
  slrc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_maxGamma_cut_hdr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrc_sim_tmp) = NULL
  slrc_mg_hdr_sims_list[[i]] = data.table::data.table(slrc_sim_tmp)
  
  ###
  
  # slr1sc - 1minusGamma - cut
  slr1sc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_1sc_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr1sc_sim_tmp) = NULL
  slr1sc_sims_list[[i]] = data.table::data.table(slr1sc_sim_tmp)
  
  # slr1sc - maxGammaMinusGamma - cut
  slr1sc_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_1sc_maxGamma_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr1sc_sim_tmp) = NULL
  slr1sc_mg_sims_list[[i]] = data.table::data.table(slr1sc_sim_tmp)
}


# metrics boxplots
#
classo_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(classo_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo")
###
slrk_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrk_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrk")
#
slrk_am_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrk_am_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrk-am")
#
slrk_hdr_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrk_hdr_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrk-hdr")
###
slrk_mg_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrk_mg_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrk-mg")
#
slrk_mg_am_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrk_mg_am_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrk-mg-am")
#
slrk_mg_hdr_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrk_mg_hdr_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrk-mg-hdr")
###
###
slrc_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrc_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrc")
slrc_am_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrc_am_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrc-am")
#
slrc_hdr_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrc_hdr_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrc-hdr")
###
slrc_mg_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrc_mg_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrc-mg")
#
slrc_mg_am_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrc_mg_am_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrc-mg-am")
#
slrc_mg_hdr_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrc_mg_hdr_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slrc-mg-hdr")
###
###
slr1sc_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr1sc_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr1sc")
#
slr1sc_mg_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr1sc_mg_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr1sc-mg")
###

data.gg = rbind(
  classo_sims.gg,
  slrk_sims.gg,
  slrk_am_sims.gg,
  slrk_hdr_sims.gg,
  slrk_mg_sims.gg,
  slrk_mg_am_sims.gg,
  slrk_mg_hdr_sims.gg,
  slrc_sims.gg, 
  slrc_am_sims.gg,
  slrc_hdr_sims.gg,
  slrc_mg_sims.gg,
  slrc_mg_am_sims.gg,
  slrc_mg_hdr_sims.gg,
  slr1sc_sims.gg,
  slr1sc_mg_sims.gg
) %>%
  mutate(
    Metric = factor(
      Metric, levels = c(
        "PEtr", "PEte", 
        "EA1", "EA2", "EAInfty", 
        "EA1Active", "EA2Active", "EAInftyActive",
        "EA1Inactive", "EA2Inactive", "EAInftyInactive",
        "FP", "FN", "TPR", "precision", "Fscore",
        "FP+", "FN+", "TPR+", 
        "FP-", "FN-", "TPR-", 
        "betasparsity", "logratios", "time"
        ))
  )

data.gg_main = data.gg %>% 
  dplyr::filter(
    Metric %in% c(
      "PEtr", "PEte", 
      "EA1", "EA2", "EAInfty",
      "FP", "FN", "TPR", "precision", 
      "Fscore", "time"
    )
  # ) %>%
  # mutate(
  #   Method = factor(Method, levels = c(
  #     "slr", "slr-s1",
  #     "slr-am", "slr-am-s1",
  #     "slr-ap", "slr-ap-s1",
  #     "slr-am-ap", "slr-am-ap-s1"
  #   ), labels = c(
  #     "slr", "slr-s1",
  #     "slr-am", "slr-am-s1",
  #     "slr-ap", "slr-ap-s1",
  #     "slr-am-ap", "slr-am-ap-s1"
  #   ))
  )
plt_main = ggplot(
  data.gg_main, 
  aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(Metric), scales = "free_y") +
  geom_boxplot() +
  stat_summary(
    fun = mean, geom = "point", shape = 4, size = 1.5,
    color = "red") +
  # scale_x_discrete(limits = c(
  #   "slr", "slr-s1",
  #   "slr-am", "slr-am-s1",
  #   "slr-ap", "slr-ap-s1",
  #   "slr-am-ap", "slr-am-ap-s1"
  # )) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank())
plt_main
ggsave(
  filename = paste0(
    "20220427",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 8, height = 6, units = c("in")
)
