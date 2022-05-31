rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 2/12/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_binary"

source("Kristyn/Functions/util.R")

library(tidyverse)
library(reshape2)

numSims = 100 

# Settings to toggle with
sigma.settings = "latentVarModel_binary"
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
# sigma_eps1 = 0.1
sigma_eps2 = 0.1
# SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 4 # 2, 4
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0

file.end0 = paste0(
  "_", sigma.settings,
  "_", paste0(
    paste(which(SBP.true == 1), collapse = ""), "v", 
    paste(which(SBP.true == -1), collapse = "")),
  "_dim", n, "x", p, 
  "_noisex", sigma_eps2,
  "_b0", b0, 
  "_b1", b1, 
  "_a0", a0, 
  "_theta", theta.value)

################################################################################
# plot metrics

# import metrics
classo_sims_list = list()
slr_0.05_sims_list = list()
slr_0.01_sims_list = list()
slrscreen_sims_list = list()
selbal_sims_list = list()
codacore_sims_list = list()
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
  
  # slr - alpha = 0.05
  slr_alpha0.05_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_alpha0.05_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_alpha0.05_sim_tmp) = NULL
  slr_0.05_sims_list[[i]] = data.table::data.table(slr_alpha0.05_sim_tmp)
  
  # slr - alpha = 0.01
  slr_alpha0.01_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_alpha0.01_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_alpha0.01_sim_tmp) = NULL
  slr_0.01_sims_list[[i]] = data.table::data.table(slr_alpha0.01_sim_tmp)
  
  # slr - screen
  slrscreen_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrscreen_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrscreen_sim_tmp) = NULL
  slrscreen_sims_list[[i]] = data.table::data.table(slrscreen_sim_tmp)
  ###
  
  # selbal
  slbl_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/selbal_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slbl_sim_tmp) = NULL
  selbal_sims_list[[i]] = data.table::data.table(slbl_sim_tmp)
  
  # codacore
  cdcr_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/codacore_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cdcr_sim_tmp) = NULL
  codacore_sims_list[[i]] = data.table::data.table(cdcr_sim_tmp)
}

# metrics boxplots
classo_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(classo_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo")
###
slr_0.05_sims.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(slr_0.05_sims_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-0.05")
slr_0.01_sims.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(slr_0.01_sims_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-0.01")
slrscreen_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrscreen_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr")
###
selbal_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(selbal_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "selbal")
codacore_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(codacore_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "codacore")
###
data.gg = rbind(
  classo_sims.gg,
  # slr_0.05_sims.gg,
  # slr_0.01_sims.gg,
  slrscreen_sims.gg,
  selbal_sims.gg, 
  codacore_sims.gg
) %>%
  mutate(
    Metric = factor(
      Metric, levels = c(
        "AUCtr", "AUCte", 
        "EA1", "EA2", "EAInfty", 
        "TPR", "FPR", "Fscore",
        "betasparsity", "logratios", "adhoc", "time"
      ))
  )

data.gg_main = data.gg %>%
  dplyr::filter(
    Metric %in% c(
      "AUCtr", "AUCte",
      "EA1", "EA2", "EAInfty",
      "TPR", "FPR", "Fscore",
      "time"
    )
  )
plt_main = ggplot(
  data.gg_main, 
  aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(Metric), scales = "free_y") +
  geom_boxplot() +
  stat_summary(
    fun = mean, geom = "point", shape = 4, size = 1.5,
    color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank())
plt_main
ggsave(
  filename = paste0(
    "20220530",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 6, height = 5, units = c("in")
)
# data.gg %>% filter(Metric == "adhoc") %>% 
#   group_by(Metric, Method) %>%
#   summarize(percentage = mean(value, na.rm = TRUE)) %>%
#   ungroup()