rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 4/27/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_slrs"

source("Kristyn/Functions/util.R")

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
b1 = 0.5 # 1, 0.5, 0.25
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
slr_sims_list = list()
slr_am_sims_list = list()
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
  
  # slr
  slr_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_sim_tmp) = NULL
  slr_sims_list[[i]] = data.table::data.table(slr_sim_tmp)
  
  # slr
  slr_am_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_amini_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_am_sim_tmp) = NULL
  slr_am_sims_list[[i]] = data.table::data.table(slr_am_sim_tmp)
}


# metrics boxplots
#
classo_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(classo_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo")
###
slr_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr")
#
slr_am_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr_am_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-am")

data.gg = rbind(
  classo_sims.gg,
  slr_sims.gg,
  slr_am_sims.gg
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
    "20220502",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 8, height = 6, units = c("in")
)
