rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 5/24/2022

################################################################################
# libraries and settings

output_dir = "slr_analyses/Experiments/outputs/metrics_missing"

source("slr_analyses/Functions/util.R")

library(tidyverse)
library(reshape2)

numSims = 100

sigma.settings = "latentVarModel_missing"
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_eps1 = 0.01
sigma_eps2 = 0.1
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
# (b1 = 0.5, theta.value = 0.5, a0 = 0, prop.missing = 0.75, ulimit = 0.5)
# (b1 = 1, theta.value = 0.3, a0 = 0, prop.missing = 0.70, 0.75, ulimit = 0.5)
b0 = 0 # 0
b1 = 0.5 # 0.5
theta.value = 0.5 # weight on a1 -- 1, 0.75, 0.5
a0 = 0 # 0
prop.missing = 0.5 # 0.75, 0.80

file.end0 = paste0(
  "_", sigma.settings,
  "_", prop.missing,
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
slr_spec_sims_list = list()
slr_hier_sims_list = list()
semislr_spec_sims_list = list()
semislr_spec_f_sims_list = list()
semislr_hier_sims_list = list()
semislr_hier_f_sims_list = list()
semislr_spec_nocv_sims_list = list()
semislr_hier_nocv_sims_list = list()
selbal_sims_list = list()
codacore_sims_list = list()
lrlasso_sims_list = list()
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
  
  # slr - spectral
  slr_spec_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_spectral_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_spec_sim_tmp) = NULL
  slr_spec_sims_list[[i]] = data.table::data.table(slr_spec_sim_tmp)
  
  # slr - hierarchical
  slr_hier_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hierarchical_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hier_sim_tmp) = NULL
  slr_hier_sims_list[[i]] = data.table::data.table(slr_hier_sim_tmp)
  
  # semislr - spectral
  semislr_spec_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/semislr_spectral_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(semislr_spec_sim_tmp) = NULL
  semislr_spec_sims_list[[i]] = data.table::data.table(semislr_spec_sim_tmp)
  
  # semislr - spectral
  semislr_spec_f_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/semislr_spectral_foldx2_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(semislr_spec_f_sim_tmp) = NULL
  semislr_spec_f_sims_list[[i]] = data.table::data.table(semislr_spec_f_sim_tmp)
  
  # semislr - hierarchical
  semislr_hier_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/semislr_hierarchical_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(semislr_hier_sim_tmp) = NULL
  semislr_hier_sims_list[[i]] = data.table::data.table(semislr_hier_sim_tmp)
  
  # semislr - hierarchical
  semislr_hier_f_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/semislr_hierarchical_foldx2_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(semislr_hier_f_sim_tmp) = NULL
  semislr_hier_f_sims_list[[i]] = data.table::data.table(semislr_hier_f_sim_tmp)
  
  # slr - spectral - no CV
  slr_spec_nocv_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/semislr_spectral_nocvx2_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_spec_nocv_sim_tmp) = NULL
  semislr_spec_nocv_sims_list[[i]] = data.table::data.table(slr_spec_nocv_sim_tmp)
  
  # slr - hierarchical - no CV
  slr_hier_nocv_sim_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/semislr_hierarchical_nocvx2_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr_hier_nocv_sim_tmp) = NULL
  semislr_hier_nocv_sims_list[[i]] = data.table::data.table(slr_hier_nocv_sim_tmp)
  
  ###
  
  # # selbal
  # slbl_sim_tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/selbal_metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(slbl_sim_tmp) = NULL
  # selbal_sims_list[[i]] = data.table::data.table(slbl_sim_tmp)
  
  # # codacore
  # cdcr_sim_tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/codacore_metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(cdcr_sim_tmp) = NULL
  # codacore_sims_list[[i]] = data.table::data.table(cdcr_sim_tmp)
  
  # # log-ratio lasso
  # lrl_sim_tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/lrlasso_metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(lrl_sim_tmp) = NULL
  # lrlasso_sims_list[[i]] = data.table::data.table(lrl_sim_tmp)
}

# metrics boxplots
classo_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(classo_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo")
###
slr_spec_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr_spec_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-spec")
slr_hier_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr_hier_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-hier")
#
semislr_spec_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(semislr_spec_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "semislr-spec")
semislr_spec_f_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(semislr_spec_f_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "semislr-spec-f")
semislr_hier_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(semislr_hier_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "semislr-hier")
semislr_hier_f_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(semislr_hier_f_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "semislr-hier-f")
semislr_spec_nocv_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(semislr_spec_nocv_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "semislr-spec-nocv")
semislr_hier_nocv_sims.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(semislr_hier_nocv_sims_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "semislr-hier-nocv")
###
# selbal_sims.gg = 
#   pivot_longer(as.data.frame(data.table::rbindlist(selbal_sims_list)), 
#                cols = everything(),
#                names_to = "Metric") %>%
#   mutate("Method" = "selbal")
# codacore_sims.gg = 
#   pivot_longer(as.data.frame(data.table::rbindlist(codacore_sims_list)), 
#                cols = everything(),
#                names_to = "Metric") %>%
#   mutate("Method" = "codacore")
# lrlasso_sims.gg =
#   pivot_longer(as.data.frame(data.table::rbindlist(lrlasso_sims_list)),
#                cols = everything(),
#                names_to = "Metric") %>%
#   mutate("Method" = "lrlasso")

###
data.gg = rbind(
  classo_sims.gg,
  slr_spec_sims.gg,
  slr_hier_sims.gg,
  semislr_spec_sims.gg,
  semislr_hier_sims.gg,
  semislr_spec_f_sims.gg,
  semislr_hier_f_sims.gg ,
  semislr_spec_nocv_sims.gg,
  semislr_hier_nocv_sims.gg #,
  # selbal_sims.gg, 
  # codacore_sims.gg,
  # lrlasso_sims.gg
) %>%
  dplyr::filter(
    Metric %in% c(
      # "PEtr", 
      "PEte",
      # "EA1", 
      "EA2", 
      # "EAInfty",
      "randindex", "adjrandindex",
      "TPR", "FPR", "Fscore",
      "time"
    )
  ) %>% 
  # mutate(
  #   value = ifelse(
  #     Metric %in% c("time", "EA1", "EA2", "EAInfty"), log(value), value)
  # ) %>%
  mutate(
    Metric = factor(
      Metric, 
      levels = c(
        "PEtr", "PEte",
        "EA1", "EA2", "EAInfty",
        "time",
        "randindex", "adjrandindex",
        "TPR", "FPR", "Fscore"
      ), 
      labels = c(
        "PEtr", "MSE",
        "EA1", "EA2", "EAInfty",
        "Timing",
        "RandIdx", "AdjRandIdx",
        "TPR", "FPR", "F1"
      ))
  ) %>% 
  mutate(
    Method = factor(
      Method, 
      levels = c(
        "selbal", "classo", "codacore", "lrlasso", 
        "slr-spec", "slr-hier", 
        "semislr-spec", "semislr-hier", 
        "semislr-spec-f", "semislr-hier-f", 
        "semislr-spec-nocv", "semislr-hier-nocv"
      )
    )
  )
data.gg_main = data.gg
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
    "20220726",
    file.end0,
    "_", "metrics", ".png"),
  plot = plt_main,
  width = 6, height = 6, units = c("in")
)
