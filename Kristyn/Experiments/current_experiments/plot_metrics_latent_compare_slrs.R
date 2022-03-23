rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 2/12/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs/metrics_slrs"

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100

# Settings to toggle with
sigma.settings = "latentVarModel"
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_eps1 = 0.1 # MUST BE LESS THAN OR EQUAL TO sigma_eps2
sigma_eps2 = 0.1
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0
b1 = 0.25 # 1, 0.5, 0.25
a0 = 0
theta.value = 1 # weight on a1

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
slr0_sims_list = list()
slr0approx_sims_list = list()
slrmult0_sims_list = list()
slrmult0approx_sims_list = list()
hslr0_sims_list = list()
hslr0approx_sims_list = list()
slrmultcv0_sims_list = list()
slrmultcv0approx_sims_list = list()
for(i in 1:numSims){
  print(i)
  
  # plain slr - approx
  slr0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr0a.sim.tmp) = NULL
  slr0approx_sims_list[[i]] = data.table(slr0a.sim.tmp)
  
  # plain slr - no approx
  slr0.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slr0.sim.tmp) = NULL
  slr0_sims_list[[i]] = data.table(slr0.sim.tmp)
  
  # mult slr - rsq - approx
  slrmult0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrmult_rsq_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrmult0a.sim.tmp) = NULL
  slrmult0approx_sims_list[[i]] = data.table(slrmult0a.sim.tmp)
  
  # mult slr - rsq - no approx
  slrmult0.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrmult_rsq_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrmult0.sim.tmp) = NULL
  slrmult0_sims_list[[i]] = data.table(slrmult0.sim.tmp)
  
  # hslr - approx
  hslr0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/hslr_rsq_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(hslr0a.sim.tmp) = NULL
  hslr0approx_sims_list[[i]] = data.table(hslr0a.sim.tmp)
  
  # hslr - no approx
  hslr0.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/hslr_rsq_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(hslr0.sim.tmp) = NULL
  hslr0_sims_list[[i]] = data.table(hslr0.sim.tmp)
  
  # mult slr - cv - approx
  slrmultcv0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrmult_cv_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrmultcv0a.sim.tmp) = NULL
  slrmultcv0approx_sims_list[[i]] = data.table(slrmultcv0a.sim.tmp)
  
  # mult slr - cv - no approx
  slrmultcv0.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrmult_cv_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrmultcv0.sim.tmp) = NULL
  slrmultcv0_sims_list[[i]] = data.table(slrmultcv0.sim.tmp)
  
}

slr0_sims = as.data.frame(rbindlist(slr0_sims_list))
slr0approx_sims = as.data.frame(rbindlist(slr0approx_sims_list))
slrmult0_sims = as.data.frame(rbindlist(slrmult0_sims_list))
slrmult0approx_sims = as.data.frame(rbindlist(slrmult0approx_sims_list))
hslr0_sims = as.data.frame(rbindlist(hslr0_sims_list))
hslr0approx_sims = as.data.frame(rbindlist(hslr0approx_sims_list))
slrmultcv0_sims = as.data.frame(rbindlist(slrmultcv0_sims_list))
slrmultcv0approx_sims = as.data.frame(rbindlist(slrmultcv0approx_sims_list))

# summary stats
slr0_summaries = data.frame(
  "mean" = apply(slr0_sims, 2, mean), 
  "sd" = apply(slr0_sims, 2, sd), 
  "se" =  apply(slr0_sims, 2, sd) / sqrt(numSims)
)
# ... and so on...

# metrics boxplots
slr0_sims.gg = reshape2::melt(slr0_sims)
slr0_sims.gg$Method = "slr"
slr0approx_sims.gg = reshape2::melt(slr0approx_sims)
slr0approx_sims.gg$Method = "slr-appr"
slrmult0_sims.gg = reshape2::melt(slrmult0_sims)
slrmult0_sims.gg$Method = "mslr-rsq"
slrmult0approx_sims.gg = reshape2::melt(slrmult0approx_sims)
slrmult0approx_sims.gg$Method = "mslr-rsq-appr"
hslr0_sims.gg = reshape2::melt(hslr0_sims)
hslr0_sims.gg$Method = "hslr-rsq"
hslr0approx_sims.gg = reshape2::melt(hslr0approx_sims)
hslr0approx_sims.gg$Method = "hslr-rsq-appr"
slrmultcv0_sims.gg = reshape2::melt(slrmultcv0_sims)
slrmultcv0_sims.gg$Method = "mslr-cv"
slrmultcv0approx_sims.gg = reshape2::melt(slrmultcv0approx_sims)
slrmultcv0approx_sims.gg$Method = "mslr-cv-appr"

data.gg = rbind(
  slr0_sims.gg,
  slr0approx_sims.gg, 
  slrmult0_sims.gg, 
  slrmult0approx_sims.gg, 
  hslr0_sims.gg, 
  hslr0approx_sims.gg, 
  slrmultcv0_sims.gg, 
  slrmultcv0approx_sims.gg)

data.gg_main = data.gg %>% 
  dplyr::filter(
    variable %in% c(
      "PEte", 
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
    "20220323",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 8, height = 6, units = c("in")
)

data.gg_main2 = data.gg_main %>%
  dplyr::filter(
    !(Method %in% c("mslr-cv", "mslr-cv-appr") & variable == "time")
  )
plt_main2 = ggplot(
  data.gg_main2, 
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
plt_main2
ggsave(
  filename = paste0(
    "20220323",
    file.end0,
    "_", "metrics_exclude", ".pdf"),
  plot = plt_main2,
  width = 8, height = 6, units = c("in")
)

data.gg_pos = data.gg %>% dplyr::filter(
  variable %in% c(
    "FP+", "FN+", "TPR+"
  )
)
plt_pos = ggplot(
  data.gg_pos, 
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
data.gg_neg = data.gg %>% dplyr::filter(
  variable %in% c(
    "FP-", "FN-", "TPR-"
  )
)
plt_neg = ggplot(
  data.gg_neg, 
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
ggarrange(plt_pos, plt_neg, nrow = 2)
ggsave(
  filename = paste0(
    "20220315",
    file.end0, 
    "_", "metrics_posneg", ".pdf"),
  plot = last_plot(),
  width = 8, height = 5, units = c("in")
)
