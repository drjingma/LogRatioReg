rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 2/12/2022

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
slrmax = 10
sigma_eps1 = 0.1
sigma_eps2 = 0.1
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.5 # 1, 0.5, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0
rho_alrXj = 0.2

file.end0 = paste0(
  "_", sigma.settings,
  "_", paste0(
    paste(which(SBP.true == 1), collapse = ""), "v", 
    paste(which(SBP.true == -1), collapse = "")),
  "_dim", n, "x", p, 
  "_slrmax", slrmax,
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
slr0_sims_list = list()
slr0approx_sims_list = list()
slrcv0_sims_list = list()
slrcv0approx_sims_list = list()
hslrcv0_sims_list = list()
hslrcv0approx_sims_list = list()
hslr1sccv0_sims_list = list()
hslr1sccv0approx_sims_list = list()
for(i in 1:numSims){
  print(i)
  
  # classo
  cl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cl.sim.tmp) = NULL
  classo_sims_list[[i]] = data.table(cl.sim.tmp)
  
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
  
  # cv.slr - approx
  slrcv0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrcv_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrcv0a.sim.tmp) = NULL
  slrcv0approx_sims_list[[i]] = data.table(slrcv0a.sim.tmp)
  
  # cv.slr - no approx
  slrcv0.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrcv_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrcv0.sim.tmp) = NULL
  slrcv0_sims_list[[i]] = data.table(slrcv0.sim.tmp)
  
  # cv.hslr - approx
  hslrcv0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/hslrcv_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(hslrcv0a.sim.tmp) = NULL
  hslrcv0approx_sims_list[[i]] = data.table(hslrcv0a.sim.tmp)
  
  # cv.hslr - no approx
  hslr0cv.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/hslrcv_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(hslr0cv.sim.tmp) = NULL
  hslrcv0_sims_list[[i]] = data.table(hslr0cv.sim.tmp)
  
  # cv.hslr1sc - approx
  hslr1sc0a.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/hslr1sccv_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(hslr1sc0a.sim.tmp) = NULL
  hslr1sccv0approx_sims_list[[i]] = data.table(hslr1sc0a.sim.tmp)
  
  # cv.hslr1sc - no approx
  hslr1sc0.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/hslr1sccv_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(hslr1sc0.sim.tmp) = NULL
  hslr1sccv0_sims_list[[i]] = data.table(hslr1sc0.sim.tmp)
  
}
classo_sims = as.data.frame(rbindlist(classo_sims_list))
slr0_sims = as.data.frame(rbindlist(slr0_sims_list))
slr0approx_sims = as.data.frame(rbindlist(slr0approx_sims_list))
slrcv0_sims = as.data.frame(rbindlist(slrcv0_sims_list))
slrcv0approx_sims = as.data.frame(rbindlist(slrcv0approx_sims_list))
hslrcv0_sims = as.data.frame(rbindlist(hslrcv0_sims_list))
hslrcv0approx_sims = as.data.frame(rbindlist(hslrcv0approx_sims_list))
hslr1sccv0_sims = as.data.frame(rbindlist(hslr1sccv0_sims_list))
hslr1sccv0approx_sims = as.data.frame(rbindlist(hslr1sccv0approx_sims_list))

# summary stats
slr0_summaries = data.frame(
  "mean" = apply(slr0_sims, 2, mean), 
  "sd" = apply(slr0_sims, 2, sd), 
  "se" =  apply(slr0_sims, 2, sd) / sqrt(numSims)
)
# ... and so on...

# metrics boxplots
classo_sims.gg = reshape2::melt(classo_sims)
classo_sims.gg$Method = "classo"
slr0_sims.gg = reshape2::melt(slr0_sims)
slr0_sims.gg$Method = "slr"
slr0approx_sims.gg = reshape2::melt(slr0approx_sims)
slr0approx_sims.gg$Method = "slr-a"
slrcv0_sims.gg = reshape2::melt(slrcv0_sims)
slrcv0_sims.gg$Method = "cv-slr"
slrcv0approx_sims.gg = reshape2::melt(slrcv0approx_sims)
slrcv0approx_sims.gg$Method = "cv-slr-a"
hslrcv0_sims.gg = reshape2::melt(hslrcv0_sims)
hslrcv0_sims.gg$Method = "cv-hslr"
hslrcv0approx_sims.gg = reshape2::melt(hslrcv0approx_sims)
hslrcv0approx_sims.gg$Method = "cv-hslr-a"
hslr1sccv0_sims.gg = reshape2::melt(hslr1sccv0_sims)
hslr1sccv0_sims.gg$Method = "cv-hslr1sc"
hslr1sccv0approx_sims.gg = reshape2::melt(hslr1sccv0approx_sims)
hslr1sccv0approx_sims.gg$Method = "cv-hslr1sc-a"

data.gg = rbind(
  classo_sims.gg,
  slr0_sims.gg,
  slr0approx_sims.gg, 
  slrcv0_sims.gg, 
  slrcv0approx_sims.gg, 
  hslrcv0_sims.gg, 
  hslrcv0approx_sims.gg, 
  hslr1sccv0_sims.gg, 
  hslr1sccv0approx_sims.gg)

data.gg_main = data.gg %>% 
  dplyr::filter(
    variable %in% c(
      "PEtr", "PEte", 
      "EA1", "EA2", "EAInfty",
      "FP", "FN", "TPR", "precision", 
      "Fscore", "logratios", "time"
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
