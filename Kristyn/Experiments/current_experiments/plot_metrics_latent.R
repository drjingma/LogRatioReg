rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 2/12/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/current_experiments/outputs"

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
cl_sims_list = list()
slrnew_sims_list = list()
slrnew2_sims_list = list()
slbl_sims_list = list()
cdcr_sims_list = list()
for(i in 1:numSims){
  print(i)
  
  # classo
  cl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/classo_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cl.sim.tmp) = NULL
  cl_sims_list[[i]] = data.table(cl.sim.tmp)
  
  # slr
  slrnew.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrnew.sim.tmp) = NULL
  slrnew_sims_list[[i]] = data.table(slrnew.sim.tmp)
  
  # slr with rank 1 approximation
  slrnew2.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_approx_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrnew2.sim.tmp) = NULL
  slrnew2_sims_list[[i]] = data.table(slrnew2.sim.tmp)
  
  # selbal
  slbl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/selbal_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slbl.sim.tmp) = NULL
  slbl_sims_list[[i]] = data.table(slbl.sim.tmp)
  
  # codacore
  cdcr.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/codacore_metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cdcr.sim.tmp) = NULL
  cdcr_sims_list[[i]] = data.table(cdcr.sim.tmp)
}

cl_sims = as.data.frame(rbindlist(cl_sims_list))
slrnew_sims = as.data.frame(rbindlist(slrnew_sims_list))
slrnew2_sims = as.data.frame(rbindlist(slrnew2_sims_list))
slbl_sims = as.data.frame(rbindlist(slbl_sims_list))
cdcr_sims = as.data.frame(rbindlist(cdcr_sims_list))

# summary stats
cl_summaries = data.frame(
  "mean" = apply(cl_sims, 2, mean), 
  "sd" = apply(cl_sims, 2, sd), 
  "se" =  apply(cl_sims, 2, sd) / sqrt(numSims)
)
slrnew_summaries = data.frame(
  "mean" = apply(slrnew_sims, 2, mean), 
  "sd" = apply(slrnew_sims, 2, sd), 
  "se" =  apply(slrnew_sims, 2, sd) / sqrt(numSims)
)
slrnew2_summaries = data.frame(
  "mean" = apply(slrnew2_sims, 2, mean), 
  "sd" = apply(slrnew2_sims, 2, sd), 
  "se" =  apply(slrnew2_sims, 2, sd) / sqrt(numSims)
)
slbl_summaries = data.frame(
  "mean" = apply(slbl_sims, 2, mean),
  "sd" = apply(slbl_sims, 2, sd),
  "se" =  apply(slbl_sims, 2, sd) / sqrt(numSims)
)
cdcr_summaries = data.frame(
  "mean" = apply(cdcr_sims, 2, mean),
  "sd" = apply(cdcr_sims, 2, sd),
  "se" =  apply(cdcr_sims, 2, sd) / sqrt(numSims)
)

# metrics boxplots
cl.sims.gg = reshape2::melt(cl_sims)
cl.sims.gg$Method = "classo"
slrnew.sims.gg = reshape2::melt(slrnew_sims)
slrnew.sims.gg$Method = "slr"
slrnew2.sims.gg = reshape2::melt(slrnew2_sims)
slrnew2.sims.gg$Method = "slr-approx"
slbl.sims.gg = reshape2::melt(slbl_sims)
slbl.sims.gg$Method = "selbal"
cdcr.sims.gg = reshape2::melt(cdcr_sims)
cdcr.sims.gg$Method = "codacore"

data.gg = rbind(
  cl.sims.gg,
  slrnew.sims.gg, 
  slrnew2.sims.gg, 
  slbl.sims.gg, 
  cdcr.sims.gg)

data.gg_main = data.gg %>% 
  dplyr::filter(
    variable %in% c(
      "PEtr", "PEte", 
      "EA1", "EA2", "EAInfty",
      "FP", "FN", "TPR", "precision", 
      "Fscore", "time"
    )
  ) %>%
  dplyr::filter(
    !(Method %in% c("selbal", "codacore") & variable == "time")
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
    "20220302",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
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
    "20220302",
    file.end0, 
    "_", "metrics_posneg", ".pdf"),
  plot = last_plot(),
  width = 8, height = 5, units = c("in")
)
