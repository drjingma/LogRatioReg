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
sigma.settings = "diagSigma"
theta.value = 1
intercept = TRUE
K = 10
n = 100
p = 30
scaling = TRUE
linkage = "average"
tol = 1e-4
nlam = 100
neta = p
#################
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
#################
rho = 0.2 #
desired_Rsquared = 0.6 #
# sigma_eps1 = get_sigma_eps(
#   sbp = SBP.true, ilr.trans.constant = ilrtrans.true$const, theta = theta.value, 
#   Rsq = desired_Rsquared, rho = rho)
sigma_eps1 = 0.1
sigma_eps2 = 0.1

file.end0 = paste0(
  "_", sigma.settings,
  "_", paste0(
    paste(which(SBP.true == 1), collapse = ""), "v", 
    paste(which(SBP.true == -1), collapse = "")),
  "_dim", n, "x", p, 
  # "_rho", rho, 
  "_noisey", sigma_eps1, 
  "_noisex", sigma_eps2)

################################################################################
# plot metrics

# import metrics
slrhsc_thresh_lasso_sims_list = list()
cl_sims_list = list()
slrnew_sims_list = list()
slrnew2_sims_list = list()
pr_sims_list = list()
# slbl_sims_list = list()
for(i in 1:numSims){
  print(i)
  
  # slr hsc thresh lasso
  slrhsc_thresh_lasso.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hsc_thresh_lasso_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhsc_thresh_lasso.sim.tmp) = NULL
  slrhsc_thresh_lasso_sims_list[[i]] = data.table(slrhsc_thresh_lasso.sim.tmp)
  
  # classo
  cl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/classo_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cl.sim.tmp) = NULL
  cl_sims_list[[i]] = data.table(cl.sim.tmp)
  
  # slrnew
  slrnew.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_new_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrnew.sim.tmp) = NULL
  slrnew_sims_list[[i]] = data.table(slrnew.sim.tmp)
  
  # slrnew - no rank 1 approximation
  slrnew2.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_new_noapprox_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrnew2.sim.tmp) = NULL
  slrnew2_sims_list[[i]] = data.table(slrnew2.sim.tmp)
  
  # propr
  pr.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/propr_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(pr.sim.tmp) = NULL
  pr_sims_list[[i]] = data.table(pr.sim.tmp)
  
  # # selbal
  # slbl.sim.tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/metrics", "/selbal_", "metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(slbl.sim.tmp) = NULL
  # slbl_sims_list[[i]] = data.table(slbl.sim.tmp)
}

slrhsc_thresh_lasso_sims = as.data.frame(rbindlist(slrhsc_thresh_lasso_sims_list))
cl_sims = as.data.frame(rbindlist(cl_sims_list))
slrnew_sims = as.data.frame(rbindlist(slrnew_sims_list))
slrnew2_sims = as.data.frame(rbindlist(slrnew2_sims_list))
pr_sims = as.data.frame(rbindlist(pr_sims_list))
# slbl_sims = as.data.frame(rbindlist(slbl_sims_list))

# summary stats
slrhsc_thresh_lasso_summaries = data.frame(
  "mean" = apply(slrhsc_thresh_lasso_sims, 2, mean), 
  "sd" = apply(slrhsc_thresh_lasso_sims, 2, sd), 
  "se" =  apply(slrhsc_thresh_lasso_sims, 2, sd) / sqrt(numSims)
)
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
pr_summaries = data.frame(
  "mean" = apply(pr_sims, 2, mean),
  "sd" = apply(pr_sims, 2, sd),
  "se" =  apply(pr_sims, 2, sd) / sqrt(numSims)
)
# slbl_summaries = data.frame(
#   "mean" = apply(slbl_sims, 2, mean), 
#   "sd" = apply(slbl_sims, 2, sd), 
#   "se" =  apply(slbl_sims, 2, sd) / sqrt(numSims)
# )

# metrics boxplots
slrhsc_thresh_lasso.sims.gg = reshape2::melt(slrhsc_thresh_lasso_sims)
slrhsc_thresh_lasso.sims.gg$Method = "slr-thresh-lasso"
cl.sims.gg = reshape2::melt(cl_sims)
cl.sims.gg$Method = "classo"
slrnew.sims.gg = reshape2::melt(slrnew_sims)
slrnew.sims.gg$Method = "slr-new-approx"
slrnew2.sims.gg = reshape2::melt(slrnew2_sims)
slrnew2.sims.gg$Method = "slr-new"
pr.sims.gg = reshape2::melt(pr_sims)
pr.sims.gg$Method = "propr"
# slbl.sims.gg = reshape2::melt(slbl_sims)
# slbl.sims.gg$Method = "selbal"

data.gg = rbind(
  slrhsc_thresh_lasso.sims.gg, 
  cl.sims.gg,
  slrnew.sims.gg, 
  slrnew2.sims.gg, 
  pr.sims.gg)
# pr.sims.gg, 
# or.sims.gg)#, 
# slbl.sims.gg)

data.gg_main = data.gg %>% dplyr::filter(
  variable %in% c(
    "PEtr", "PEte", 
    # "EA1", "EA2", "EAInfty", 
    "FP", "FN", "TPR", "precision", 
    "Fscore", "time"
  )
)
plt_main = ggplot(
  data.gg_main, 
  aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(variable), scales = "free_y") +
  geom_boxplot() +
  # stat_summary(
  #   fun = mean, fun.min = mean, fun.max = mean,
  #   geom = "errorbar", width = 0.75,
  #   linetype = "dashed") +
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
# ggsave(
#   filename = paste0(
#     "20220222",
#     file.end0, 
#     "_", "metrics", ".pdf"),
#   plot = plt_main,
#   width = 8, height = 6, units = c("in")
# )

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
  # stat_summary(
  #   fun = mean, fun.min = mean, fun.max = mean,
  #   geom = "errorbar", width = 0.75,
  #   linetype = "dashed") +
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
  # stat_summary(
  #   fun = mean, fun.min = mean, fun.max = mean,
  #   geom = "errorbar", width = 0.75,
  #   linetype = "dashed") +
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
    "20220222",
    file.end0, 
    "_", "metrics_posneg", ".pdf"),
  plot = last_plot(),
  width = 8, height = 8, units = c("in")
)
