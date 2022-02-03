rm(list=ls())
# Purpose: demonstrate hierarchical spectral clustering with a threshold
#   explore various sigma_eps & rho values to get specified Rsquared values
# Date: 1/3/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_hsclust_experiments/outputs"

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100

# Settings to toggle with
sigma.settings = "expdecaySigma"
values.theta = 1
linkage = "average"
tol = 1e-4
nlam = 100
neta = 50
intercept = TRUE
K = 10
n = 100
p = 30
scaling = TRUE
#################
# if rho = 0, 
#   sigma_eps = sqrt(2/3) => R^2 = 0.6
#   sigma_eps = sqrt(1/4) => R^2 = 0.8
# if rho = 0.2, 
#   sigma_eps = sqrt(0.7125333) => R^2 = 0.6
#   sigma_eps = sqrt(0.2672) => R^2 = 0.8
# if rho = 0.5, 
#   sigma_eps = sqrt(0.808333) => R^2 = 0.6
#   sigma_eps = sqrt(0.303125) => R^2 = 0.8
get_sigma_eps = function(theta_val, Rsq_val, rho_val){
  sigma_eps_sq.tmp = theta_val^2 * (1 - Rsq_val) / Rsq_val + 
    theta_val^2 * (1 - Rsq_val) * (rho_val^3 + 2 * rho_val^2 + 3 * rho_val) / 
    (10 * Rsq_val)
  return(sqrt(sigma_eps_sq.tmp))
}
rho = 0.2 #
desired_Rsquared = 0.6 #
sigma_eps = get_sigma_eps(
  theta_val = values.theta, Rsq_val = desired_Rsquared, rho_val = rho)

file.end0 = paste0(
  "_", sigma.settings,
  "_dim", n, "x", p, 
  "_Rsq", desired_Rsquared,
  "_rho", rho)

################################################################################
# plot metrics

metric_names = c(
  "PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR", "precision",
  "Fscore", "betaSparsity", "Rsq", "time")

# import metrics
slrhc_sims_list = list()
slrhc2_sims_list = list()
slrhc_distal_sims_list = list()
slrhsc_sims_list = list()
slrhsc2_sims_list = list()
# slr2hsc2_sims_list = list()
slrhsc_natstop_sims_list = list()
slrhsc_ngmstop_sims_list = list()
cl_sims_list = list()
pr_sims_list = list()
or_sims_list = list()
# slbl_sims_list = list()
for(i in 1:numSims){
  print(i)
  # slr hc
  slrhc.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hc_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhc.sim.tmp) = NULL
  slrhc_sims_list[[i]] = data.table(slrhc.sim.tmp)
  # slr hc - eta
  slrhc2.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hc_eta_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhc2.sim.tmp) = NULL
  slrhc2_sims_list[[i]] = data.table(slrhc2.sim.tmp)
  # slr hc - eta
  slrhc_distal.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hc_distal_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhc_distal.sim.tmp) = NULL
  slrhc_distal_sims_list[[i]] = data.table(slrhc_distal.sim.tmp)
  # slr hsc
  slrhsc.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hsc_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhsc.sim.tmp) = NULL
  slrhsc_sims_list[[i]] = data.table(slrhsc.sim.tmp)
  # slr hsc - eta
  slrhsc2.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hsc_eta_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhsc2.sim.tmp) = NULL
  slrhsc2_sims_list[[i]] = data.table(slrhsc2.sim.tmp)
  # # slr2 hsc - eta
  # slr2hsc2.sim.tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/metrics", "/slr2_hsc_eta_", "metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(slr2hsc2.sim.tmp) = NULL
  # slr2hsc2_sims_list[[i]] = data.table(slr2hsc2.sim.tmp)
  # slr hsc - eta - natural stop
  slrhsc_natstop.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hsc_eta_natstop_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhsc_natstop.sim.tmp) = NULL
  slrhsc_natstop_sims_list[[i]] = data.table(slrhsc_natstop.sim.tmp)
  # # slr hsc - eta - NGM stop
  slrhsc_ngmstop.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hsc_eta_ngmstop_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhsc_ngmstop.sim.tmp) = NULL
  slrhsc_ngmstop_sims_list[[i]] = data.table(slrhsc_ngmstop.sim.tmp)
  # classo
  cl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/classo_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(cl.sim.tmp) = NULL
  cl_sims_list[[i]] = data.table(cl.sim.tmp)
  # propr
  pr.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/propr_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(pr.sim.tmp) = NULL
  pr_sims_list[[i]] = data.table(pr.sim.tmp)
  # oracle
  or.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/oracle_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(or.sim.tmp) = NULL
  or_sims_list[[i]] = data.table(or.sim.tmp)
  # # selbal
  # slbl.sim.tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/metrics", "/selbal_", "metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(slbl.sim.tmp) = NULL
  # slbl_sims_list[[i]] = data.table(slbl.sim.tmp)
}

slrhc_sims = as.data.frame(rbindlist(slrhc_sims_list))
slrhc2_sims = as.data.frame(rbindlist(slrhc2_sims_list))
slrhc_distal_sims = as.data.frame(rbindlist(slrhc_distal_sims_list))
slrhsc_sims = as.data.frame(rbindlist(slrhsc_sims_list))
slrhsc2_sims = as.data.frame(rbindlist(slrhsc2_sims_list))
slrhsc_natstop_sims = as.data.frame(rbindlist(slrhsc_natstop_sims_list))
slrhsc_ngmstop_sims = as.data.frame(rbindlist(slrhsc_ngmstop_sims_list))
cl_sims = as.data.frame(rbindlist(cl_sims_list))
pr_sims = as.data.frame(rbindlist(pr_sims_list))
or_sims = as.data.frame(rbindlist(or_sims_list))
slbl_sims = as.data.frame(rbindlist(slbl_sims_list))

# summary stats
slrhc_summaries = data.frame(
  "mean" = apply(slrhc_sims, 2, mean), 
  "sd" = apply(slrhc_sims, 2, sd), 
  "se" =  apply(slrhc_sims, 2, sd) / sqrt(numSims)
)
slrhc2_summaries = data.frame(
  "mean" = apply(slrhc2_sims, 2, mean), 
  "sd" = apply(slrhc2_sims, 2, sd), 
  "se" =  apply(slrhc2_sims, 2, sd) / sqrt(numSims)
)
slrhc_distal_summaries = data.frame(
  "mean" = apply(slrhc_distal_sims, 2, mean), 
  "sd" = apply(slrhc_distal_sims, 2, sd), 
  "se" =  apply(slrhc_distal_sims, 2, sd) / sqrt(numSims)
)
slrhsc_summaries = data.frame(
  "mean" = apply(slrhsc_sims, 2, mean), 
  "sd" = apply(slrhsc_sims, 2, sd), 
  "se" =  apply(slrhsc_sims, 2, sd) / sqrt(numSims)
)
slrhsc2_summaries = data.frame(
  "mean" = apply(slrhsc2_sims, 2, mean), 
  "sd" = apply(slrhsc2_sims, 2, sd), 
  "se" =  apply(slrhsc2_sims, 2, sd) / sqrt(numSims)
)
slrhsc_natstop_summaries = data.frame(
  "mean" = apply(slrhsc_natstop_sims, 2, mean), 
  "sd" = apply(slrhsc_natstop_sims, 2, sd), 
  "se" =  apply(slrhsc_natstop_sims, 2, sd) / sqrt(numSims)
)
slrhsc_ngmstop_summaries = data.frame(
  "mean" = apply(slrhsc_ngmstop_sims, 2, mean),
  "sd" = apply(slrhsc_ngmstop_sims, 2, sd),
  "se" =  apply(slrhsc_ngmstop_sims, 2, sd) / sqrt(numSims)
)
cl_summaries = data.frame(
  "mean" = apply(cl_sims, 2, mean), 
  "sd" = apply(cl_sims, 2, sd), 
  "se" =  apply(cl_sims, 2, sd) / sqrt(numSims)
)
pr_summaries = data.frame(
  "mean" = apply(pr_sims, 2, mean), 
  "sd" = apply(pr_sims, 2, sd), 
  "se" =  apply(pr_sims, 2, sd) / sqrt(numSims)
)
or_summaries = data.frame(
  "mean" = apply(or_sims, 2, mean), 
  "sd" = apply(or_sims, 2, sd), 
  "se" =  apply(or_sims, 2, sd) / sqrt(numSims)
)
slbl_summaries = data.frame(
  "mean" = apply(slbl_sims, 2, mean), 
  "sd" = apply(slbl_sims, 2, sd), 
  "se" =  apply(slbl_sims, 2, sd) / sqrt(numSims)
)

# metrics boxplots
slrhc.sims.gg = reshape2::melt(slrhc_sims)
slrhc.sims.gg$Method = "slr-hc"
slrhc2.sims.gg = reshape2::melt(slrhc2_sims)
slrhc2.sims.gg$Method = "slr-hc-eta"
slrhc_distal.sims.gg = reshape2::melt(slrhc_distal_sims)
slrhc_distal.sims.gg$Method = "slr-hc-distal"
slrhsc.sims.gg = reshape2::melt(slrhsc_sims)
slrhsc.sims.gg$Method = "slr-hsc"
slrhsc2.sims.gg = reshape2::melt(slrhsc2_sims)
slrhsc2.sims.gg$Method = "slr-hsc-eta"
slrhsc_natstop.sims.gg = reshape2::melt(slrhsc_natstop_sims)
slrhsc_natstop.sims.gg$Method = "slr-hsc-eta-nat"
slrhsc_ngmstop.sims.gg = reshape2::melt(slrhsc_ngmstop_sims)
slrhsc_ngmstop.sims.gg$Method = "slr-hsc-eta-ngm"
cl.sims.gg = reshape2::melt(cl_sims)
cl.sims.gg$Method = "classo"
pr.sims.gg = reshape2::melt(pr_sims)
pr.sims.gg$Method = "propr"
or.sims.gg = reshape2::melt(or_sims)
or.sims.gg$Method = "oracle (hc)"
slbl.sims.gg = reshape2::melt(slbl_sims)
slbl.sims.gg$Method = "selbal"

data.gg0 = rbind(
  slrhc.sims.gg, 
  # slrhc2.sims.gg, # slr-hc-eta -- not good
  slrhc_distal.sims.gg,
  slrhsc.sims.gg, # slr-hsc -- not good, if not thresholding
  slrhsc2.sims.gg, 
  slrhsc_natstop.sims.gg, 
  slrhsc_ngmstop.sims.gg,
  cl.sims.gg, 
  pr.sims.gg, 
  or.sims.gg, 
  slbl.sims.gg)
# levels.gg = c(
#   "slr-hc", "slr-hsc", "slr-hc-eta", "slr-hsc-eta", "classo", "propr")
data.gg = data.gg0
if(!is.null(metric_names)){
  data.gg = dplyr::filter(data.gg, variable %in% metric_names)
}
# data.gg$Method = factor(data.gg$Method, levels = levels.gg)
data.gg = data.gg %>% dplyr::filter(
  variable %in% c(
    "PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR", "precision", 
    "Fscore", "time"
  )
) %>% dplyr::filter(
  Method != "selbal"
)

ggplot(data.gg, aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(variable), scales = "free_y") +
  geom_boxplot() +
  stat_summary(
    fun = mean, fun.min = mean, fun.max = mean,
    geom = "errorbar", width = 0.75,
    linetype = "dashed") +
  stat_summary(
    fun = mean, geom = "point", shape = 17, size = 2,
    color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank())

ggsave(
  filename = paste0(
    "20220125",
    "_Rsq", desired_Rsquared,
    "_rho", rho, 
    "_", "metrics", ".pdf"),
  plot = last_plot(),
  width = 8, height = 5, units = c("in")
)


