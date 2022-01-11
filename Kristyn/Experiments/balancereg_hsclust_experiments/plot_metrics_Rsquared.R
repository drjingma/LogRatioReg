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
slrhsc_sims_list = list()
slrhsc2_sims_list = list()
slrhsc_natstop_sims_list = list()
# slrhsc_ngmstop_sims_list = list()
cl_sims_list = list()
pr_sims_list = list()
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
  # slr hsc - eta - natural stop
  slrhsc_natstop.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_hsc_eta_natstop_", "metrics", file.end0,
    "_sim", i, ".rds"
  ))))
  rownames(slrhsc_natstop.sim.tmp) = NULL
  slrhsc_natstop_sims_list[[i]] = data.table(slrhsc_natstop.sim.tmp)
  # # slr hsc - eta - NGM stop
  # slrhsc_ngmstop.sim.tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/metrics", "/slr_hsc_eta_ngmstop_", "metrics", file.end0,
  #   "_sim", i, ".rds"
  # ))))
  # rownames(slrhsc_ngmstop.sim.tmp) = NULL
  # slrhsc_ngmstop_sims_list[[i]] = data.table(slrhsc_ngmstop.sim.tmp)
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
}

slrhc_sims = as.data.frame(rbindlist(slrhc_sims_list))
slrhc2_sims = as.data.frame(rbindlist(slrhc2_sims_list))
slrhsc_sims = as.data.frame(rbindlist(slrhsc_sims_list))
slrhsc2_sims = as.data.frame(rbindlist(slrhsc2_sims_list))
slrhsc_natstop_sims = as.data.frame(rbindlist(slrhsc_natstop_sims_list))
# slrhsc_ngmstop_sims = as.data.frame(rbindlist(slrhsc_ngmstop_sims_list))
cl_sims = as.data.frame(rbindlist(cl_sims_list))
pr_sims = as.data.frame(rbindlist(pr_sims_list))

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
# slrhsc_ngmstop_summaries = data.frame(
#   "mean" = apply(slrhsc_ngmstop_sims, 2, mean), 
#   "sd" = apply(slrhsc_ngmstop_sims, 2, sd), 
#   "se" =  apply(slrhsc_ngmstop_sims, 2, sd) / sqrt(numSims)
# )
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

# metrics boxplots
slrhc.sims.gg = reshape2::melt(slrhc_sims)
slrhc.sims.gg$Method = "slr-hc"
slrhc2.sims.gg = reshape2::melt(slrhc2_sims)
slrhc2.sims.gg$Method = "slr-hc-eta"
slrhsc.sims.gg = reshape2::melt(slrhsc_sims)
slrhsc.sims.gg$Method = "slr-hsc"
slrhsc2.sims.gg = reshape2::melt(slrhsc2_sims)
slrhsc2.sims.gg$Method = "slr-hsc-eta"
slrhsc_natstop.sims.gg = reshape2::melt(slrhsc_natstop_sims)
slrhsc_natstop.sims.gg$Method = "slr-hsc-eta-nat"
# slrhsc_ngmstop.sims.gg = reshape2::melt(slrhsc_ngmstop_sims)
# slrhsc_ngmstop.sims.gg$Method = "slr-hsc-eta-ngm"
cl.sims.gg = reshape2::melt(cl_sims)
cl.sims.gg$Method = "classo"
pr.sims.gg = reshape2::melt(pr_sims)
pr.sims.gg$Method = "propr"

data.gg0 = rbind(
  slrhc.sims.gg, 
  # slrhc2.sims.gg, # slr-hc-eta -- not good
  slrhsc.sims.gg, # slr-hsc -- not good, if not thresholding
  slrhsc2.sims.gg, 
  slrhsc_natstop.sims.gg, 
  # slrhsc_ngmstop.sims.gg,
  cl.sims.gg, 
  pr.sims.gg)
# levels.gg = c(
#   "slr-hc", "slr-hsc", "slr-hc-eta", "slr-hsc-eta", "classo", "propr")
data.gg = data.gg0
if(!is.null(metric_names)){
  data.gg = dplyr::filter(data.gg, variable %in% metric_names)
}
data.gg$Method = factor(data.gg$Method)#, levels = levels.gg)
data.gg = data.gg %>% filter(
  variable %in% c(
    "PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR", "precision", 
    "Fscore", "time"
  )
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
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.title.y = element_blank())

ggsave(
  filename = paste0(
    "20220111",
    "_Rsq", desired_Rsquared,
    "_rho", rho, 
    "_", "metrics", ".pdf"),
  plot = last_plot(),
  width = 8, height = 5, units = c("in")
)



# ################################################################################
# # plot roc curves
# 
# # import roc curves and organize TPR, S.hat, lambda information
# # cl
# cl.roc.list = list()
# # each row corresponds to a lambda in the lambda sequence (different in ea. sim)
# # each column corresponds to a different simulation
# cl.TPR.mat = matrix(NA, nlam, numSims) 
# cl.S.hat.mat = matrix(NA, nlam, numSims)
# cl.TP.mat = matrix(NA, nlam, numSims)
# cl.prec.mat = matrix(NA, nlam, numSims)
# # slr
# slr.roc.list = list()
# slr.TPR.mat = matrix(NA, nlam, numSims)
# slr.S.hat.mat = matrix(NA, nlam, numSims)
# slr.TP.mat = matrix(NA, nlam, numSims)
# slr.prec.mat = matrix(NA, nlam, numSims)
# if(has.oracle){
#   # oracle
#   or.roc.list = list()
#   or.TPR.mat = matrix(NA, nlam, numSims)
#   or.S.hat.mat = matrix(NA, nlam, numSims)
#   or.TP.mat = matrix(NA, nlam, numSims)
#   or.prec.mat = matrix(NA, nlam, numSims)
# }
# if(has.coat){
#   # coat
#   coat.roc.list = list()
#   coat.TPR.mat = matrix(NA, nlam, numSims)
#   coat.S.hat.mat = matrix(NA, nlam, numSims)
#   coat.TP.mat = matrix(NA, nlam, numSims)
#   coat.prec.mat = matrix(NA, nlam, numSims)
# }
# if(has.propr){
#   # coat
#   pr.roc.list = list()
#   pr.TPR.mat = matrix(NA, nlam, numSims)
#   pr.S.hat.mat = matrix(NA, nlam, numSims)
#   pr.TP.mat = matrix(NA, nlam, numSims)
#   pr.prec.mat = matrix(NA, nlam, numSims)
# }
# for(i in 1:numSims){
#   # cl
#   cl.sim.tmp = readRDS(paste0(
#     output_dir, "/roccurves", "/classo_roc", file.end, "_sim", i, ".rds"
#   ))
#   cl.roc.list[[i]] = cl.sim.tmp
#   cl.TPR.mat[, i] = cl.sim.tmp["tpr", ]
#   cl.S.hat.mat[, i] = cl.sim.tmp["S_hat", ]
#   cl.TP.mat[, i] = cl.sim.tmp["TP", ]
#   cl.prec.mat[, i] = cl.sim.tmp["precision", ]
#   # slr
#   slr.sim.tmp = readRDS(paste0(
#     output_dir, "/roccurves", "/slr_roc", file.end, "_sim", i, ".rds"
#   ))
#   slr.roc.list[[i]] = slr.sim.tmp
#   slr.TPR.mat[, i] = slr.sim.tmp["tpr", ]
#   slr.S.hat.mat[, i] = slr.sim.tmp["S_hat", ]
#   slr.TP.mat[, i] = slr.sim.tmp["TP", ]
#   slr.prec.mat[, i] = slr.sim.tmp["precision", ]
#   if(has.oracle){
#     # oracle
#     or.sim.tmp = readRDS(paste0(
#       output_dir, "/roccurves", "/oracle_roc", file.end, "_sim", i, ".rds"
#     ))
#     or.roc.list[[i]] = or.sim.tmp
#     or.TPR.mat[, i] = or.sim.tmp["tpr", ]
#     or.S.hat.mat[, i] = or.sim.tmp["S_hat", ]
#     or.TP.mat[, i] = or.sim.tmp["TP", ]
#     or.prec.mat[, i] = or.sim.tmp["precision", ]
#   }
#   if(has.coat){
#     # coat
#     coat.sim.tmp = readRDS(paste0(
#       output_dir, "/roccurves", "/coat_roc", file.end, "_sim", i, ".rds"
#     ))
#     coat.roc.list[[i]] = coat.sim.tmp
#     coat.TPR.mat[, i] = coat.sim.tmp["tpr", ]
#     coat.S.hat.mat[, i] = coat.sim.tmp["S_hat", ]
#     coat.TP.mat[, i] = coat.sim.tmp["TP", ]
#     coat.prec.mat[, i] = coat.sim.tmp["precision", ]
#   }
#   if(has.propr){
#     # oracle
#     pr.sim.tmp = readRDS(paste0(
#       output_dir, "/roccurves", "/propr_roc", file.end, "_sim", i, ".rds"
#     ))
#     pr.roc.list[[i]] = pr.sim.tmp
#     pr.TPR.mat[, i] = pr.sim.tmp["tpr", ]
#     pr.S.hat.mat[, i] = pr.sim.tmp["S_hat", ]
#     pr.TP.mat[, i] = pr.sim.tmp["TP", ]
#     pr.prec.mat[, i] = pr.sim.tmp["precision", ]
#   }
# }
# 
# # plot one sim
# idx = 1
# ggdat0 = data.frame(
#   Method = c(
#     rep("classo", nlam), rep("slr", nlam), rep("oracle", nlam),
#     rep("propr", nlam)
#   ),
#   S_hat = c(
#     cl.S.hat.mat[, idx], slr.S.hat.mat[, idx], or.S.hat.mat[, idx],
#     pr.S.hat.mat[, idx]),
#   TPR = c(
#     cl.TPR.mat[, idx], slr.TPR.mat[, idx], or.TPR.mat[, idx],
#     pr.TPR.mat[, idx])
# )
# ggplot(ggdat0, aes(x = S_hat, y = TPR, color = Method)) +
#   geom_line()
# 
# # average over each possible S.hat/TP value
# # stack columns so which() is more interpretable
# cl.TPR.vec = as.vector(cl.TPR.mat)
# cl.S.hat.vec = as.vector(cl.S.hat.mat)
# cl.TP.vec = as.vector(cl.TP.mat)
# cl.prec.vec = as.vector(cl.prec.mat)
# slr.TPR.vec = as.vector(slr.TPR.mat)
# slr.S.hat.vec = as.vector(slr.S.hat.mat)
# slr.TP.vec = as.vector(slr.TP.mat)
# slr.prec.vec = as.vector(slr.prec.mat)
# if(has.oracle){
#   or.TPR.vec = as.vector(or.TPR.mat)
#   or.S.hat.vec = as.vector(or.S.hat.mat)
#   or.TP.vec = as.vector(or.TP.mat)
#   or.prec.vec = as.vector(or.prec.mat)
# }
# if(has.coat){
#   coat.TPR.vec = as.vector(coat.TPR.mat)
#   coat.S.hat.vec = as.vector(coat.S.hat.mat)
#   coat.TP.vec = as.vector(coat.TP.mat)
#   coat.prec.vec = as.vector(coat.prec.mat)
# }
# if(has.propr){
#   pr.TPR.vec = as.vector(pr.TPR.mat)
#   pr.S.hat.vec = as.vector(pr.S.hat.mat)
#   pr.TP.vec = as.vector(pr.TP.mat)
#   pr.prec.vec = as.vector(pr.prec.mat)
# }
# 
# # get the averages
# S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec)))
# if(has.oracle & !has.coat & !has.propr){
#   S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec)))
# } else if(has.oracle & !has.coat & has.propr){
#   S.hat.vals = sort(unique(c(
#     cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec, pr.S.hat.vec)))
# }
# cl.TPR.avg = rep(NA, length(S.hat.vals))
# cl.TP.avg = rep(NA, length(S.hat.vals))
# cl.prec.avg = rep(NA, length(S.hat.vals))
# slr.TPR.avg = rep(NA, length(S.hat.vals))
# slr.TP.avg = rep(NA, length(S.hat.vals))
# slr.prec.avg = rep(NA, length(S.hat.vals))
# if(has.oracle){
#   or.TPR.avg = rep(NA, length(S.hat.vals))
#   or.TP.avg = rep(NA, length(S.hat.vals))
#   or.prec.avg = rep(NA, length(S.hat.vals))
# }
# if(has.coat){
#   coat.TPR.avg = rep(NA, length(S.hat.vals))
#   coat.TP.avg = rep(NA, length(S.hat.vals))
#   coat.prec.avg = rep(NA, length(S.hat.vals))
# }
# if(has.propr){
#   pr.TPR.avg = rep(NA, length(S.hat.vals))
#   pr.TP.avg = rep(NA, length(S.hat.vals))
#   pr.prec.avg = rep(NA, length(S.hat.vals))
# }
# for(i in 1:length(S.hat.vals)){
#   val.tmp = S.hat.vals[i]
#   # classo
#   cl.which.idx.tmp = which(cl.S.hat.vec == val.tmp)
#   cl.TPR.avg[i] = mean(cl.TPR.vec[cl.which.idx.tmp])
#   cl.TP.avg[i] = mean(cl.TP.vec[cl.which.idx.tmp])
#   cl.prec.avg[i] = mean(cl.prec.vec[cl.which.idx.tmp])
#   # slr
#   slr.which.idx.tmp = which(slr.S.hat.vec == val.tmp)
#   slr.TPR.avg[i] = mean(slr.TPR.vec[slr.which.idx.tmp])
#   slr.TP.avg[i] = mean(slr.TP.vec[slr.which.idx.tmp])
#   slr.prec.avg[i] = mean(slr.prec.vec[slr.which.idx.tmp])
#   if(has.oracle){
#     # oracle
#     or.which.idx.tmp = which(or.S.hat.vec == val.tmp)
#     or.TPR.avg[i] = mean(or.TPR.vec[or.which.idx.tmp])
#     or.TP.avg[i] = mean(or.TP.vec[or.which.idx.tmp])
#     or.prec.avg[i] = mean(or.prec.vec[or.which.idx.tmp])
#   }
#   if(has.coat){
#     # coat
#     coat.which.idx.tmp = which(coat.S.hat.vec == val.tmp)
#     coat.TPR.avg[i] = mean(coat.TPR.vec[coat.which.idx.tmp])
#     coat.TP.avg[i] = mean(coat.TP.vec[coat.which.idx.tmp])
#     coat.prec.avg[i] = mean(coat.prec.vec[coat.which.idx.tmp])
#   }
#   if(has.propr){
#     # propr
#     pr.which.idx.tmp = which(pr.S.hat.vec == val.tmp)
#     pr.TPR.avg[i] = mean(pr.TPR.vec[pr.which.idx.tmp])
#     pr.TP.avg[i] = mean(pr.TP.vec[pr.which.idx.tmp])
#     pr.prec.avg[i] = mean(pr.prec.vec[pr.which.idx.tmp])
#   }
# }
# 
# # plot
# data.gg = rbind(
#   data.frame(
#     S_hat = S.hat.vals, TPR = cl.TPR.avg, TP = cl.TP.avg, 
#     precision = cl.prec.avg, Method = "classo"), 
#   data.frame(
#     S_hat = S.hat.vals, TPR = slr.TPR.avg, TP = slr.TP.avg, 
#     precision = slr.prec.avg, Method = "slr")
# )
# if(has.oracle){
#   data.gg = rbind(
#     data.gg, 
#     data.frame(
#       S_hat = S.hat.vals, TPR = or.TPR.avg, TP = or.TP.avg, 
#       precision = or.prec.avg, Method = "oracle"))
# }
# if(has.coat){
#   data.gg = rbind(
#     data.gg, 
#     data.frame(
#       S_hat = S.hat.vals, TPR = coat.TPR.avg, TP = coat.TP.avg, 
#       precision = coat.prec.avg, Method = "coat"))
# }
# if(has.propr){
#   data.gg = rbind(
#     data.gg, 
#     data.frame(
#       S_hat = S.hat.vals, TPR = pr.TPR.avg, TP = pr.TP.avg, 
#       precision = pr.prec.avg, Method = "propr"))
# }
# data.gg$Method = factor(data.gg$Method, levels = levels.gg)
# tp_roc = ggplot(
#   data.gg[!is.na(data.gg$TP),], 
#   aes(x = S_hat, y = TP, color = Method, linetype = Method)) + 
#   geom_line(alpha = 0.75, size = 1, na.rm = TRUE) +
#   # geom_point(alpha = 0.5, na.rm = TRUE) +
#   theme_bw()
# tpr_roc = ggplot(
#   data.gg[!is.na(data.gg$TPR),], 
#   aes(x = S_hat, y = TPR, color = Method, linetype = Method)) + 
#   geom_line(alpha = 0.75, size = 1, na.rm = TRUE) +
#   # geom_point(alpha = 0.5, na.rm = TRUE) +
#   theme_bw()
# prec_roc = ggplot(
#   data.gg[!is.na(data.gg$precision),], 
#   aes(x = S_hat, y = precision, color = Method, linetype = Method)) + 
#   geom_line(alpha = 0.75, size = 1, na.rm = TRUE) +
#   # geom_point(alpha = 0.5, na.rm = TRUE) +
#   theme_bw()
# tpr_roc
# # ggarrange(tp_roc, tpr_roc, prec_roc, nrow = 1)
# 
# ggsave(
#   filename = paste0(
#     "20211202_", 
#     sigma.settings, "_noise", sigma_eps, 
#     "_", theta.settings, 
#     "_val", values.theta,
#     "_rocs.pdf"),
#   plot = last_plot(),
#   width = 6, height = 5, units = c("in")
# )
# 
# 
