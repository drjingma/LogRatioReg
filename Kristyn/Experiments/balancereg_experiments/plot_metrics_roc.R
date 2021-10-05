rm(list=ls())
# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Date: 10/04/2021

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_experiments/outputs"

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100
rng.seed = 123

# Settings to toggle with
sigma.settings = "10blockSigma" # 2blockSigma, 4blockSigma, 10blockSigma, lin14Sigma
rho.type = "square" # 1 = "absolute value", 2 = "square"
theta.settings = "pairperblock" # "dense", "sparse", "both", "multsparse"
# if "4blockSigma", then "2blocks", "1block", "2blocks2contrasts"
# if "2blockSigma" then "dense"
# if "10blockSigma", then "pairperblock"
# if "lin14Sigma" then "sparse" or "dense"
linkage = "average"
tol = 1e-4
nlam = 200
intercept = TRUE
K = 10
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
cor_ij = 0.2 # 0.2, 0.5
scaling = TRUE
sigma_eps = 0.5 # 0.1, 0.5

if(sigma.settings == "lin14Sigma"){
  file.end = paste0( # for old simulations
    "_dim", n, "x", p,
    "_", sigma.settings,
    "_", theta.settings,
    "_noise", sigma_eps,
    "_rho", rho,
    "_int", intercept,
    "_scale", scaling,
    "_K", K,
    "_seed", rng.seed,
    ".rds")
} else{ # for block-diagonal Sigma, either "2blockSigma" or "4blockSigma"
  file.end = paste0(
    "_dim", n, "x", p, 
    "_", sigma.settings,
    "_", theta.settings, 
    "_noise", sigma_eps,
    "_cor", cor_ij, 
    "_int", intercept,
    "_scale", scaling,
    "_K", K,
    "_seed", rng.seed,
    ".rds")
}

has.selbal = FALSE
has.coat = FALSE
has.oracle = TRUE
if(FALSE){
  has.selbal = TRUE
}
if(sigma.settings == "10blockSigma"){
  has.coat = TRUE
}
# if(theta.settings %in% c("block", "pairblock")){
#   has.oracle = TRUE
# }

################################################################################
# plot metrics

metric_names = NULL
# metrics_names = c(
#   "PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR", "timing", 
#   "betaSparsity")
# metrics_names = c(
#   "PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR", "timing")

# import metrics
cl.sims.list = list()
slr.sims.list = list()
if(has.selbal) selbal.sims.list = list()
if(has.oracle) or.sims.list = list()
if(has.coat) coat.sims.list = list()
for(i in 1:numSims){
  # classo
  cl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo_metrics", i, file.end
  ))))
  rownames(cl.sim.tmp) = NULL
  cl.sims.list[[i]] = data.table(cl.sim.tmp)
  # slr
  slr.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_metrics", i, file.end
  ))))
  rownames(slr.sim.tmp) = NULL
  slr.sims.list[[i]] = data.table(slr.sim.tmp)
  if(has.selbal){
    # selbal
    selbal.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/selbal_metrics", i, file.end
    ))))
    rownames(selbal.sim.tmp) = NULL
    or.sims.list[[i]] = data.table(selbal.sim.tmp)
  }
  if(has.oracle){
    # oracle
    or.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/oracle_metrics", i, file.end
    ))))
    rownames(or.sim.tmp) = NULL
    or.sims.list[[i]] = data.table(or.sim.tmp)
  }
  if(has.coat){
    # coat
    coat.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/oracle_metrics", i, file.end
    ))))
    rownames(coat.sim.tmp) = NULL
    coat.sims.list[[i]] = data.table(coat.sim.tmp)
  }
}
cl.sims = as.data.frame(rbindlist(cl.sims.list))
slr.sims = as.data.frame(rbindlist(slr.sims.list))
if(has.selbal) selbal.sims = as.data.frame(rbindlist(selbal.sims.list))
if(has.oracle) or.sims = as.data.frame(rbindlist(or.sims.list))
if(has.oracle) coat.sims = as.data.frame(rbindlist(coat.sims.list))

# summary stats for classo metrics
cl.eval.means = apply(cl.sims, 2, mean)
cl.eval.sds = apply(cl.sims, 2, sd)
cl.eval.ses = cl.eval.sds / sqrt(numSims)
cl.summaries = data.frame(
  "mean" = cl.eval.means, "sd" = cl.eval.sds, "se" = cl.eval.ses)
# print(cl.summaries[metrics, c("mean", "se")])

# summary stats for slr metrics
slr.eval.means = apply(slr.sims, 2, mean)
slr.eval.sds = apply(slr.sims, 2, sd)
slr.eval.ses = slr.eval.sds / sqrt(numSims)
slr.summaries = data.frame(
  "mean" = slr.eval.means, "sd" = slr.eval.sds, "se" = slr.eval.ses)
# print(slr.summaries[metrics, c("mean", "se")])

if(has.selbal){
  # summary stats for selbal metrics
  selbal.eval.means = apply(selbal.sims, 2, mean)
  selbal.eval.sds = apply(selbal.sims, 2, sd)
  selbal.eval.ses = selbal.eval.sds / sqrt(numSims)
  selbal.summaries = data.frame(
    "mean" = selbal.eval.means, "sd" = selbal.eval.sds, "se" = selbal.eval.ses)
  # print(selbal.summaries[metrics, c("mean", "se")])
}

if(has.oracle){
  # summary stats for oracle metrics
  or.eval.means = apply(or.sims, 2, mean)
  or.eval.sds = apply(or.sims, 2, sd)
  or.eval.ses = or.eval.sds / sqrt(numSims)
  or.summaries = data.frame(
    "mean" = or.eval.means, "sd" = or.eval.sds, "se" = or.eval.ses)
  # print(or.summaries[metrics, c("mean", "se")])
}

if(has.coat){
  # summary stats for oracle metrics
  coat.eval.means = apply(coat.sims, 2, mean)
  coat.eval.sds = apply(coat.sims, 2, sd)
  coat.eval.ses = coat.eval.sds / sqrt(numSims)
  coat.summaries = data.frame(
    "mean" = coat.eval.means, "sd" = coat.eval.sds, "se" = coat.eval.ses)
  # print(coat.summaries[metrics, c("mean", "se")])
}

# boxplots for the slr and classo metrics
cl.sims.gg = reshape2::melt(cl.sims)
cl.sims.gg$Method = "classo"
slr.sims.gg = reshape2::melt(slr.sims)
slr.sims.gg$Method = "slr"
if(has.selbal){
  selbal.sims.gg = reshape2::melt(selbal.sims)
  selbal.sims.gg$Method = "selbal"
}
if(has.oracle){
  or.sims.gg = reshape2::melt(or.sims)
  or.sims.gg$Method = "oracle"
}
if(has.coat){
  coat.sims.gg = reshape2::melt(coat.sims)
  coat.sims.gg$Method = "coat"
}
data.gg = rbind(cl.sims.gg, slr.sims.gg)
levels.gg = c("classo", "slr")
if(has.oracle){
  data.gg = rbind(data.gg, or.sims.gg)
  levels.gg = c(levels.gg, "oracle")
}
if(has.selbal){
  data.gg = rbind(data.gg, selbal.sims.gg)
  levels.gg = c(levels.gg, "selbal")
}
if(has.coat){
  data.gg = rbind(data.gg, coat.sims.gg)
  levels.gg = c(levels.gg, "coat")
}

if(!is.null(metric_names)){
  data.gg = dplyr::filter(data.gg, variable %in% metrics_names)
}
data.gg$Method = factor(data.gg$Method, levels = levels.gg)

ggplot(data.gg, aes(x = Method, y = value, color = Method)) + 
  facet_wrap(vars(variable), scales = "free_y") + 
  geom_boxplot() + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank())
if(sigma.settings == "lin14Sigma"){
  ggsave(
    filename = paste0(
      "20211004_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_metrics.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
} else{
  ggsave(
    filename = paste0(
      "20211004_", 
      sigma.settings, "_noise", sigma_eps,
      "_", theta.settings, "_metrics.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
}

# zoom in to PEtr, PEte
# plot PEtr
PEtr.gg = data.gg[data.gg$variable == "PEtr", ]
# PEtr.gg$value = log(PEtr.gg$value)
plt.PEtr = ggplot(PEtr.gg, aes(x = type, y = value, color = type)) + 
  geom_boxplot() +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")
# plot PEte
PEte.gg = data.gg[data.gg$variable == "PEte", ]
PEte.gg$value = log(PEte.gg$value)
plt.PEte = ggplot(PEte.gg, aes(x = type, y = value, color = type)) + 
  geom_boxplot() +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")
# plot both
# ggarrange(plt.PEtr, plt.PEte)



################################################################################
# plot roc curves

# import roc curves and organize TPR, S.hat, lambda information
# cl
cl.roc.list = list()
# each row corresponds to a lambda in the lambda sequence (different in ea. sim)
# each column corresponds to a different simulation
cl.TPR.mat = matrix(NA, nlam, numSims) 
cl.S.hat.mat = matrix(NA, nlam, numSims)
# slr
slr.roc.list = list()
slr.TPR.mat = matrix(NA, nlam, numSims)
slr.S.hat.mat = matrix(NA, nlam, numSims)
if(has.oracle){
  # oracle
  or.roc.list = list()
  or.TPR.mat = matrix(NA, nlam, numSims)
  or.S.hat.mat = matrix(NA, nlam, numSims)
}
for(i in 1:numSims){
  # cl
  cl.sim.tmp = readRDS(paste0(
    output_dir, "/classo_roc", i, file.end
  ))
  cl.roc.list[[i]] = cl.sim.tmp
  cl.TPR.mat[, i] = cl.sim.tmp["tpr", ]
  cl.S.hat.mat[, i] = cl.sim.tmp["S_hat", ]
  # slr
  slr.sim.tmp = readRDS(paste0(
    output_dir, "/slr_roc", i, file.end
  ))
  slr.roc.list[[i]] = slr.sim.tmp
  slr.TPR.mat[, i] = slr.sim.tmp["tpr", ]
  slr.S.hat.mat[, i] = slr.sim.tmp["S_hat", ]
  if(has.oracle){
    # oracle
    or.sim.tmp = readRDS(paste0(
      output_dir, "/oracle_roc", i, file.end
    ))
    or.roc.list[[i]] = or.sim.tmp
    or.TPR.mat[, i] = or.sim.tmp["tpr", ]
    or.S.hat.mat[, i] = or.sim.tmp["S_hat", ]
  }
}

# average over each possible S.hat value
# stack columns so which() is more interpretable
cl.TPR.vec = as.vector(cl.TPR.mat)
cl.S.hat.vec = as.vector(cl.S.hat.mat)
slr.TPR.vec = as.vector(slr.TPR.mat)
slr.S.hat.vec = as.vector(slr.S.hat.mat)
if(has.oracle){
  or.TPR.vec = as.vector(or.TPR.mat)
  or.S.hat.vec = as.vector(or.S.hat.mat)
}

# get the averages
if(has.oracle){
  S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec)))
} else{
  S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec)))
}
cl.tpr.avg = rep(NA, length(S.hat.vals))
slr.tpr.avg = rep(NA, length(S.hat.vals))
if(has.oracle) or.tpr.avg = rep(NA, length(S.hat.vals))
for(i in 1:length(S.hat.vals)){
  val.tmp = S.hat.vals[i]
  # classo
  cl.which.idx.tmp = which(cl.S.hat.vec == val.tmp)
  cl.tpr.avg[i] = mean(cl.TPR.vec[cl.which.idx.tmp])
  # slr
  slr.which.idx.tmp = which(slr.S.hat.vec == val.tmp)
  slr.tpr.avg[i] = mean(slr.TPR.vec[slr.which.idx.tmp])
  # oracle
  or.which.idx.tmp = which(or.S.hat.vec == val.tmp)
  or.tpr.avg[i] = mean(or.TPR.vec[or.which.idx.tmp])
}

# plot
# plot(S.hat.vals, cl.tpr.avg, type = "l", col = 2)
# lines(S.hat.vals, slr.tpr.avg, col = 3)
data.gg = rbind(
  data.frame(S_hat = S.hat.vals, TPR = cl.tpr.avg, Method = "classo"), 
  data.frame(S_hat = S.hat.vals, TPR = slr.tpr.avg, Method = "slr")
)
if(has.oracle){
  data.gg = rbind(
    data.gg, 
    data.frame(S_hat = S.hat.vals, TPR = or.tpr.avg, Method = "oracle"))
}
data.gg$Method = factor(data.gg$Method, levels = levels.gg)
ggplot(data.gg[!is.na(data.gg$TPR),], aes(x = S_hat, y = TPR, color = Method)) + 
  geom_line(alpha = 0.5, na.rm = TRUE) +
  geom_point(alpha = 0.5, na.rm = TRUE) +
  # xlim(0, 40) +
  theme_bw()
if(sigma.settings == "lin14Sigma"){
  ggsave(
    filename = paste0(
      "20211004_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_roc.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
} else{
  ggsave(
    filename = paste0(
      "20211004_", 
      sigma.settings, "_noise", sigma_eps, 
      "_", theta.settings, "_roc.pdf"),
    plot = last_plot(),
    width = 8, height = 5, units = c("in")
  )
}


