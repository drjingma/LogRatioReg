
rm(list=ls())
# Purpose: Simulate data from balance regression model to compare
#   compositional lasso and supervised log-ratios methods
# Date: 10/11/2021

values.theta = 10 # NULL, 5, 10

################################################################################
# libraries and settings

output_dir = "Kristyn/Experiments/balancereg_hsclust_kmeans_experiments/outputs"

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100
rng.seed = 123

# Settings to toggle with
sigma.settings = "10blockSigma" # 2blockSigma, 4blockSigma, 10blockSigma, lin14Sigma
rho.type = "square" # 1 = "absolute value", 2 = "square"
theta.settings = "1blockpair4halves" # "dense" or "sparse"
# if "2blockSigma" then "dense"
# if "4blockSigma", then "1blockpair"
# if "10blockSigma", then "pairperblock" or "1blockpair4halves"
# if "lin14Sigma" then "dense" or "sparse"
linkage = "average"
tol = 1e-4
nlam = 200
intercept = TRUE
K = 10
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
scaling = TRUE
sigma_eps = 0.01  # 0.01, 0.1, 0.5

if(!is.null(values.theta)){
  theta.settings2 = paste0(theta.settings, "_val",values.theta)
} else{
  theta.settings2 = theta.settings
}

file.end = paste0(
  "_", sigma.settings,
  "_", theta.settings2, 
  "_dim", n, "x", p, 
  "_noise", sigma_eps,
  "_rho", rho, 
  "_int", intercept,
  "_scale", scaling
)

has.selbal = FALSE
has.coat = FALSE
has.oracle = TRUE
has.propr = TRUE

################################################################################
# plot metrics

metric_names = NULL
if(theta.settings == "dense"){
  metric_names = c(
    "PEtr", "PEte", "EA1", "EA2", "EAInfty", "FP", "FN", "TPR", "precision",
    "Fscore", "timing", "betaSparsity")
}
metrics.file = "metrics"

# import metrics
cl.sims.list = list()
slr.sims.list = list()
if(has.selbal) selbal.sims.list = list()
if(has.oracle) or.sims.list = list()
if(has.coat) coat.sims.list = list()
if(has.propr) pr.sims.list = list()
for(i in 1:numSims){
  # classo
  cl.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/classo_", metrics.file, file.end,
    "_sim", i, ".rds"
  ))))
  rownames(cl.sim.tmp) = NULL
  cl.sims.list[[i]] = data.table(cl.sim.tmp)
  # slr
  slr.sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/metrics", "/slr_", metrics.file, file.end,
    "_sim", i, ".rds"
  ))))
  rownames(slr.sim.tmp) = NULL
  slr.sims.list[[i]] = data.table(slr.sim.tmp)
  if(has.selbal){
    # selbal
    selbal.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/metrics", "/selbal_", metrics.file, file.end,
      "_sim", i, ".rds"
    ))))
    rownames(selbal.sim.tmp) = NULL
    selbal.sims.list[[i]] = data.table(selbal.sim.tmp)
  }
  if(has.oracle){
    # oracle
    or.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/metrics", "/oracle_", metrics.file, file.end,
      "_sim", i, ".rds"
    ))))
    rownames(or.sim.tmp) = NULL
    or.sims.list[[i]] = data.table(or.sim.tmp)
  }
  if(has.coat){
    # coat
    coat.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/metrics", "/coat_", metrics.file, file.end,
      "_sim", i, ".rds"
    ))))
    rownames(coat.sim.tmp) = NULL
    coat.sims.list[[i]] = data.table(coat.sim.tmp)
  }
  if(has.propr){
    # propr
    pr.sim.tmp = t(data.frame(readRDS(paste0(
      output_dir, "/metrics", "/propr_", metrics.file, file.end,
      "_sim", i, ".rds"
    ))))
    rownames(pr.sim.tmp) = NULL
    pr.sims.list[[i]] = data.table(pr.sim.tmp)
  }
}

cl.sims = as.data.frame(rbindlist(cl.sims.list))
slr.sims = as.data.frame(rbindlist(slr.sims.list))
if(has.selbal) selbal.sims = as.data.frame(rbindlist(selbal.sims.list))
if(has.oracle) or.sims = as.data.frame(rbindlist(or.sims.list))
if(has.coat) coat.sims = as.data.frame(rbindlist(coat.sims.list))
if(has.propr) pr.sims = as.data.frame(rbindlist(pr.sims.list))

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

if(has.propr){
  # summary stats for oracle metrics
  pr.eval.means = apply(pr.sims, 2, mean)
  pr.eval.sds = apply(pr.sims, 2, sd)
  pr.eval.ses = pr.eval.sds / sqrt(numSims)
  pr.summaries = data.frame(
    "mean" = pr.eval.means, "sd" = pr.eval.sds, "se" = pr.eval.ses)
  # print(pr.summaries[metrics, c("mean", "se")])
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
if(has.propr){
  pr.sims.gg = reshape2::melt(pr.sims)
  pr.sims.gg$Method = "propr"
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
if(has.propr){
  data.gg = rbind(data.gg, pr.sims.gg)
  levels.gg = c(levels.gg, "propr")
}

if(!is.null(metric_names)){
  data.gg = dplyr::filter(data.gg, variable %in% metric_names)
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

ggsave(
  filename = paste0(
    "20211202_",
    sigma.settings, "_noise", sigma_eps,
    "_", theta.settings, 
    "_val", values.theta,
    "_", metrics.file, ".pdf"),
  plot = last_plot(),
  width = 8, height = 5, units = c("in")
)



################################################################################
# plot roc curves

# import roc curves and organize TPR, S.hat, lambda information
# cl
cl.roc.list = list()
# each row corresponds to a lambda in the lambda sequence (different in ea. sim)
# each column corresponds to a different simulation
cl.TPR.mat = matrix(NA, nlam, numSims) 
cl.S.hat.mat = matrix(NA, nlam, numSims)
cl.TP.mat = matrix(NA, nlam, numSims)
cl.prec.mat = matrix(NA, nlam, numSims)
# slr
slr.roc.list = list()
slr.TPR.mat = matrix(NA, nlam, numSims)
slr.S.hat.mat = matrix(NA, nlam, numSims)
slr.TP.mat = matrix(NA, nlam, numSims)
slr.prec.mat = matrix(NA, nlam, numSims)
if(has.oracle){
  # oracle
  or.roc.list = list()
  or.TPR.mat = matrix(NA, nlam, numSims)
  or.S.hat.mat = matrix(NA, nlam, numSims)
  or.TP.mat = matrix(NA, nlam, numSims)
  or.prec.mat = matrix(NA, nlam, numSims)
}
if(has.coat){
  # coat
  coat.roc.list = list()
  coat.TPR.mat = matrix(NA, nlam, numSims)
  coat.S.hat.mat = matrix(NA, nlam, numSims)
  coat.TP.mat = matrix(NA, nlam, numSims)
  coat.prec.mat = matrix(NA, nlam, numSims)
}
if(has.propr){
  # coat
  pr.roc.list = list()
  pr.TPR.mat = matrix(NA, nlam, numSims)
  pr.S.hat.mat = matrix(NA, nlam, numSims)
  pr.TP.mat = matrix(NA, nlam, numSims)
  pr.prec.mat = matrix(NA, nlam, numSims)
}
for(i in 1:numSims){
  # cl
  cl.sim.tmp = readRDS(paste0(
    output_dir, "/roccurves", "/classo_roc", file.end, "_sim", i, ".rds"
  ))
  cl.roc.list[[i]] = cl.sim.tmp
  cl.TPR.mat[, i] = cl.sim.tmp["tpr", ]
  cl.S.hat.mat[, i] = cl.sim.tmp["S_hat", ]
  cl.TP.mat[, i] = cl.sim.tmp["TP", ]
  cl.prec.mat[, i] = cl.sim.tmp["precision", ]
  # slr
  slr.sim.tmp = readRDS(paste0(
    output_dir, "/roccurves", "/slr_roc", file.end, "_sim", i, ".rds"
  ))
  slr.roc.list[[i]] = slr.sim.tmp
  slr.TPR.mat[, i] = slr.sim.tmp["tpr", ]
  slr.S.hat.mat[, i] = slr.sim.tmp["S_hat", ]
  slr.TP.mat[, i] = slr.sim.tmp["TP", ]
  slr.prec.mat[, i] = slr.sim.tmp["precision", ]
  if(has.oracle){
    # oracle
    or.sim.tmp = readRDS(paste0(
      output_dir, "/roccurves", "/oracle_roc", file.end, "_sim", i, ".rds"
    ))
    or.roc.list[[i]] = or.sim.tmp
    or.TPR.mat[, i] = or.sim.tmp["tpr", ]
    or.S.hat.mat[, i] = or.sim.tmp["S_hat", ]
    or.TP.mat[, i] = or.sim.tmp["TP", ]
    or.prec.mat[, i] = or.sim.tmp["precision", ]
  }
  if(has.coat){
    # coat
    coat.sim.tmp = readRDS(paste0(
      output_dir, "/roccurves", "/coat_roc", file.end, "_sim", i, ".rds"
    ))
    coat.roc.list[[i]] = coat.sim.tmp
    coat.TPR.mat[, i] = coat.sim.tmp["tpr", ]
    coat.S.hat.mat[, i] = coat.sim.tmp["S_hat", ]
    coat.TP.mat[, i] = coat.sim.tmp["TP", ]
    coat.prec.mat[, i] = coat.sim.tmp["precision", ]
  }
  if(has.propr){
    # oracle
    pr.sim.tmp = readRDS(paste0(
      output_dir, "/roccurves", "/propr_roc", file.end, "_sim", i, ".rds"
    ))
    pr.roc.list[[i]] = pr.sim.tmp
    pr.TPR.mat[, i] = pr.sim.tmp["tpr", ]
    pr.S.hat.mat[, i] = pr.sim.tmp["S_hat", ]
    pr.TP.mat[, i] = pr.sim.tmp["TP", ]
    pr.prec.mat[, i] = pr.sim.tmp["precision", ]
  }
}

# plot one sim
idx = 1
ggdat0 = data.frame(
  Method = c(
    rep("classo", nlam), rep("slr", nlam), rep("oracle", nlam),
    rep("propr", nlam)
  ),
  S_hat = c(
    cl.S.hat.mat[, idx], slr.S.hat.mat[, idx], or.S.hat.mat[, idx],
    pr.S.hat.mat[, idx]),
  TPR = c(
    cl.TPR.mat[, idx], slr.TPR.mat[, idx], or.TPR.mat[, idx],
    pr.TPR.mat[, idx])
)
ggplot(ggdat0, aes(x = S_hat, y = TPR, color = Method)) +
  geom_line()

# average over each possible S.hat/TP value
# stack columns so which() is more interpretable
cl.TPR.vec = as.vector(cl.TPR.mat)
cl.S.hat.vec = as.vector(cl.S.hat.mat)
cl.TP.vec = as.vector(cl.TP.mat)
cl.prec.vec = as.vector(cl.prec.mat)
slr.TPR.vec = as.vector(slr.TPR.mat)
slr.S.hat.vec = as.vector(slr.S.hat.mat)
slr.TP.vec = as.vector(slr.TP.mat)
slr.prec.vec = as.vector(slr.prec.mat)
if(has.oracle){
  or.TPR.vec = as.vector(or.TPR.mat)
  or.S.hat.vec = as.vector(or.S.hat.mat)
  or.TP.vec = as.vector(or.TP.mat)
  or.prec.vec = as.vector(or.prec.mat)
}
if(has.coat){
  coat.TPR.vec = as.vector(coat.TPR.mat)
  coat.S.hat.vec = as.vector(coat.S.hat.mat)
  coat.TP.vec = as.vector(coat.TP.mat)
  coat.prec.vec = as.vector(coat.prec.mat)
}
if(has.propr){
  pr.TPR.vec = as.vector(pr.TPR.mat)
  pr.S.hat.vec = as.vector(pr.S.hat.mat)
  pr.TP.vec = as.vector(pr.TP.mat)
  pr.prec.vec = as.vector(pr.prec.mat)
}

# get the averages
S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec)))
if(has.oracle & !has.coat & !has.propr){
  S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec)))
} else if(has.oracle & !has.coat & has.propr){
  S.hat.vals = sort(unique(c(
    cl.S.hat.vec, slr.S.hat.vec, or.S.hat.vec, pr.S.hat.vec)))
}
cl.TPR.avg = rep(NA, length(S.hat.vals))
cl.TP.avg = rep(NA, length(S.hat.vals))
cl.prec.avg = rep(NA, length(S.hat.vals))
slr.TPR.avg = rep(NA, length(S.hat.vals))
slr.TP.avg = rep(NA, length(S.hat.vals))
slr.prec.avg = rep(NA, length(S.hat.vals))
if(has.oracle){
  or.TPR.avg = rep(NA, length(S.hat.vals))
  or.TP.avg = rep(NA, length(S.hat.vals))
  or.prec.avg = rep(NA, length(S.hat.vals))
}
if(has.coat){
  coat.TPR.avg = rep(NA, length(S.hat.vals))
  coat.TP.avg = rep(NA, length(S.hat.vals))
  coat.prec.avg = rep(NA, length(S.hat.vals))
}
if(has.propr){
  pr.TPR.avg = rep(NA, length(S.hat.vals))
  pr.TP.avg = rep(NA, length(S.hat.vals))
  pr.prec.avg = rep(NA, length(S.hat.vals))
}
for(i in 1:length(S.hat.vals)){
  val.tmp = S.hat.vals[i]
  # classo
  cl.which.idx.tmp = which(cl.S.hat.vec == val.tmp)
  cl.TPR.avg[i] = mean(cl.TPR.vec[cl.which.idx.tmp])
  cl.TP.avg[i] = mean(cl.TP.vec[cl.which.idx.tmp])
  cl.prec.avg[i] = mean(cl.prec.vec[cl.which.idx.tmp])
  # slr
  slr.which.idx.tmp = which(slr.S.hat.vec == val.tmp)
  slr.TPR.avg[i] = mean(slr.TPR.vec[slr.which.idx.tmp])
  slr.TP.avg[i] = mean(slr.TP.vec[slr.which.idx.tmp])
  slr.prec.avg[i] = mean(slr.prec.vec[slr.which.idx.tmp])
  if(has.oracle){
    # oracle
    or.which.idx.tmp = which(or.S.hat.vec == val.tmp)
    or.TPR.avg[i] = mean(or.TPR.vec[or.which.idx.tmp])
    or.TP.avg[i] = mean(or.TP.vec[or.which.idx.tmp])
    or.prec.avg[i] = mean(or.prec.vec[or.which.idx.tmp])
  }
  if(has.coat){
    # coat
    coat.which.idx.tmp = which(coat.S.hat.vec == val.tmp)
    coat.TPR.avg[i] = mean(coat.TPR.vec[coat.which.idx.tmp])
    coat.TP.avg[i] = mean(coat.TP.vec[coat.which.idx.tmp])
    coat.prec.avg[i] = mean(coat.prec.vec[coat.which.idx.tmp])
  }
  if(has.propr){
    # propr
    pr.which.idx.tmp = which(pr.S.hat.vec == val.tmp)
    pr.TPR.avg[i] = mean(pr.TPR.vec[pr.which.idx.tmp])
    pr.TP.avg[i] = mean(pr.TP.vec[pr.which.idx.tmp])
    pr.prec.avg[i] = mean(pr.prec.vec[pr.which.idx.tmp])
  }
}

# plot
data.gg = rbind(
  data.frame(
    S_hat = S.hat.vals, TPR = cl.TPR.avg, TP = cl.TP.avg, 
    precision = cl.prec.avg, Method = "classo"), 
  data.frame(
    S_hat = S.hat.vals, TPR = slr.TPR.avg, TP = slr.TP.avg, 
    precision = slr.prec.avg, Method = "slr")
)
if(has.oracle){
  data.gg = rbind(
    data.gg, 
    data.frame(
      S_hat = S.hat.vals, TPR = or.TPR.avg, TP = or.TP.avg, 
      precision = or.prec.avg, Method = "oracle"))
}
if(has.coat){
  data.gg = rbind(
    data.gg, 
    data.frame(
      S_hat = S.hat.vals, TPR = coat.TPR.avg, TP = coat.TP.avg, 
      precision = coat.prec.avg, Method = "coat"))
}
if(has.propr){
  data.gg = rbind(
    data.gg, 
    data.frame(
      S_hat = S.hat.vals, TPR = pr.TPR.avg, TP = pr.TP.avg, 
      precision = pr.prec.avg, Method = "propr"))
}
data.gg$Method = factor(data.gg$Method, levels = levels.gg)
tp_roc = ggplot(
  data.gg[!is.na(data.gg$TP),], 
  aes(x = S_hat, y = TP, color = Method, linetype = Method)) + 
  geom_line(alpha = 0.75, size = 1, na.rm = TRUE) +
  # geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
tpr_roc = ggplot(
  data.gg[!is.na(data.gg$TPR),], 
  aes(x = S_hat, y = TPR, color = Method, linetype = Method)) + 
  geom_line(alpha = 0.75, size = 1, na.rm = TRUE) +
  # geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
prec_roc = ggplot(
  data.gg[!is.na(data.gg$precision),], 
  aes(x = S_hat, y = precision, color = Method, linetype = Method)) + 
  geom_line(alpha = 0.75, size = 1, na.rm = TRUE) +
  # geom_point(alpha = 0.5, na.rm = TRUE) +
  theme_bw()
tpr_roc
# ggarrange(tp_roc, tpr_roc, prec_roc, nrow = 1)

ggsave(
  filename = paste0(
    "20211202_", 
    sigma.settings, "_noise", sigma_eps, 
    "_", theta.settings, 
    "_val", values.theta,
    "_rocs.pdf"),
  plot = last_plot(),
  width = 6, height = 5, units = c("in")
)


