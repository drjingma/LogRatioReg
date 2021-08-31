# Method: Simulation study for compositional Lasso
# Purpose: Simulate data, fit supervised log-ratios method to the data
# Date: 08/30/2021

# directory
# aim31: 
output_dir = "Kristyn/Experiments/complasso_simulations/output_aim31_metrics_roc"
# kristyn: 
# output_dir = "Kristyn/Experiments/complasso_simulations/output_kristyn_metrics_roc"

################################################################################
# libraries and settings

library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

numSims = 100
rng.seed = 123

# Settings to toggle with
rho.type = "square" # 1 = "absolute value", 2 = "square"
beta.settings = "new"
linkage = "average"
tol = 1e-4
nlam = 100
intercept = TRUE
K = 10
n = 100
p = 200
rho = 0.5 # 0.2, 0.5
scaling = TRUE

file.end = paste0(
  "_dim", n, "x", p, 
  "_", beta.settings, 
  "_rho", rho, 
  "_int", intercept,
  "_scale", scaling,
  "_K", K,
  "_seed", rng.seed,
  ".rds")

################################################################################
# plot metrics
metrics_names = c("PEtr", "PEte", "EA1", "EA2", "EAInfty", 
            "FP", "FN", "TPR", "betaSparsity")


# import metrics
cl.sims.list = list()
slr.sims.list = list()
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
}
cl.sims = as.data.frame(rbindlist(cl.sims.list))
slr.sims = as.data.frame(rbindlist(slr.sims.list))

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

# boxplots for the slr and classo metrics
cl.sims.gg = reshape2::melt(cl.sims)
cl.sims.gg$type = "classo"
slr.sims.gg = reshape2::melt(slr.sims)
slr.sims.gg$type = "slr"
data.gg = rbind(cl.sims.gg, slr.sims.gg)
data.gg = dplyr::filter(data.gg, variable %in% metrics_names)
data.gg$type = factor(data.gg$type, levels = c("classo", "slr"))
ggplot(data.gg, aes(x = type, y = value, color = type)) + 
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
ggarrange(plt.PEtr, plt.PEte)



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
}

# average over each possible S.hat value
# stack columns so which() is more interpretable
cl.TPR.vec = as.vector(cl.TPR.mat)
cl.S.hat.vec = as.vector(cl.S.hat.mat)
slr.TPR.vec = as.vector(slr.TPR.mat)
slr.S.hat.vec = as.vector(slr.S.hat.mat)
# get the averages
S.hat.vals = sort(unique(c(cl.S.hat.vec, slr.S.hat.vec)))
cl.tpr.avg = rep(NA, length(S.hat.vals))
slr.tpr.avg = rep(NA, length(S.hat.vals))
for(i in 1:length(S.hat.vals)){
  val.tmp = S.hat.vals[i]
  # classo
  cl.which.idx.tmp = which(cl.S.hat.vec == val.tmp)
  cl.tpr.avg[i] = mean(cl.TPR.vec[cl.which.idx.tmp])
  # slr
  slr.which.idx.tmp = which(slr.S.hat.vec == val.tmp)
  slr.tpr.avg[i] = mean(slr.TPR.vec[slr.which.idx.tmp])
}

# plot
plot(S.hat.vals, cl.tpr.avg, type = "l", col = 2)
lines(S.hat.vals, slr.tpr.avg, col = 3)



################################################################################
# plot roc curves (same lambda sequence)

# import roc curves and organize TPR, S.hat, lambda information
# cl
cl.roc.list2 = list()
# each row corresponds to a lambda in the lambda sequence (different in ea. sim)
# each column corresponds to a different simulation
cl.TPR.mat2 = matrix(NA, nlam, numSims) 
cl.S.hat.mat2 = matrix(NA, nlam, numSims)
# slr
slr.roc.list2 = list()
slr.TPR.mat2 = matrix(NA, nlam, numSims)
slr.S.hat.mat2 = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  # cl
  cl.sim.tmp = readRDS(paste0(
    output_dir, "/classo_roc_samelam", i, file.end
  ))
  cl.roc.list2[[i]] = cl.sim.tmp
  cl.TPR.mat2[, i] = cl.sim.tmp["tpr", ]
  cl.S.hat.mat2[, i] = cl.sim.tmp["S_hat", ]
  # slr
  slr.sim.tmp = readRDS(paste0(
    output_dir, "/slr_roc", i, file.end
  ))
  slr.roc.list2[[i]] = slr.sim.tmp
  slr.TPR.mat2[, i] = slr.sim.tmp["tpr", ]
  slr.S.hat.mat2[, i] = slr.sim.tmp["S_hat", ]
}


# average across different lambdas (same lambda seq in each sim, so it's okay)
# classo
cl.S.hat.avg2 = apply(cl.S.hat.mat2, 1, mean, na.rm = TRUE)
cl.TPR.avg2 = apply(cl.TPR.mat2, 1, mean, na.rm = TRUE)
# slr
slr.S.hat.avg2 = apply(slr.S.hat.mat2, 1, mean, na.rm = TRUE)
slr.TPR.avg2 = apply(slr.TPR.mat2, 1, mean, na.rm = TRUE)

# plot
plot(cl.S.hat.avg2, cl.TPR.avg2, type = "l", col = 2)
lines(slr.S.hat.avg2, slr.TPR.avg2, col = 3)

# # complasso stuff
# cl.gg.complete = data.frame(
#   "S.hat" = S.hat.cl.avg,
#   "TPR" = TPR.cl.avg)
# cl.gg.complete$Type = "classo"
# # slr stuff
# slr.gg.complete = data.frame(
#   "S.hat" = S.hat.slr.avg,
#   "TPR" = TPR.slr.avg)
# slr.gg.complete$Type = "slr"
# # slr0.5 stuff
# slr0.5.gg.complete = data.frame(
#   "S.hat" = S.hat.slr0.5.avg,
#   "TPR" = TPR.slr0.5.avg)
# slr0.5.gg.complete$Type = "slr0.5"
# # slr1 stuff
# slr1.gg.complete = data.frame(
#   "S.hat" = S.hat.slr1.avg,
#   "TPR" = TPR.slr1.avg)
# slr1.gg.complete$Type = "slr1"
# # ggplot
# gg.complete = rbind(
#   cl.gg.complete, slr.gg.complete, slr0.5.gg.complete, slr1.gg.complete)
# gg.complete$Type = factor(gg.complete$Type, 
#                           levels = c("classo", "slr", "slr0.5", "slr1"))
# ggplot(gg.complete, aes(x = S.hat, y = TPR, color = Type
#                         , shape = Type, linetype = Type
# )) +
#   geom_line(size = 1) +
#   geom_point(size = 3) +
#   xlim(0, 40) +
#   theme_bw() + 
#   theme(text = element_text(size = 20))
# 
# all.equal(slr.gg.complete, slr.gg.complete.old)

