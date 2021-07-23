library(ggplot2)
library(ggpubr)
library(data.table)
library(reshape2)

# output_dir = "Kristyn/Experiments/complasso_simulations/output_20210630"
output_dir = "Kristyn/Experiments/complasso_simulations/output_20210722"

numSims = 100
rng.seed = 123
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

# other stuff
metrics = c("PEtr", "PEte", "EA1", "EA2", "EAInfty", 
            "FP", "FN", "TPR", "betaSparsity")

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
# Compositional Lasso #
################################################################################

# import
cl.sims.list = list()
for(i in 1:numSims){
  sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/complasso_sim", i, file.end
  ))))
  rownames(sim.tmp) = NULL
  cl.sims.list[[i]] = data.table(sim.tmp)
}
cl.sims = as.data.frame(rbindlist(cl.sims.list))

# summary stats
cl.eval.means = apply(cl.sims, 2, mean)
cl.eval.sds = apply(cl.sims, 2, sd)
cl.eval.ses = cl.eval.sds / sqrt(numSims)
cl.summaries = data.frame(
  "mean" = cl.eval.means, "sd" = cl.eval.sds, "se" = cl.eval.ses)
cl.summaries

print(cl.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################

# import
slr.sims.list = list()
for(i in 1:numSims){
  sim.tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_sim", i, file.end
  ))))
  rownames(sim.tmp) = NULL
  slr.sims.list[[i]] = data.table(sim.tmp)
}
slr.sims = as.data.frame(rbindlist(slr.sims.list))

# summary stats
slr.eval.means = apply(slr.sims, 2, mean)
slr.eval.sds = apply(slr.sims, 2, sd)
slr.eval.ses = slr.eval.sds / sqrt(numSims)
slr.summaries = data.frame(
  "mean" = slr.eval.means, "sd" = slr.eval.sds, "se" = slr.eval.ses)
slr.summaries

print(slr.summaries[metrics, c("mean", "se")])

################################################################################
################################################################################
################################################################################


# plot

cl.sims.gg = reshape2::melt(cl.sims)
cl.sims.gg$type = "CompLasso"
slr.sims.gg = reshape2::melt(slr.sims)
slr.sims.gg$type = "SLR"
data.gg = rbind(cl.sims.gg, slr.sims.gg)
data.gg = dplyr::filter(data.gg, variable %in% metrics)
data.gg$type = factor(data.gg$type, levels = c("CompLasso", "SLR"))
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

