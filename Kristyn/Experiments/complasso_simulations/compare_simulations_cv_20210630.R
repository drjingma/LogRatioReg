library(ggplot2)
library(reshape2)
library(ggpubr)

output_dir = "Kristyn/Experiments/complasso_simulations/output_20210630"

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

################################################################################
# Compositional Lasso #
################################################################################
file.end = paste0(
  "_dim", n, "x", p, 
  "_", beta.settings, 
  "_rho", rho, 
  "_int", intercept,
  "_scale", scaling,
  "_K", K,
  "_seed", rng.seed,
  ".rds")

complasso.sim1 = readRDS(paste0(
  output_dir, "/complasso_sim1", file.end
))
complasso.summaries = readRDS(paste0(
  output_dir, "/complasso_summaries", file.end
))
print(complasso.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################
file.end = paste0(
  "_dim", n, "x", p, 
  "_", beta.settings, 
  "_rho", rho, 
  "_type", rho.type,
  "_int", intercept,
  "_K", K,
  "_seed", rng.seed,
  ".rds")
slr.sims = readRDS(paste0(
  output_dir, "/slr_sims", file.end
))
slr.summaries = readRDS(paste0(
  output_dir, "/slr_summaries", file.end
))
print(slr.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios with alpha -- alpha = 0.5 #
################################################################################
alpha = 0.5
file.end = paste0(
  "_dim", n, "x", p, 
  "_", beta.settings, 
  "_alpha", alpha, 
  "_rho", rho, 
  "_type", rho.type,
  "_int", intercept,
  "_K", K,
  "_seed", rng.seed,
  ".rds")
slr0.5.sims = readRDS(paste0(
  output_dir, "/slralpha_sims", file.end
))
slr0.5.summaries = readRDS(paste0(
  output_dir, "/slralpha_summaries", file.end
))
print(slr0.5.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios with alpha -- alpha = 1 #
################################################################################
alpha = 1
file.end = paste0(
  "_dim", n, "x", p, 
  "_", beta.settings, 
  "_alpha", alpha, 
  "_rho", rho, 
  "_type", rho.type,
  "_int", intercept,
  "_K", K,
  "_seed", rng.seed,
  ".rds")
slr1.sims = readRDS(paste0(
  output_dir, "/slralpha_sims", file.end
))
slr1.summaries = readRDS(paste0(
  output_dir, "/slralpha_summaries", file.end
))
print(slr1.summaries[metrics, c("mean", "se")])


################################################################################
################################################################################
################################################################################


# plot

cl.sims.gg = melt(data.frame(t(complasso.sims)))
cl.sims.gg$type = "CompLasso"
slr.sims.gg = melt(data.frame(t(slr.sims)))
slr.sims.gg$type = "SLR"
slr1.sims.gg = melt(data.frame(t(slr1.sims)))
slr1.sims.gg$type = "SLRalpha1"
slr0.5.sims.gg = melt(data.frame(t(slr0.5.sims)))
slr0.5.sims.gg$type = "SLRalpha0.5"
data.gg = rbind(cl.sims.gg, slr.sims.gg, slr1.sims.gg, slr0.5.sims.gg)
data.gg = dplyr::filter(data.gg, variable %in% metrics)
data.gg$type = factor(data.gg$type, levels = c("CompLasso", "SLR", "SLRalpha1", "SLRalpha0.5"))
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

