# last updated: 04/19/2021
# compare different methods from SigmaW tree simulations

library(ggplot2)
library(reshape2)
library(ggpubr)

output_dir = "Kristyn/Experiments/SigmaW_simulations/output"
rng.seed = 123
n = 100
p = 200
rho = 0 # 0.2, 0.5
# indices.theta = 1
indices.theta = p - 1
# indices.theta = c(159, 179, 14, 195, 170)
# indices.theta = c(197, 198, 199)
intercept = TRUE

################################################################################
# Compositional Lasso #
################################################################################
complasso.sims = readRDS(paste0(
  output_dir,
  "/complasso_cv_sims", 
  "_SigmaW", 
  "_theta_", paste(indices.theta, collapse = "_"),
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
complasso.summaries = readRDS(paste0(
  output_dir,
  "/complasso_cv_summaries", 
  "_SigmaW", 
  "_theta_", paste(indices.theta, collapse = "_"),
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(complasso.summaries[, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################
slr.sims = readRDS(paste0(
  output_dir,
  "/slr_cv_sims", 
  "_SigmaW", 
  "_theta_", paste(indices.theta, collapse = "_"),
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
slr.summaries = readRDS(paste0(
  output_dir,
  "/slr_cv_summaries", 
  "_SigmaW", 
  "_theta_", paste(indices.theta, collapse = "_"),
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(slr.summaries[, c("mean", "se")])

################################################################################
# Oracle #
################################################################################
or.sims = readRDS(paste0(
  output_dir,
  "/oracle_cv_sims",
  "_SigmaW", 
  "_theta_", paste(indices.theta, collapse = "_"),
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
or.summaries = readRDS(paste0(
  output_dir,
  "/oracle_cv_summaries", 
  "_SigmaW", 
  "_theta_", paste(indices.theta, collapse = "_"),
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(or.summaries[, c("mean", "se")])

################################################################################
# Selbal #
################################################################################
# selbal.sims = readRDS(paste0(
#   output_dir,
#   "/selbal_cv_sims", 
#   "_PBA", 
#   "_theta_", paste(indices.theta, collapse = "_"),
#   "_dim", n, "x", p, 
#   "_rho", rho, 
#   "_int", intercept,
#   "_seed", rng.seed,
#   ".rds"))
# selbal.summaries = readRDS(paste0(
#   output_dir,
#   "/selbal_cv_summaries", 
#   "_PBA", 
#   "_theta_", paste(indices.theta, collapse = "_"),
#   "_dim", n, "x", p, 
#   "_rho", rho, 
#   "_int", intercept,
#   "_seed", rng.seed,
#   ".rds"))
# print(selbal.summaries[, c("mean", "se")])


# plot

cl.sims.gg = melt(data.frame(t(complasso.sims)))
cl.sims.gg$type = "CompLasso"
slr.sims.gg = melt(data.frame(t(slr.sims)))
slr.sims.gg$type = "SLR"
or.sims.gg = melt(data.frame(t(or.sims)))
or.sims.gg$type = "Oracle"
# sel.sims.gg = melt(data.frame(t(selbal.sims)))
# sel.sims.gg$type = "selbal"
# data.gg = rbind(cl.sims.gg, slr.sims.gg, or.sims.gg, sel.sims.gg)
# data.gg$type = factor(data.gg$type, levels = c("CompLasso", "SLR", "Oracle", "selbal"))

data.gg = rbind(cl.sims.gg, slr.sims.gg, or.sims.gg)
data.gg$type = factor(data.gg$type, levels = c("CompLasso", "SLR", "Oracle"))

ggplot(data.gg, aes(x = type, y = value, color = type)) + 
  facet_wrap(vars(variable), scales = "free_y") + geom_boxplot() + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# zoom in to PEtr, PEte
PEtr.gg = data.gg[data.gg$variable == "PEtr", ]
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
# plt.PEtr
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
# plt.PEte
ggarrange(plt.PEtr, plt.PEte)
