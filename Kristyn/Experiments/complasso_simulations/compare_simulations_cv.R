library(ggplot2)
library(reshape2)

output_dir = "Kristyn/Experiments/complasso_simulations/output"
rng.seed = 123
n = 50
p = 30
rho = 0.2 # 0.2, 0.5
intercept = TRUE

################################################################################
# Compositional Lasso #
################################################################################
complasso.sims = readRDS(paste0(
  output_dir,
  "/complasso_cv_sims", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
complasso.summaries = readRDS(paste0(
  output_dir,
  "/complasso_cv_summaries", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
print(complasso.summaries[, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################
slr.sims = readRDS(paste0(
  output_dir,
  "/slr_cv_sims", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
slr.summaries = readRDS(paste0(
  output_dir,
  "/slr_cv_summaries", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
print(slr.summaries[, c("mean", "se")])




# plot

cl.sims.gg = melt(data.frame(t(complasso.sims)))
cl.sims.gg$type = "CompLasso"
slr.sims.gg = melt(data.frame(t(slr.sims)))
slr.sims.gg$type = "SLR"
data.gg = rbind(cl.sims.gg, slr.sims.gg)
data.gg$type = factor(data.gg$type, levels = c("CompLasso", "SLR"))
ggplot(data.gg, aes(x = type, y = value, color = type)) + 
  facet_wrap(vars(variable), scales = "free_y") + geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())
