library(ggplot2)
library(reshape2)

output_dir = "Kristyn/Experiments/slr_simulations/output"
rng.seed = 123
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
generate.theta = 2
intercept = TRUE

################################################################################
# Compositional Lasso #
################################################################################
complasso.sims = readRDS(paste0(
  output_dir,
  "/complasso_cv_sims", 
  "_PBA", 
  "_theta", generate.theta,
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
complasso.summaries = readRDS(paste0(
  output_dir,
  "/complasso_cv_summaries", 
  "_PBA", 
  "_theta", generate.theta,
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(complasso.summaries[, c("mean", "se")])

cl.sims.gg = melt(data.frame(t(complasso.sims)))
ggplot(cl.sims.gg) + 
  geom_boxplot(aes(x = variable, y = value))

################################################################################
# Supervised Log Ratios #
################################################################################
slr.sims = readRDS(paste0(
  output_dir,
  "/slr_cv_sims", 
  "_PBA", 
  "_theta", generate.theta,
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
slr.summaries = readRDS(paste0(
  output_dir,
  "/slr_cv_summaries", 
  "_PBA", 
  "_theta", generate.theta,
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(slr.summaries[, c("mean", "se")])

slr.sims.gg = melt(data.frame(t(slr.sims)))
ggplot(slr.sims.gg) + 
  geom_boxplot(aes(x = variable, y = value))
