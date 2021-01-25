output_home = "Kristyn/Experiments/output/"
rng.seed = 123
n = 100
p = 1000
rho = 0.5 # 0.2, 0.5
get_lambda = "original"
################################################################################
# Compositional Lasso #
################################################################################

complasso = readRDS(paste0(
  "Kristyn/Experiments/output",
  "/complasso_simulations_summary", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_seed", rng.seed,
  ".rds"
))
print(complasso[, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################
slr = readRDS(paste0(
  output_home, 
  "/slr_prediction", 
  "_", get_lambda, 
  "_K", cv.K, 
  "_seed", rng.seed,
  ".rds"
))
dim(slr)
slr.mse = colMeans(slr)
print(paste0("mean prediction error: ", mean(slr.mse)))
print(paste0("standard deviation: ", (sd(slr.mse))))
print(paste0("standard error: ", (sd(slr.mse)) / sqrt(rep.n)))
print(paste0(
  "95% CI: (", 
  mean(slr.mse) - 2 * (sd(slr.mse)) / sqrt(rep.n), 
  ", ",
  mean(slr.mse) + 2 * (sd(slr.mse)) / sqrt(rep.n), ")"
))
