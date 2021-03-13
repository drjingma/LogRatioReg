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

complasso = readRDS(paste0(
  output_dir,
  "/complasso_cv_simulations", 
  "_PBA", 
  "_theta", generate.theta,
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(complasso[, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################
slr.lr = readRDS(paste0(
  output_dir,
  "/slr_cv_simulations", 
  "_PBA", 
  "_theta", generate.theta,
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_seed", rng.seed,
  ".rds"))
print(slr.lr[, c("mean", "se")])

