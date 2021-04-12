output_dir = "Kristyn/Experiments/complasso_simulations/output"
rng.seed = 123
n = 50
p = 30
rho = 0.2 # 0.2, 0.5
intercept = FALSE

################################################################################
# Compositional Lasso #
################################################################################

complasso = readRDS(paste0(
  output_dir,
  "/complasso_gic_simulations", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
print(complasso[, c("mean", "se")])

################################################################################
# Supervised Log Ratios - GIC on Log-Ratios #
################################################################################
slr.lr = readRDS(paste0(
  output_dir,
  "/slr_gicLR_simulations", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
print(slr.lr[, c("mean", "se")])

################################################################################
# Supervised Log Ratios - GIC on Log-Contrasts #
################################################################################
slr.lc = readRDS(paste0(
  output_dir,
  "/slr_gicLC_simulations", 
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept, 
  "_seed", rng.seed,
  ".rds"
))
print(slr.lc[, c("mean", "se")])
