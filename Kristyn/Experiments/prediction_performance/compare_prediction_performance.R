output_home = "Kristyn/Experiments/output/"
rng.seed = 123
rep.n = 100
cv.K = 3
get_lambda = "original"
################################################################################
# Compositional Lasso #
################################################################################

complasso = readRDS(paste0(
  output_home, 
  "complasso_prediction", 
  "_K", cv.K, 
  "_seed", rng.seed,
  ".rds"
))
dim(complasso)
complasso.mse = colMeans(complasso)
print(paste0("mean prediction error: ", mean(complasso.mse)))
print(paste0("standard deviation: ", (sd(complasso.mse))))
print(paste0("standard error: ", (sd(complasso.mse)) / sqrt(rep.n)))
print(paste0(
  "95% CI: (", 
  mean(complasso.mse) - 2 * (sd(complasso.mse)) / sqrt(rep.n), 
  ", ",
  mean(complasso.mse) + 2 * (sd(complasso.mse)) / sqrt(rep.n), ")"
  ))

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
