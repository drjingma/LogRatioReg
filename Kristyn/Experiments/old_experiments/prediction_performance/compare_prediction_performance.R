output_home = "Kristyn/Experiments/prediction_performance/output/"
rng.seed = 123
cv.K = 5
intercept = FALSE
# for some reason, results for slr when intercept = TRUE are very bad, 
#   whereas when intercept = FALSE not that different compared to complasso
################################################################################
# Compositional Lasso #
################################################################################

complasso = readRDS(paste0(
  output_home, 
  "complasso_prediction", 
  "_int", intercept,
  "_K", cv.K, 
  "_seed", rng.seed,
  ".rds"
))
complasso.mse = colMeans(complasso)
complasso.mse.mean = mean(complasso.mse)
complasso.mse.sd = sd(complasso.mse)
complasso.mse.se = complasso.mse.sd / sqrt(ncol(complasso))
complasso.summary = data.frame(
  mean = complasso.mse.mean, sd = complasso.mse.sd, 
  se = complasso.mse.se, 
  lower = complasso.mse.mean - 2 * complasso.mse.se, 
  upper = complasso.mse.mean + 2 * complasso.mse.se)

################################################################################
# Supervised Log Ratios #
################################################################################
slr = readRDS(paste0(
  output_home, 
  "/slr_prediction", 
  "_int", intercept,
  "_K", cv.K, 
  "_seed", rng.seed,
  ".rds"
))
slr.mse = colMeans(slr)
slr.mse.mean = mean(slr.mse)
slr.mse.sd = sd(slr.mse)
slr.mse.se = slr.mse.sd / sqrt(ncol(slr))
slr.summary = data.frame(
  mean = slr.mse.mean, sd = slr.mse.sd, 
  se = slr.mse.se, 
  lower = slr.mse.mean - 2 * slr.mse.se, 
  upper = slr.mse.mean + 2 * slr.mse.se)

################################################################################
# Compare #
################################################################################
complasso.summary
slr.summary
