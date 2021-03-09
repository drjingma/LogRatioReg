output_home = "Kristyn/Experiments/selection_probabilities/output/"
rng.seed = 123
bs.n = 500
intercept = TRUE
################################################################################
# Compositional Lasso #
################################################################################

complasso = readRDS(paste0(
  output_home,
  "complasso_selection_refit", 
  "_int", intercept, 
  "_B", bs.n, 
  "_seed", rng.seed,
  ".rds"
))
# length(complasso$selection_percentages)
# dim(complasso$selected_variables)
# all(!is.na(complasso$selected_variables)) # no NAs
head(sort(complasso$selection_percentages, decreasing = TRUE), 5)

################################################################################
# Supervised Log Ratios - refit Log Ratios #
################################################################################

slrLR = readRDS(paste0(
  output_home,
  "slr_selection", 
  "_refitLRs",
  "_int", intercept, 
  "_B", bs.n, 
  "_seed", rng.seed,
  ".rds"
))
# length(slrLR$selection_percentages)
# dim(slrLR$selected_variables)
# all(!is.na(slrLR$selected_variables)) # no NAs
head(sort(slrLR$selection_percentages, decreasing = TRUE), 5)

################################################################################
# Supervised Log Ratios - refit Log Contrasts #
################################################################################

slrLC = readRDS(paste0(
  output_home,
  "slr_selection", 
  "_refitLCs",
  "_int", intercept, 
  "_B", bs.n, 
  "_seed", rng.seed,
  ".rds"
))
# length(slrLC$selection_percentages)
# dim(slrLC$selected_variables)
# all(!is.na(slrLC$selected_variables)) # no NAs
head(sort(slrLC$selection_percentages, decreasing = TRUE), 5)

