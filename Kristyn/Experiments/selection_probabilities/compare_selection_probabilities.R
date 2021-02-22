output_home = "Kristyn/Experiments/selection_probabilities/output/"
################################################################################
# Compositional Lasso #
################################################################################

complasso = readRDS(paste0(output_home, "complasso_selection_refit_seed123.rds"))
length(complasso$selection_percentages)
dim(complasso$selected_variables)
all(!is.na(complasso$selected_variables)) # no NAs
sort(complasso$selection_percentages)

################################################################################
# Supervised Log Ratios - refit Log Ratios #
################################################################################

slrLR = readRDS(paste0(output_home, "slr_selection_refitLRs_seed123.rds"))
length(slrLR$selection_percentages)
dim(slrLR$selected_variables)
all(!is.na(slrLR$selected_variables)) # no NAs
sort(slrLR$selection_percentages)

################################################################################
# Supervised Log Ratios - refit Log Contrasts #
################################################################################

slrLC = readRDS(paste0(output_home, "slr_selection_refitLCs_seed123.rds"))
length(slrLC$selection_percentages)
dim(slrLC$selected_variables)
all(!is.na(slrLC$selected_variables)) # no NAs
sort(slrLC$selection_percentages)

