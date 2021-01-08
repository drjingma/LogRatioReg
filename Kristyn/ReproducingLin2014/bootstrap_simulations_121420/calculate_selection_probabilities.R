getwd()

################################################################################
# read in output from bootstrap_parallelized_121420.R
################################################################################
res0 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/bootstraps_121420_safekeep.rds"))
dim(res0$selected_variables)
all(!is.na(res0$selected_variables)) # no NAs
sort(res0$selection_percentages)
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megamonas 
# 50 
# Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Oscillibacter 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Zymophilus 
# 53 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella 
# 78 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium 
# 80 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus 
# 85 

################################################################################
# read in output from bootstrap_parallelized_121420_manual.R
################################################################################
res1 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/bootstraps_121420_safekeep.rds"))
dim(res1$selected_variables)
all(!is.na(res1$selected_variables)) # no NAs
sort(res1$selection_percentages)
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megamonas 
# 50 
# Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Oscillibacter 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Zymophilus 
# 53 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella 
# 78 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium 
# 80 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus 
# 85 
# they match the results above... but why?

################################################################################
# read in output from bootstrap_parallelized_121420.R
################################################################################
res0 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/bootstraps_121420_safekeep.rds"))
dim(res0$selected_variables)
all(!is.na(res0$selected_variables)) # no NAs
sort(res0$selection_percentages)
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megamonas 
# 50 
# Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Oscillibacter 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Zymophilus 
# 53 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella 
# 78 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium 
# 80 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus 
# 85 

################################################################################
# read in output from bootstrap_notparallelized_121420_manual.R
################################################################################
res2 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/bootstraps_121420_notparallel_manual.rds"))
dim(res2$selected_variables)
all(!is.na(res2$selected_variables)) # no NAs
sort(res2$selection_percentages)
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megamonas 
# 52 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea 
# 55 
# Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes 
# 57 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella 
# 81 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus 
# 83 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium 
# 90 

################################################################################
# read in output from bootstrap_notparallelized_121420.R
################################################################################
res3 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/bootstraps_121420_notparallel.rds"))
dim(res3$selected_variables)
all(!is.na(res3$selected_variables)) # no NAs
sort(res3$selection_percentages)
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megamonas 
# 53 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea 
# 56 
# Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes 
# 58 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella 
# 82 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus 
# 84 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium 
# 91

################################################################################
# read in output from using_subcompositions/bootstrap_parallelized_121420.R
################################################################################
res00 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420/using_subcompositions",
  "/bootstraps_sub_121420.rds"))
dim(res00$selected_variables)
all(!is.na(res00$selected_variables)) # no NAs
sort(res00$selection_percentages)
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megamonas 
# 50 
# Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Oscillibacter 
# 51 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Zymophilus 
# 53 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella 
# 78 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium 
# 80 
# Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus 
# 85