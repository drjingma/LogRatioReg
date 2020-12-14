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

