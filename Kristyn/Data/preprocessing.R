################################################################################
# Shi, et al. 2016 data -- took out counts that were too sparse 
#   (otu.prop too small), acc. to paper
cleaned_counts <- read.csv("Kristyn/Data/cleaned_BMI_Counts.csv")
View(cleaned_counts)
dim(cleaned_counts)
names(cleaned_counts)

# compare to previous dataset
load(paste0("Data/", "BMI.rda"))
dim(X)

# full data set?
combo <- read.csv("Kristyn/Data/combo_ffq_adj.txt", sep = " ")
View(combo)
dim(combo)
names(combo)

