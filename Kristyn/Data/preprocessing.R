################################################################################
# Shi, et al. 2016 data -- took out counts that were too sparse 
#   (otu.prop too small), acc. to paper
cleaned_counts <- read.csv("Kristyn/Data/cleaned_BMI_Counts.csv")
# View(cleaned_counts)
dim(cleaned_counts)
names(cleaned_counts)
# genera names
cleaned_counts_genera_names = names(cleaned_counts)[-c(1, 2)]
# subject indices
cleaned_counts_indices = cleaned_counts[, 1]
# design matrix of counts
X.shi = as.matrix(cleaned_counts[, cleaned_counts_genera_names])
X.shi[which(X.shi==0)] = 0.5
dim(X.shi)

# compare to previous dataset
raw_counts <- read.csv("Data/BMI_and_Counts.txt",sep=" ")
# subject indices
raw_counts_indices = raw_counts[, 1]
# is shi et al contained in this set?
cleaned_counts_indices %in% raw_counts_indices # yes
# get the rows in raw_counts_indices that correspond to cleaned_counts_indices
raw_rows = rep(NA, length(cleaned_counts_indices))
which_row = function(x){
  which_rows = which(raw_counts_indices == x)
  if(length(which_rows) == 1){
    return(which_rows)
  } else if(length(which_rows) == 0){
    return(NA)
  } else{
    warning("repeated indices??")
  }
}
raw_rows = sapply(cleaned_counts_indices, FUN = which_row)
# make X from previous data set match X.shi
load(paste0("Data/", "BMI.rda"))
dim(X)
X.match = X[raw_rows, cleaned_counts_genera_names]
dim(X.match)
all.equal(X.shi, X.match)

# full data set?
combo <- read.csv("Kristyn/Data/combo_ffq_adj.txt", sep = " ")
# View(combo)
dim(combo)
names(combo)
# assuming that combo is in the same order as X, get other covariates for the 
#   subset of subjects included in X.shi
combo.match = combo[raw_rows, c("calor", "tfat")]

# see why these were chosen...
cor_y = function(x){
  cor(x, y)
}
cor.y = apply(combo, MARGIN = 2, FUN = cor_y)
sort(cor.y) #...






