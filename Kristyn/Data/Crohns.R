# Date: 3/7/2022

source("Kristyn/Functions/supervisedlogratios.R")

# Crohn: a data set in selbal package
#   n = 975 samples, 
#   p = 48 taxa (counts for microbial taxa at genus level), 
#   1 response (y - binary)
X = selbal::Crohn[, 1:48]
y = ifelse(selbal::Crohn[, 49] == "CD", 1, 0)

# getSlrMatrix(y = y, X = X, type = "similarity")
# need to deal with 0's in X, first -- log(Xij) = Inf if Xij = 0.