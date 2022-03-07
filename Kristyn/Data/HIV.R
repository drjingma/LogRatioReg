# Date: 3/7/2022

# HIV: a data set in selbal package
#   n = 155 samples, 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 covariate (MSM), 
#   1 response (HIV_Status - binary)
X_hiv1 = selbal::HIV[, 1:60]


# sCD14: a data set in selbal package
#   n = 151 samples (a subset from HIV), 
#   p = 60 taxa (counts for microbial taxa at genus level), 
#   1 response (sCD14 - continuous)
X_hiv2 = selbal::sCD14[, 1:60]
