# Data for predicting BMI 
raw_data <- read.csv("Data/BMI_and_Counts.txt",sep=" ")
y <- raw_data$BMI
X <- raw_data[,-c(1,2)]
X <- apply(X,2,as.integer)
prop_zero_per_sample <- apply(X, 1, function(a) sum(a==0)/ncol(X))
prop_zero_per_otu <- apply(X, 2, function(a) sum(a==0)/nrow(X))

## The covariate matrix is sparse. Because the logarithm of zero does not exist, we replace all zeros with a pseudocount of 0.5
## before transforming X into compositional data. 
X[which(X==0)] <- 0.5
X.prop <- sweep(X,MARGIN = 1, STATS = rowSums(X), FUN = "/")
save(y,X,X.prop,raw_data,file='Data/BMI.rda')
