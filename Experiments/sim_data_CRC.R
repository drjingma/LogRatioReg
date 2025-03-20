# update 3/20/2025

set.seed(1314)

library(selbal)
library(balance)

lowCor <- F
highCor <- T

s1 <- 0.01
s2 <- 0.01
d <- 50 # can be 50 or 375
eigen.id <- 3

filename <- 'CRC'
filename <- ifelse(highCor, paste0(filename,"_highCor"), paste0(filename,"_lowCor"))
filename <- paste0(filename, '_PC', eigen.id)

file.end <- paste0(
  "/",filename,
  "_sparsity", s1*100, 
  "_ratio", max(s1,s2)/min(s1,s2), 
  '_dim', d)

load(paste0("../Data/Yachida_CRC_2019.rda"))

# zero replacement
X = as.matrix(selbal::cmultRepl2(data.list$X, zero.rep = "bayes"))

# learn the empirical correlation
X.clr <- apply(X,1,function(s) log(s) - mean(log(s)))
rho <- cor(t(X.clr))

v <- eigen(rho)$vectors[,eigen.id]

# determine active variables based on the largest positive and largest 
## negative values in the selected PC loading
if (highCor){
  # high correlation case
  sbp <- as.matrix((v <= quantile(v,s1)) - (v >= quantile(v,1-s2)))
}
if (lowCor){
  # low correlation case
  sbp <- as.matrix((v <= quantile(v,0.5+s1)) * (v >= quantile(v,0.5-s2)))
  sbp <- sbp * (2*rbinom(length(sbp),1,s1/(s1+s2)) - 1)
}
rho.sub <- rho[sbp!=0,sbp!=0]
summary(rho.sub[lower.tri(rho.sub)])
summary(abs(rho)[sbp!=0,sbp!=0])
#pheatmap::pheatmap(rho[sbp!=0,sbp!=0])
print(table(sbp))
rownames(sbp) <- colnames(data.list$X)

# partition data into training and test data
n <- nrow(X)
split.perc = 0.7
train.idx <- sample(1:n,split.perc*n) 

X.train <- X[train.idx,]
X.test <- X[-train.idx,]

# randomly add inactive variables
variables2add <- sample(which(sbp==0),d)
X.train.select <- X.train[,c(which(sbp!=0), variables2add)]
X.test.select  <- X.test[,c(which(sbp!=0), variables2add)]
identical(colnames(X.train.select), colnames(X.test.select))

sbp <- sbp[match(colnames(X.train.select),rownames(sbp),nomatch = 0),,drop=FALSE]

saveRDS(list(sbp=sbp,X = X.train.select, X.test=X.test.select), 
        file = paste0("../Data/",file.end,".rds"))


