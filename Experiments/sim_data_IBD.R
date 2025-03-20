# update 3/20/2025

set.seed(1314)

library(selbal)
library(balance)

s1 <- 0.01
s2 <- 0.01
d <- 50 # can be 50 or 400
eigen.id <- 1

filename <- 'IBD_highCor'
filename <- paste0(filename, '_PC', eigen.id)

file.end <- paste0(
  "/",filename,
  "_sparsity", s1*100, 
  "_ratio", max(s1,s2)/min(s1,s2), 
  '_dim', d)

load(paste0("../Data/Franzosa_PRISM_UC.rda"))
load(paste0("../Data/Franzosa_Validation_UC.rda"))

# zero replacement
X = as.matrix(selbal::cmultRepl2(Franzosa_PRISM$X, zero.rep = "bayes"))
X.clr <- apply(X,1,function(s) log(s) - mean(log(s)))

# learn the empirical correlation
rho <- cor(t(X.clr))
v <- eigen(rho)$vectors[,eigen.id]

# determine active variables based on the largest positive and largest 
## negative values in the selected PC loading
sbp <- as.matrix((v <= quantile(v,s1)) - (v >= quantile(v,1-s2)))
summary(abs(rho)[sbp!=0,sbp!=0])
#pheatmap::pheatmap(rho[sbp!=0,sbp!=0])
print(table(sbp))
rownames(sbp) <- colnames(Franzosa_PRISM$X)

# randomly add inactive variables
variables2add <- sample(which(sbp==0),d)
counts <- Franzosa_PRISM$X[,c(which(sbp!=0), variables2add)]
counts.test <- Franzosa_Validation$X[,c(which(sbp!=0), variables2add)]

# adjust zeros if necessary
if (sum(counts==0)>0){
  X <- as.matrix(selbal::cmultRepl2(counts, zero.rep = "bayes"))
}
if (sum(counts.test==0)>0){
  X.test <- as.matrix(selbal::cmultRepl2(counts.test, zero.rep = "bayes"))
}
sbp <- sbp[match(colnames(X),rownames(sbp),nomatch = 0),,drop=FALSE]

# save data
saveRDS(list(sbp=sbp,X=X, X.test=X.test), 
        file = paste0("../Data/",file.end,".rds"))





