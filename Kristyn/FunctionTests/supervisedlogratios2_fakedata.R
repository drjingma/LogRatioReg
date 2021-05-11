getwd()

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)
library(balance) # for sbp.fromHclust()
library(rare)

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# Dr. Ma sources
library(Matrix)
library(glmnet)
library(compositions)
library(stats)
source("RCode/func_libs.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "supervisedlogratios2.R"))

# Method Settings
tol = 1e-4
nlam = 200
intercept = TRUE
K = 10
rho.type = "square"

# Simulation settings
numSims = 100
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
# which beta?
beta.settings = "new"
if(beta.settings == "old" | beta.settings == "linetal2014"){
  beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
} else{
  beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
}
non0.beta = (beta != 0)
sigma_eps = 0.5
seed = 1
muW = c(
  rep(log(p), 5), 
  rep(0, p - 5)
)
SigmaW = matrix(0, p, p)
for(i in 1:p){
  for(j in 1:p){
    SigmaW[i, j] = rho^abs(i - j)
  }
}

################################################################################
# Simulated Data #
################################################################################
set.seed(123)
# set.seed(1)
# set.seed(16) # slr0 is more sparse

# simulate training data #
# generate W
W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
# let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
V = exp(W)
rowsumsV = apply(V, 1, sum)
X = V / rowsumsV
epsilon = rnorm(n, 0, sigma_eps)
Z = log(X)
# generate Y
Y = Z %*% beta + epsilon

# simulate test data #
# simulate independent test set of size n
# generate W
W.test = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
# let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
V.test = exp(W.test)
rowsumsV.test = apply(V.test, 1, sum)
X.test = V.test / rowsumsV.test
epsilon.test = rnorm(n, 0, sigma_eps)
Z.test = log(X.test)
# generate Y
Y.test = Z.test %*% beta + epsilon.test

# get foldid like cv.glmnet does it
#   https://github.com/cran/glmnet/blob/master/R/cv.glmnet.R
nfolds = K
foldid = sample(rep(seq(nfolds), length = n))

################################################################################
# Applying SLR #
################################################################################

# apply the old slr #
slr0 = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
             rho.type = rho.type
             , foldid = foldid
             )

# choose lambda
lam.min.idx0 = which.min(slr0$cvm)
lam.min0 = slr0$lambda[lam.min.idx0]

# plot cvm vs lambda
# plot(slr0$lambda, slr0$cvm, #xlim = c(range(slr$lambda) + c(-0.1, 0.2)), 
#      col = rgb(0,0,0, alpha = 0.25))
ggplot(data.frame(lambda = slr0$lambda, cvm = slr0$cvm), 
       aes(x = lambda, y = cvm)) + 
  geom_path() + 
  geom_point(alpha = 0.25) + 
  theme_classic()
# what is the minimizing lambda?
lam.min0

# betahat?
thetahat0 = slr0$bet[, lam.min.idx0]
betahat0 = getBeta(thetahat0, btree = slr0$btree)
sum(betahat0 != 0)

################################################

# apply the new slr, with alpha = 1 -- should be equivalent to old slr #
slr = cvSLR2(y = Y, X = X, nlam = nlam, nfolds = K, alpha = 1
             , foldid = foldid
             , lambda = slr0$lambda
             )

# choose lambda
lam.min.idx = which.min(slr$cvm)
lam.min = slr$lambda[lam.min.idx]

# plot cvm vs lambda
# plot(slr$lambda, slr$cvm, col = rgb(0,0,0, alpha = 0.25))
ggplot(data.frame(lambda = slr$lambda, cvm = slr$cvm), 
       aes(x = lambda, y = cvm)) + 
  geom_path() + 
  geom_point(alpha = 0.25) + 
  theme_classic()
# what is the minimizing lambda?
lam.min

################################################

# apply the new slr with scaling = TRUE, with alpha = 1
slr2 = cvSLR2(y = Y, X = X, nlam = nlam, nfolds = K, alpha = 1
             , foldid = foldid
             , lambda = slr0$lambda
             , scaling = TRUE
)

# choose lambda
lam.min.idx2 = which.min(slr2$cvm)
lam.min2 = slr2$lambda[lam.min.idx2]

# plot cvm vs lambda
# plot(slr2$lambda, slr2$cvm, col = rgb(0,0,0, alpha = 0.25))
ggplot(data.frame(lambda = slr2$lambda, cvm = slr2$cvm), 
       aes(x = lambda, y = cvm)) + 
  geom_path() + 
  geom_point(alpha = 0.25) + 
  theme_classic()
# what is the minimizing lambda?
lam.min2

##############################################

# plot both together
ggplot(data.frame(lambda = c(slr0$lambda, slr$lambda, slr2$lambda), 
                  cvm = c(slr0$cvm, slr$cvm, slr2$cvm), 
                  type = c(rep("SLR", length(slr0$lambda)),
                           rep("SLR2", length(slr$lambda)), 
                           rep("SLR2scale", length(slr2$lambda)))),
       aes(x = lambda, y = cvm, color = type)) + 
  geom_path() + 
  geom_point(alpha = 0.25) + 
  theme_classic() + 
  geom_vline(xintercept = c(lam.min0, lam.min, lam.min2), alpha = 0.5)

# betahat?
betahat = slr$bet[[1]][, lam.min.idx]
thetahat = slr$thet[[1]][, lam.min.idx]
sum(betahat != 0)

# cbind(betahat0, betahat)
# sum(betahat0 != 0)
# sum(betahat != 0)

info = data.frame(
  "lambda min" = c(lam.min0, lam.min), 
  "min cvm" = c(min(slr0$cvm), min(slr$cvm)), 
  "S.hat" = c(sum(betahat0 != 0), sum(betahat != 0)), 
  "S0" = rep(sum(beta != 0), 2))
rownames(info) = c("SLR", "SLR2")
info
# library(knitr)
# kable(info, "latex")




