################################################################################
# same simulated data is used to test 
#   classic lasso
#   compositional lasso
################################################################################

setwd("/home/kristyn/Documents/research/supervisedlogratios/lasso")
source("functions/classic_lasso.R")

library(mvtnorm)

# generate data, n < p
n = 100
p = 200

# settings
rho = 0.2
beta = c(1, 0.4, 1.2, -1.5, -0.8, 0.3, rep(0, p - 6))
sigma_epsilon = 0.5
seed = 1

# simulate data
set.seed(seed)
# sample W from log-normal distribution
# if V has normal distribution, then W = exp(V) has log-normal distribution
muW = c(
  rep(log(p), 5), 
  rep(0, p - 5)
)
SigmaW = matrix(0, p, p)
for(j in 1:p){
  for(k in j:p){
    SigmaW[j, k] = rho^(k - j)
  }
}
SigmaW = SigmaW + t(SigmaW) - diag(diag(SigmaW))
V = rmvnorm(n = 100, mean = muW, sigma = SigmaW)
W = exp(V)
rowsumsW = apply(W, 1, sum)
X = W / rowsumsW
epsilon = rnorm(n, 0, sigma_epsilon)
Z = log(X)
Y = Z %*% beta + epsilon

# Lasso
# y = Z beta* + epsilon
# beta_p^* = sum_{i = 1}^{p - 1} beta^*
# so that sum_{i = 1}^p beta^* == 1

stdZY = standardizeXY(Z, Y)

# lambda_max gives beta = 0
lambda_max = max(abs(crossprod(stdZY$Xtilde, stdZY$Ytilde) / n))
# 
lambda = lambda_max
testLasso = fitLASSOstandardized(
  stdZY$Xtilde, stdZY$Ytilde, lambda)
plot(testLasso$objectives ~ c(1:length(testLasso$objectives)), type = "l")
testLasso$beta
# slightly smaller lambda
lambda = lambda_max * 0.3
testLasso = fitLASSOstandardized(
  stdZY$Xtilde, stdZY$Ytilde, lambda)
plot(testLasso$objectives ~ c(1:length(testLasso$objectives)), type = "l")

# back scaling and centering to get original intercept and coefficient vector
tildebeta = testLasso$beta
betahat = diag(stdZY$weights) %*% tildebeta
betahat0 = stdZY$Ymean - crossprod(stdZY$Xmeans, betahat)
beta[1:10]
betahat[1:10]
betahat0 # why?? because Y ~ log(X) instead of Y ~ X?

# error of back-transformed beta is slightly smaller
sum(abs(beta - betahat)) # 2.269
sum(abs(beta - tildebeta)) # 2.577
sum(betahat)

################################################################################

# test fitLASSOstandardized_seq for lambda = 0 (LS soln)
testLASSOstd_seq = fitLASSOstandardized_seq(stdZY$Xtilde, stdZY$Ytilde)
numNonzerosTot = colSums(testLASSOstd_seq$beta_mat == 0) # total number of non-zero elements in beta
plot(numNonzerosTot ~ testLASSOstd_seq$lambda_seq, type = "l") # checks out

# test fitLASSO
testLASSO = fitLASSO(Z, Y)
numNonzerosTot = colSums(testLASSO$beta_mat == 0) # total number of non-zero elements in beta
plot(numNonzerosTot ~ testLASSO$lambda_seq, type = "l") # checks out

# test cvLASSO
test_cvLASSO = cvLASSO(Z ,Y, lambda_seq = NULL, n_lambda = 30, k = 10)
plot(test_cvLASSO$cvm ~ test_cvLASSO$lambda_seq, type = "l")

# Use CV lambda
betahatclassic = test_cvLASSO$beta_mat[, test_cvLASSO$cvm_idx]
betahat0classic = test_cvLASSO$beta0_vec[test_cvLASSO$cvm_idx]

beta[1:10]
round(betahat0classic, 2)
round(betahatclassic[1:10], 2)
