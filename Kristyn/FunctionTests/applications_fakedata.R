setwd("/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg")

# Kristyn sources
source("Kristyn/Functions/classic_lasso.R")
source("Kristyn/Functions/compositional_lasso.R")
source("Kristyn/Functions/supervisedlogratios.R")

# Dr. Ma sources
source("RCode/func_libs.R")
source("COAT-master/coat.R")

# libraries
library(mvtnorm)

# settings to generate data, n < p
n = 100
p = 200
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

################################################################################
################################################################################
################################################################################

################################################################################
#  compositional lasso
################################################################################

# y = log(X) beta + epsilon
# subject to constraint sum_{i = 1}^p beta_i == 1

# betahat from compositional Lasso
test_cvLASSO = cvCompositionalLASSO(Z ,Y, lambda_seq = NULL, n_lambda = 30, k = 10)
# plot(test_cvLASSO$cvm ~ test_cvLASSO$lambda_seq, type = "l")
betahatcompositional = test_cvLASSO$beta_mat[, test_cvLASSO$cvm_idx]
betahat0compositional = test_cvLASSO$beta0_vec[test_cvLASSO$cvm_idx]

# evaluate betahat
sum(beta - betahatcompositional)^2
beta[1:10]
round(betahat0compositional, 2)
round(betahatcompositional[1:10], 2)


compositionalLassofit = function(x){
  betahat0compositional + log(x) %*% betahatcompositional
}

################################################################################
# classic lasso
################################################################################

# betahat from classic Lasso
test_cvLASSO = cvLASSO(Z ,Y, lambda_seq = NULL, n_lambda = 30, k = 10)
# plot(test_cvLASSO$cvm ~ test_cvLASSO$lambda_seq, type = "l")
betahatclassic = test_cvLASSO$beta_mat[, test_cvLASSO$cvm_idx]
betahat0classic = test_cvLASSO$beta0_vec[test_cvLASSO$cvm_idx]

# evaluate betahat
sum(beta[1:10] - betahatclassic)^2
beta[1:10]
round(betahat0classic, 2)
round(betahatclassic[1:10], 2)

classicLassofit = function(x){
  betahat0classic + log(x) %*% betahatclassic
}
