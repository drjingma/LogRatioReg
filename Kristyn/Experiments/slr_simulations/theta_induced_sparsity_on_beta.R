# understand how the index of the nonzero element of theta (where all other
# elements = 0) determines the sparsity of beta

getwd()
output_dir = "Kristyn/Experiments/slr_simulations/output"

# libraries
library(mvtnorm)
library(stats) # for hclust()
library(balance) # for sbp.fromHclust()

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

# Method Settings
tol = 1e-4
nlam = 200
intercept = TRUE
K = 5
rho.type = "square"

# Simulation settings
numSims = 100
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
# generate.theta = 1 # 1 = sparse beta, 2 = not-sparse beta ####################
sigma_eps = 0.2
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
# Simulations #
################################################################################

set.seed(1)

# simulate training data #
# generate W
W = rmvnorm(n = n, mean = muW, sigma = SigmaW) # n x p
# let X = exp(w_ij) / (sum_k=1:p w_ik) ~ Logistic Normal (the covariates)
V = exp(W)
rowsumsV = apply(V, 1, sum)
X = V / rowsumsV
sbp = sbp.fromPBA(X) # contrasts matrix, a.k.a. sbp matrix
U = getU(sbp = sbp) # U
epsilon = rnorm(n, 0, sigma_eps)
Xb = log(X) %*% U # ilr(X)
# each theta has one nonzero element
theta_mat = diag(p - 1)
beta_mat = apply(theta_mat, 2, function(x) getBeta(x, sbp = sbp))
non0beta_vec = apply(beta_mat, 2, function(x) sum(x != 0))

library(ggplot2)
ggplot(data = data.frame(x = 1:(p - 1), y = non0beta_vec), aes(x = x, y = y)) + 
  geom_path() + theme_bw() + labs(x = "index of nonzero theta_j", y = "|beta|")




