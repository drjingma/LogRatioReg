# last updated: 03/24/2021
# simulate a data set, fit Complasso and SLR, and plot the solution path

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
numSims = 12
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
# indices.theta = c(25, 100, 150, 199) # some index between 1 and p - 1
## maybe try different magnitudes? different signs? ############################
indices.theta = sample(1:(p - 1), 5, replace = FALSE)
values.theta = NULL

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

# set.seed(1)
sims = foreach(
  b = 1:numSims
) %dorng% {
  library(limSolve)
  library(mvtnorm)
  library(Matrix)
  library(glmnet)
  library(compositions)
  library(stats)
  library(balance) # for sbp.fromHclust()
  
  source("RCode/func_libs.R")
  source(paste0(functions_path, "supervisedlogratios.R"))
  
  nlam = 200
  
  # simulate some data #
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
  # get theta
  theta = rep(0, p - 1)
  if(is.null(values.theta)){
    theta[indices.theta] = 1
  } else{
    if(length(indices.theta) == length(values.theta)){
      stop("indices.theta does not have same length as values.theta")
    }
    theta[indices.theta] = values.theta
  }
  theta = as.matrix(theta)
  # get beta
  beta = getBeta(theta, sbp = sbp)
  # generate Y
  Y = Xb %*% theta + epsilon
  
  # apply compositional lasso
  complasso = cv.func(
    method="ConstrLasso", y = Y, x = log(X), Cmat = matrix(1, p, 1), 
    nlam = nlam, nfolds = K, tol = tol, intercept = intercept)
  
  TPR.cl = rep(NA, nlam)
  S.hat.cl = rep(NA, nlam)
  for(l in 1:nlam){
    a0 = complasso$int[l]
    betahat = complasso$bet[, l]
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR.cl[l] = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat.cl[l] = sum(non0.betahat)
  }
  
  # apply supervised log-ratios
  slr = cvSLR(y = Y, X = X, nlam = nlam, nfolds = K, intercept = intercept, 
              rho.type = rho.type)
  
  btree.slr = slr$btree
  # calculate TPR for supervised log-ratios
  nlam.slr = length(slr$lambda)
  TPR.slr = rep(NA, nlam.slr)
  S.hat.slr = rep(NA, nlam.slr)
  for(l in 1:nlam.slr){
    a0 = slr$int[l]
    thetahat = slr$bet[, l]
    betahat = getBeta(thetahat, btree.slr)
    # TPR
    non0.beta = (beta != 0)
    non0.betahat = (betahat != 0)
    TPR.slr[l] = sum((non0.beta == non0.betahat) & non0.betahat) / sum(non0.beta)
    S.hat.slr[l] = sum(non0.betahat)
  }
  
  list(fit.cl = complasso, TPR.cl = TPR.cl, S.hat.cl = S.hat.cl, 
       fit.slr = slr, TPR.slr = TPR.slr, S.hat.slr = S.hat.slr)
}


# plot one simulation results #

sim.idx = 5
sim.tmp = sims[[sim.idx]]

library(ggplot2)
slr.gg = data.frame(
  lambda = sim.tmp$fit.slr$lambda,
  S.hat = sim.tmp$S.hat.slr,
  TPR = sim.tmp$TPR.slr,
  ID = "SLR")
slr.gg = slr.gg[order(slr.gg$S.hat), ]
cl.gg = data.frame(
  lambda = sim.tmp$fit.cl$lambda,
  S.hat = sim.tmp$S.hat.cl,
  TPR = sim.tmp$TPR.cl,
  ID = "CompLasso")
cl.gg = cl.gg[order(cl.gg$S.hat), ]
data.gg = rbind(slr.gg, cl.gg)
ggplot(data.gg, aes(x = S.hat, y = TPR)) +
  geom_line(aes(color = ID), size = 1) +
  theme_bw()



# note: lambda values are different in different simulations
#   this makes sense, because simulated data are different!
#   but that's okay, we are only interested in S.hat and TPR...



# organize the simulation results to have more easily-accessable matrices #
# cl
fit.cl.list = list()
TPR.cl.mat = matrix(NA, nlam, numSims)
S.hat.cl.mat = matrix(NA, nlam, numSims)
# slr
fit.slr.list = list()
TPR.slr.mat = matrix(NA, nlam, numSims)
S.hat.slr.mat = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  sim.tmp = sims[[i]]
  # cl
  fit.cl.list[[i]] = sim.tmp$fit.cl
  TPR.cl.mat[, i] = sim.tmp$TPR.cl
  S.hat.cl.mat[, i] = sim.tmp$S.hat.cl
  # slr
  fit.slr.list[[i]] = sim.tmp$fit.slr
  TPR.slr.mat[, i] = sim.tmp$TPR.slr
  S.hat.slr.mat[, i] = sim.tmp$S.hat.slr
}
# length(unique(as.list(as.data.frame(S.hat.slr))))
# View(S.hat.slr)



# plot each simulations' results #
library(reshape2)
# cl.stuff
S.hat.cl.mlt = melt(S.hat.cl.mat)
colnames(S.hat.cl.mlt) = c("LambdaIndex", "SimIndex", "S.hat")
TPR.cl.mlt = melt(TPR.cl.mat)
colnames(TPR.cl.mlt) = c("LambdaIndex", "SimIndex", "TPR")
# all.equal(S.hat.cl.mlt$SimIndex, TPR.cl.mlt$SimIndex)
cl.gg.complete = data.frame(
  "SimIndex" = S.hat.cl.mlt$SimIndex, 
  "S.hat" = S.hat.cl.mlt$S.hat, 
  "TPR" = TPR.cl.mlt$TPR)
cl.gg.complete$Type = "CompLasso"
# slr stuff
S.hat.slr.mlt = melt(S.hat.slr.mat)
colnames(S.hat.slr.mlt) = c("LambdaIndex", "SimIndex", "S.hat")
TPR.slr.mlt = melt(TPR.slr.mat)
colnames(TPR.slr.mlt) = c("LambdaIndex", "SimIndex", "TPR")
# all.equal(S.hat.slr.mlt$SimIndex, TPR.slr.mlt$SimIndex)
slr.gg.complete = data.frame(
  "SimIndex" = S.hat.slr.mlt$SimIndex, 
  "S.hat" = S.hat.slr.mlt$S.hat, 
  "TPR" = TPR.slr.mlt$TPR)
slr.gg.complete$Type = "SLR"
gg.complete = rbind(slr.gg.complete, cl.gg.complete)
ggplot(gg.complete, aes(x = S.hat, y = TPR)) + 
  facet_wrap(vars(SimIndex)) +
  geom_line(aes(color = Type), size = 1) +
  theme_bw()


##### not sure how to average over the simulations ... #####
# # for each possible value of S.hat (from 0 to p - 1), 
# #   average over the corresponding values of TPR
# 
# # first average over the corresponding values of TPR within each sim
# # this will give all simulations the same x-axis value: S.hat.vals
# # note: can expect some NAs
# S.hat.vals = 0:(p - 1)
# TPR.simwise.avgs = matrix(NA, length(S.hat.vals), numSims)
# for(i in 1:numSims){
#   
# }
# 
# 
# 
# 
# S.hat.vals = 0:(p - 1)
# TPR.avg.vals.cl = rep(NA, p)
# TPR.avg.vals.slr = rep(NA, p)
# idx = 1
# for(i in S.hat.vals){
#   # cl
#   # get the locations of the S.hat.vals[i]'s in S.hat.cl matrix
#   match.cl = (S.hat.cl == S.hat.vals[i])
#   which.cl = apply(match.cl, 1, FUN = which)
#   # find the corresponding TPRs
#   TPR.vals.cl.i = c() # append
#   for(j in 1:length(which.cl)){
#     TPR.vals.cl.i = c(TPR.vals.cl.i, TPR.cl.mat[j, which.cl[[j]]])
#   }
#   # average them
#   TPR.avg.vals.cl[idx] = mean(TPR.vals.cl.i)
#   # slr
#   # get the locations of the S.hat.vals[i]'s in S.hat.cl matrix
#   match.slr = (S.hat.slr == S.hat.vals[i])
#   which.slr = apply(match.slr, 1, FUN = which)
#   # find the corresponding TPRs
#   TPR.vals.slr.i = c() # append
#   for(j in 1:length(which.cl)){
#     TPR.vals.slr.i = c(TPR.vals.slr.i, TPR.slr.mat[j, which.slr[[j]]])
#   }
#   # average them
#   TPR.avg.vals.slr[idx] = mean(TPR.vals.slr.i)
#   idx = idx + 1
# }
