# last updated: 04/26/2021
# simulate a data set using settings in Lin et al 2014, 
# fit Complasso and SLR, and plot the solution path.

getwd()
output_dir = "Kristyn/Experiments/complasso_simulations/output"

# libraries
library(limSolve) # for constrained lm
library(mvtnorm)
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

# ggplot
library(ggplot2)

# Method Settings
tol = 1e-4
nlam = 3 # for testing
intercept = TRUE
K = 10
rho.type = "square"

# Simulation settings
numSims = 100
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
# which beta?
beta.settings = "old"
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
# if simulations are saved, read them in

nlam = 200
if(beta.settings == "old" | beta.settings == "linetal2014"){
  sims3 = readRDS(paste0(output_dir,
                         "/solpaths_old",
                         "_dim", n, "x", p,
                         "_rho", rho,
                         "_int", intercept,
                         "_K", K,
                         "_seed", rng.seed,
                         "_numSims", numSims,
                         ".rds"))
} else{
  sims3 = readRDS(paste0(output_dir,
                         "/solpaths",
                         "_dim", n, "x", p,
                         "_rho", rho,
                         "_int", intercept,
                         "_K", K,
                         "_seed", rng.seed,
                         "_numSims", numSims,
                         ".rds"))
}


# organize the simulation results to have more easily-accessable matrices ######
# cl
fit.cl.list = list()
TPR.cl.mat = matrix(NA, nlam, numSims)
S.hat.cl.mat = matrix(NA, nlam, numSims)
lambda.cl.mat = matrix(NA, nlam, numSims)
# slr
fit.slr.list = list()
TPR.slr.mat = matrix(NA, nlam, numSims)
S.hat.slr.mat = matrix(NA, nlam, numSims)
lambda.slr.mat = matrix(NA, nlam, numSims)
for(i in 1:numSims){
  sim.tmp = sims3[[i]]
  # cl
  fit.cl.list[[i]] = sim.tmp$fit.cl
  TPR.cl.mat[, i] = sim.tmp$TPR.cl
  S.hat.cl.mat[, i] = sim.tmp$S.hat.cl
  lambda.cl.mat[, i] = sim.tmp$fit.cl$lambda
  # slr
  fit.slr.list[[i]] = sim.tmp$fit.slr
  TPR.slr.mat[, i] = sim.tmp$TPR.slr
  S.hat.slr.mat[, i] = sim.tmp$S.hat.slr
  lambda.slr.mat[, i] = sim.tmp$fit.slr$lambda
}



# average TPR and S.hat ########################################################

dim(S.hat.cl.mat)
S.hat.cl.avg = apply(S.hat.cl.mat, 1, mean, na.rm = TRUE)
TPR.cl.avg = apply(TPR.cl.mat, 1, mean, na.rm = TRUE)
S.hat.slr.avg = apply(S.hat.slr.mat, 1, mean, na.rm = TRUE)
TPR.slr.avg = apply(TPR.slr.mat, 1, mean, na.rm = TRUE)

# complasso stuff
cl.gg.complete = data.frame(
  "S.hat" = S.hat.cl.avg,
  "TPR" = TPR.cl.avg)
cl.gg.complete$Type = "CompLasso"
# slr stuff
slr.gg.complete = data.frame(
  "S.hat" = S.hat.slr.avg,
  "TPR" = TPR.slr.avg)
slr.gg.complete$Type = "SLR"
# ggplot
gg.complete = rbind(slr.gg.complete, cl.gg.complete)
gg.complete$Type = factor(gg.complete$Type, levels = c("CompLasso", "SLR"))
ggplot(gg.complete, aes(x = S.hat, y = TPR)) +
  geom_line(aes(color = Type), size = 1) +
  xlim(0, 40) +
  theme_bw()


