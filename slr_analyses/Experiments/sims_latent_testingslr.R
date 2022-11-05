# Purpose: demonstrate hierarchical spectral clustering with a threshold
# Date: 11/2/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Experiments/outputs/metrics"
b = 1

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs.R")
source("slr_analyses/Functions/slrs.R")
source("slr_analyses/Functions/codalasso.R")
source("slr_analyses/Functions/util.R")

# Tuning parameters###########################################################

# Settings to toggle with
settings.name = "ContinuousResponse"
hparam = "1se"
n = 100
p = 30
K = 10
nlam = 100
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_y = 0.1 # sigma (for y)
sigma_x = 0.1 # sigma_j (for x)
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.5 # 0.5
# theta.value = 1 # weight on a1 -- 1
c.value = 1 # a1 = c.value / k+ or c.value / k- or 0
a0 = 0 # 0
ulimit = 0.5

file.end = paste0(
  "_", settings.name,
  "_", paste0(
    paste(which(SBP.true == 1), collapse = ""), "v", 
    paste(which(SBP.true == -1), collapse = "")),
  "_hparam", hparam,
  "_dim", n, "x", p, 
  "_ulimit", ulimit,
  "_noisey", sigma_y, 
  "_noisex", sigma_x, 
  "_b0", b0, 
  "_b1", b1, 
  "_a0", a0, 
  "_c", c.value,
  "_sim", b,
  ".rds")

##############################################################################
# import data
data.tmp = readRDS(paste0(output_dir, "/data", file.end))
X = data.tmp$X
Y = data.tmp$Y
X.test = data.tmp$X.test
Y.test = data.tmp$Y.test
SBP.true = data.tmp$SBP.true
llc.coefs.true = data.tmp$llc.coefs.true
llc.coefs.non0 = data.tmp$llc.coefs.non0


##############################################################################
# slr
#   screen.method = "wald"
#   cluster.method = "spectral"
#   response.type = "continuous"
#   s0.perc = 0
#   zeta = 0
#   type.measure = "mse"
# -- fits a balance regression model with one balance
##############################################################################
slrspeccv = cv.slr(
  x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  nfolds = K, type.measure = "mse",
  scale = scaling, trace.it = FALSE)

x = X
y = Y
screen.method = "wald"
clsuter.method = "spectral"
response.type = "continuous"
s0.perc = 0
zeta = 0
nfolds = K
type.measure = "mse"
scale = scaling
trace.it = FALSE
threshold = NULL
foldid = NULL
weights = NULL
trace.it = FALSE

if(hparam == "min"){
  slrspec = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0,
    threshold = slrspeccv$threshold[slrspeccv$index["min",]],
    positive.slope = TRUE)
} else if(hparam == "1se"){
  slrspec = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "spectral",
    response.type = "continuous", s0.perc = 0, zeta = 0,
    threshold = slrspeccv$threshold[slrspeccv$index["1se",]],
    positive.slope = TRUE)
} else{
  stop("invalid hparam setting (method for selecting hyperparameter(s)).")
}

slrspec.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrspec.fullSBP) = colnames(X)
slrspec.fullSBP[match(
  names(slrspec$sbp), rownames(slrspec.fullSBP))] = slrspec$sbp
slrspec.coefs = getCoefsBM(
  coefs = coefficients(slrspec$fit), sbp = slrspec.fullSBP)

# compute metrics on the selected model #
# prediction error
slrspec.Yhat.test = slrspec.coefs$a0 +
  slr.fromContrast(X.test, slrspec.fullSBP) * slrspec.coefs$bm.coefs
slrspec.MSE.test = as.vector(crossprod(Y.test - slrspec.Yhat.test) / n)
# estimation accuracy, selection accuracy #
slrspec.metrics = getMetricsBM(
  est.llc.coefs = slrspec.coefs$llc.coefs,
  true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
  true.llc.coefs = llc.coefs.true,
  metrics = c("estimation", "selection"))

##############################################################################
# slr
#   screen.method = "wald"
#   cluster.method = "hierarchical"
#   response.type = "continuous"
#   s0.perc = 0
#   zeta = 0
#   type.measure = "mse"
# -- fits a balance regression model with one balance
##############################################################################
start.time = Sys.time()
slrhiercv = cv.slr(
  x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
  response.type = "continuous", s0.perc = 0, zeta = 0,
  nfolds = K, type.measure = "mse",
  scale = scaling, trace.it = FALSE)
if(hparam == "min"){
  slrhier = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0,
    threshold = slrhiercv$threshold[slrhiercv$index["min",]],
    positive.slope = TRUE)
} else if(hparam == "1se"){
  slrhier = slr(
    x = X, y = Y, screen.method = "wald", cluster.method = "hierarchical",
    response.type = "continuous", s0.perc = 0, zeta = 0,
    threshold = slrhiercv$threshold[slrhiercv$index["1se",]],
    positive.slope = TRUE)
} else{
  stop("invalid hparam setting (method for selecting hyperparameter(s)).")
}
end.time = Sys.time()
slrhier.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

slrhier.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrhier.fullSBP) = colnames(X)
slrhier.fullSBP[match(
  names(slrhier$sbp), rownames(slrhier.fullSBP))] = slrhier$sbp
slrhier.coefs = getCoefsBM(
  coefs = coefficients(slrhier$fit), sbp = slrhier.fullSBP)

# compute metrics on the selected model #
# prediction error
slrhier.Yhat.test = slrhier.coefs$a0 +
  slr.fromContrast(X.test, slrhier.fullSBP) * slrhier.coefs$bm.coefs
slrhier.MSE.test = as.vector(crossprod(Y.test - slrhier.Yhat.test) / n)
# estimation accuracy, selection accuracy #
slrhier.metrics = getMetricsBM(
  est.llc.coefs = slrhier.coefs$llc.coefs,
  true.sbp = SBP.true, non0.true.llc.coefs = llc.coefs.non0,
  true.llc.coefs = llc.coefs.true,
  metrics = c("estimation", "selection"))
