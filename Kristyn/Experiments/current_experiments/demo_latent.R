# Purpose: look at aitchison matrix
# Date: 5/24/2022
rm(list=ls())

set.seed(1234)

################################################################################
# libraries and settings

library(mvtnorm)

library(Matrix)
library(glmnet)

library(balance)

source("RCode/func_libs_1.R") # for classo to work

source("Kristyn/Functions/slr.R")
source("Kristyn/Functions/util.R")
source("Kristyn/Functions/slrscreen.R")

# Tuning parameters###########################################################

# Settings to toggle with
sigma.settings = "latentVarModel"
n = 100
p = 30
K = 10
nlam = 100
neta = p
intercept = TRUE
scaling = TRUE
tol = 1e-4
sigma_eps1 = 0.1 # sigma (for y)
sigma_eps2 = 0.1 # sigma_j (for x)
SBP.true = matrix(c(1, 1, 1, -1, -1, -1, rep(0, p - 6)))
# SBP.true = matrix(c(1, 1, 1, 1, -1, rep(0, p - 5)))
ilrtrans.true = getIlrTrans(sbp = SBP.true, detailed = TRUE)
# ilrtrans.true$ilr.trans = transformation matrix (used to be called U) 
#   = ilr.const*c(1/k+,1/k+,1/k+,1/k-,1/k-,1/k-,0,...,0)
b0 = 0 # 0
b1 = 0.5 # 1, 0.5, 0.25
theta.value = 1 # weight on a1 -- 1
a0 = 0 # 0

##############################################################################
# generate data
# get latent variable
U.all = matrix(runif(min = -0.5, max = 0.5, 2 * n), ncol = 1)
# simulate y from latent variable
y.all = as.vector(b0 + b1 * U.all + rnorm(2 * n) * sigma_eps1)
# simulate X: 
epsj.all = matrix(rnorm(2 * n * (p - 1)), nrow = (2 * n)) * sigma_eps2
a1 = theta.value * ilrtrans.true$ilr.trans[-p] 
#   alpha1j = {
#     c1=theta*ilr.const/k+   if j \in I+
#     -c2=-theta*ilr.const/k-  if j \in I-
#     0                       o/w
#   }
alrXj.all = a0 + U.all %*% t(a1) + epsj.all #log(Xj/Xp) =alpha0j+alpha1j*U+epsj
X.all <- alrinv(alrXj.all)
colnames(X.all) = paste0('s', 1:p)

# subset out training and test sets
X = X.all[1:n, ]
X.test = X.all[-(1:n), ]
Y <- y.all[1:n]
Y.test <- y.all[-(1:n)]

# about beta
non0.beta = as.vector(SBP.true != 0)
bspars = sum(non0.beta)
# solve for beta
c1plusc2 = theta.value * sum(abs(unique(ilrtrans.true$ilr.trans)))
beta.true = (b1 / (ilrtrans.true$const * c1plusc2)) * 
  as.vector(ilrtrans.true$ilr.trans)

# aitchison similarity matrix
X.reduced <- X[, which(SBP.true != 0)] # reduced data matrix
p.reduced = sum(SBP.true != 0)
Aitchison.variation <- matrix(0, p.reduced, p.reduced)
rownames(Aitchison.variation) <- colnames(Aitchison.variation) <- colnames(
  X.reduced)
for (j in 1:p.reduced){
  for (k in 1:p.reduced){
    if (k==j){next}
    else {
      Aitchison.variation[j,k] <- stats::var(
        log(X.reduced[,j]) - log(X.reduced[,k])) # Aitchison variation
    }
  }
}
Aitchison.similarity <- max(Aitchison.variation) - Aitchison.variation 
pheatmap::pheatmap(
  Aitchison.similarity, cluster_rows = FALSE, cluster_cols = FALSE, 
  fontsize = 12)
?pheatmap
##############################################################################
# slr method with screening step
#   method = "wald"
#   response.type = "continuous"
#   s0.perc = 0
#   zeta = 0
#   type.measure = "mse"
# -- fits a balance regression model with one balance
##############################################################################
start.time = Sys.time()
slrscreen0cv = cv.slr.screen(
  x = X, y = Y, method = "wald", 
  response.type = "continuous", s0.perc = 0, zeta = 0, 
  nfolds = K, type.measure = "mse", 
  parallel = FALSE, scale = scaling, trace.it = FALSE)
slrscreen0 = slr.screen(
  x = X, y = Y, method = "wald", 
  response.type = "continuous", s0.perc = 0, zeta = 0, 
  threshold = slrscreen0cv$threshold[slrscreen0cv$index["1se",]])
end.time = Sys.time()
slrscreen0.timing = difftime(
  time1 = end.time, time2 = start.time, units = "secs")

slrscreen0.fullSBP = matrix(0, nrow = p, ncol = 1)
rownames(slrscreen0.fullSBP) = colnames(X)
slrscreen0.fullSBP[match(
  names(slrscreen0$sbp), rownames(slrscreen0.fullSBP))] = slrscreen0$sbp

slrscreen0.coefs = getCoefsBM(
  coefs = coefficients(slrscreen0$fit), sbp = slrscreen0.fullSBP)

# compute metrics on the selected model #
slrscreen0.metrics = getMetricsBM(
  y.train = Y, y.test = Y.test,
  ilrX.train = getIlrX(X, sbp = slrscreen0.fullSBP),
  ilrX.test = getIlrX(X.test, sbp = slrscreen0.fullSBP),
  n.train = n, n.test = n,
  thetahat0 = slrscreen0.coefs$a0, thetahat = slrscreen0.coefs$bm.coefs,
  betahat = slrscreen0.coefs$llc.coefs,
  true.sbp = SBP.true, non0.true.beta = non0.beta,
  true.beta = beta.true)













