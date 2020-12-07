# workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
# setwd(workdir)

# libraries
library(mvtnorm)
library(balance)
library(selbal)
library(microbenchmark)
library(ggplot2)
library(logratiolasso) # bates & tibshirani 2019

# Dr. Ma sources
source("RCode/func_libs.R")
source("COAT-master/coat.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "classic_lasso.R"))
source(paste0(functions_path, "compositional_lasso.R"))
source(paste0(functions_path, "supervisedlogratios.R"))
source(paste0(functions_path, "coat.R"))
source(paste0(functions_path, "principlebalances.R"))
source(paste0(functions_path, "selbal.R"))

# data
# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
DataFolder <- "/Data/"
load(paste0("Data/", "BMI.rda"))
log.X.prop = log(X.prop)
# dim(raw_data) # 98 x 89
# dim(X) # 98 x 87
# dim(X.prop) # 98 x 87
n = dim(X)[1]
num.genera = dim(X)[2]

# combo data
combo <- read.csv("Kristyn/Data/combo_ffq_adj.txt", sep = " ")
totalcaloricfatintake = combo[, c("calor", "tfat")]

# testing with centered and/or scaled y
# y = scale(y, center = T, scale = T)

final.selected = c(
  "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes",
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella")
beta.paper = c(-0.76, -1.35, 0.61, 1.50)

################################################################################
# CONSTRAINED LINEAR MODEL
################################################################################
par(mfrow = c(2, 3))
# constrained - no centering or scaling: 
# beta-bar = beta-hat - (X'X)^(-1) 1 [1' (X'X)^(-1) 1]^(-1) (1' beta-hat)
Xmat = log.X.prop[, final.selected]
Xmat = as.matrix(cbind(Xmat, totalcaloricfatintake))
beta.hat = solve(crossprod(Xmat), crossprod(Xmat, y))
# check that it matches lm() #
beta.hat.lm = coefficients(lm(y ~ -1 + Xmat))
all.equal(as.numeric(beta.hat), as.numeric(beta.hat.lm)) # they match, as expected
# compute beta.bar, constrained fit
Q = as.matrix(c(rep(1, dim(Xmat)[2] - dim(totalcaloricfatintake)[2]), 
                rep(0, dim(totalcaloricfatintake)[2])))
beta.bar = beta.hat - solve(crossprod(Xmat), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat), Q)), crossprod(Q, beta.hat))
# beta.bar = beta.bar[-c(5, 6)]
# sum(beta.bar)
# beta.bar
# l2 = sqrt(crossprod(beta.bar))
# dist.t = sqrt(sum((beta.bar - beta.paper)^2))
# l2
# dist.t
yhat = Xmat %*% beta.bar
xyrange = range(y, yhat)
plot(main = "no intercept", x = y, y = yhat, xlim = xyrange, ylim = xyrange); abline(0, 1, lty = 2)


# constrained - scaling
x.sd = apply(Xmat,2,sd)
Xmat.s = scale(Xmat, center=F, scale=x.sd)
beta.hat.s = solve(crossprod(Xmat.s), crossprod(Xmat.s, y))
beta.bar.s = beta.hat.s - solve(crossprod(Xmat.s), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.s), Q)), crossprod(Q, beta.hat.s))
# beta.bar.s = beta.bar.s[-c(5, 6)]
# sum(beta.bar.s)
# beta.bar.s
# l2.s = sqrt(crossprod(beta.bar.s))
# dist.t.s = sqrt(sum((beta.bar.s - beta.paper)^2))
# l2.s
# dist.t.s
yhat = Xmat.s %*% beta.bar.s
xyrange = range(y, yhat)
plot(main = "no int, scaled", x = y, y = yhat, xlim = xyrange, ylim = xyrange); abline(0, 1, lty = 2)
#
beta.bar.s2 = beta.bar.s / x.sd
# sum(beta.bar.s2) # not satisfied anymore
# beta.bar.s2
# l2.s2 = sqrt(crossprod(beta.bar.s2))
# dist.t.s2 = sqrt(sum((beta.bar.s2 - beta.paper)^2))
# l2.s2
# dist.t.s2
yhat = Xmat.s %*% beta.bar.s2
xyrange = range(y, yhat)
plot(main = "no int, back-scaled", x = y, y = yhat, xlim = xyrange, ylim = xyrange); abline(0, 1, lty = 2)

# constrained - with intercept
y.mean = mean(y)
y.c = y - y.mean
x.mean = colMeans(Xmat)
Xmat.c = scale(Xmat, center=x.mean, scale=F)
beta.hat.c = solve(crossprod(Xmat.c), crossprod(Xmat.c, y.c))
# check that it matches lm() #
beta.hat.lmint = coefficients(lm(y ~ Xmat))
# coefficients when centering match when not-centering (just not intercept, since centering gets rid of intercept)
all.equal(as.numeric(coefficients(lm(y ~ Xmat))[-1]), 
          as.numeric(coefficients(lm(y.c ~ Xmat.c))[-1]))
all.equal(as.numeric(beta.hat.c), as.numeric(beta.hat.lmint)[-1]) # they match, as expected
beta.hat.lmint0 = y.mean - as.vector(x.mean%*%beta.hat.c) # also match
# compute beta.bar, constrained fit
beta.bar.c = beta.hat.c - solve(crossprod(Xmat.c), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.c), Q)), crossprod(Q, beta.hat.c))
beta.bar.c0 <- y.mean - as.vector(x.mean%*%beta.bar.c)
# constrained - with intercept, including in Xmat
Xmat.c2 = cbind(1, Xmat)
beta.hat.c2 = solve(crossprod(Xmat.c2), crossprod(Xmat.c2, y))
Q2 = c(0, Q)
beta.bar.c2 = beta.hat.c2 - solve(crossprod(Xmat.c2), Q2) %*% 
  solve(crossprod(Q2, solve(crossprod(Xmat.c2), Q2)), crossprod(Q2, beta.hat.c2))
sum(beta.bar.c2[-c(1, 6, 7)])
beta.bar.c2 
all.equal(c(beta.bar.c0, beta.bar.c), as.numeric(beta.bar.c2)) # this matches the above
# beta.bar.c = beta.bar.c[-c(5, 6)]
# sum(beta.bar.c)
# beta.bar.c
# l2.c = sqrt(crossprod(beta.bar.c))
# dist.t.c = sqrt(sum((beta.bar.c - beta.paper)^2))
# l2.c
# dist.t.c
yhat = beta.bar.c0 + Xmat.c %*% beta.bar.c
xyrange = range(y, yhat)
plot(main = "with intercept", x = y, y = yhat, xlim = xyrange, ylim = xyrange); abline(0, 1, lty = 2)

# constrained - with intercept and scaling (standardizing, bc need to center for intercept)
Xmat.cs = scale(Xmat, center=x.mean, scale=x.sd)
beta.hat.cs = solve(crossprod(Xmat.cs), crossprod(Xmat.cs, y.c))
beta.bar.cs = beta.hat.cs - solve(crossprod(Xmat.cs), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.cs), Q)), crossprod(Q, beta.hat.cs))
# beta.bar.cs = beta.bar.cs[-c(5, 6)]
# sum(beta.bar.cs)
# beta.bar.cs
beta.bar.cs0 <- y.mean - as.vector(x.mean%*%beta.bar.cs)
beta.bar.cs0
# l2.cs = sqrt(crossprod(beta.bar.cs))
# dist.t.cs = sqrt(sum((beta.bar.cs - beta.paper)^2))
# l2.cs
# dist.t.cs
yhat = beta.bar.cs0 + Xmat.cs %*% beta.bar.cs
xyrange = range(y, yhat)
plot(main = "with int, std.", x = y, y = yhat, xlim = xyrange, ylim = xyrange); abline(0, 1, lty = 2)
#
beta.bar.cs2 = beta.bar.cs / x.sd
beta.bar.cs2
sum(beta.bar.cs2) # not satisfied anymore
beta.bar.cs02 <- y.mean - as.vector(x.mean%*%beta.bar.cs2)
beta.bar.cs02
# sum(beta.bar.cs02, beta.bar.cs2)
# beta.bar.cs2 = beta.bar.cs2[-c(5, 6)]
# sum(beta.bar.cs2) # not satisfied anymore
# l2.cs2 = sqrt(crossprod(beta.bar.cs2))
# dist.t.cs2 = sqrt(sum((beta.bar.cs2 - beta.paper)^2))
# l2.cs2
# dist.t.cs2
yhat = beta.bar.cs02 + Xmat.cs %*% beta.bar.cs2
xyrange = range(y, yhat)
plot(main = "with in, back-std.", x = y, y = yhat, xlim = xyrange, ylim = xyrange); abline(0, 1, lty = 2)

################################################################################
# CONSTRAINED LINEAR MODEL -- SUBCOMPOSITION
################################################################################
X.sel = X[, final.selected]
X.prop.sel = sweep(X.sel, MARGIN = 1, STATS = rowSums(X.sel), FUN = "/")
log.X.prop.sel = log(X.prop.sel)

# constrained - no centering or scaling: 
# beta-bar = beta-hat - (X'X)^(-1) 1 [1' (X'X)^(-1) 1]^(-1) (1' beta-hat)
Xmat.sub = log.X.prop.sel[, final.selected]
beta.hat.sub = solve(crossprod(Xmat.sub), crossprod(Xmat.sub, y))
# check that it matches lm() #
beta.hat.lm = coefficients(lm(y ~ -1 + Xmat.sub))
all.equal(as.numeric(beta.hat.sub), as.numeric(beta.hat.lm)) # they match, as expected
# compute beta.bar, constrained fit
Q = as.matrix(rep(1, dim(Xmat.sub)[2]))
beta.bar.sub = beta.hat.sub - solve(crossprod(Xmat.sub), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.sub), Q)), crossprod(Q, beta.hat.sub))
sum(beta.bar.sub)
beta.bar.sub
l2.sub = sqrt(crossprod(beta.bar.sub))
dist.t.sub = sqrt(sum((beta.bar.sub - beta.paper)^2))
l2.sub
dist.t.sub

# constrained - scaling
x.sd.sub = apply(Xmat.sub,2,sd)
Xmat.sub.s = scale(Xmat.sub, center=F, scale=x.sd.sub)
beta.hat.sub.s = solve(crossprod(Xmat.sub.s), crossprod(Xmat.sub.s, y))
beta.bar.sub.s = beta.hat.sub.s - solve(crossprod(Xmat.sub.s), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.sub.s), Q)), crossprod(Q, beta.hat.sub.s))
sum(beta.bar.sub.s)
beta.bar.sub.s
l2.sub.s = sqrt(crossprod(beta.bar.sub.s))
dist.t.sub.s = sqrt(sum((beta.bar.sub.s - beta.paper)^2))
l2.sub.s
dist.t.sub.s
beta.bar.sub.s2 = beta.bar.sub.s / x.sd.sub
sum(beta.bar.sub.s2) # not satisfied anymore
beta.bar.sub.s2
l2.sub.s2 = sqrt(crossprod(beta.bar.sub.s2))
dist.t.sub.s2 = sqrt(sum((beta.bar.sub.s2 - beta.paper)^2))
l2.sub.s2
dist.t.sub.s2

# constrained - with intercept
y.mean = mean(y)
y.c = y - y.mean
x.mean.sub = colMeans(Xmat.sub)
Xmat.sub.c = scale(Xmat.sub, center=x.mean.sub, scale=F)
beta.hat.sub.c = solve(crossprod(Xmat.sub.c), crossprod(Xmat.sub.c, y.c))
# check that it matches lm() #
beta.hat.lmint = coefficients(lm(y ~ Xmat.sub))
# coefficients when centering match when not-centering (just not intercept, since centering gets rid of intercept)
all.equal(as.numeric(coefficients(lm(y ~ Xmat.sub))[-1]), 
          as.numeric(coefficients(lm(y.c ~ Xmat.sub.c))[-1]))
all.equal(as.numeric(beta.hat.sub.c), as.numeric(beta.hat.lmint)[-1]) # they match, as expected
beta.hat.lmint0 = y.mean - as.vector(x.mean.sub%*%beta.hat.sub.c) # also match
# compute beta.bar, constrained fit
beta.bar.sub.c = beta.hat.sub.c - solve(crossprod(Xmat.sub.c), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.sub.c), Q)), crossprod(Q, beta.hat.sub.c))
sum(beta.bar.sub.c)
beta.bar.sub.c
beta.bar.sub.c0 <- y.mean - as.vector(x.mean.sub%*%beta.bar.sub.c)
# constrained - with intercept, including in Xmat.sub
Xmat.sub.c2 = cbind(1, Xmat.sub)
beta.hat.sub.c2 = solve(crossprod(Xmat.sub.c2), crossprod(Xmat.sub.c2, y))
Q2 = c(0, as.matrix(rep(1, dim(Xmat.sub)[2])))
beta.bar.sub.c2 = beta.hat.sub.c2 - solve(crossprod(Xmat.sub.c2), Q2) %*% 
  solve(crossprod(Q2, solve(crossprod(Xmat.sub.c2), Q2)), crossprod(Q2, beta.hat.sub.c2))
sum(beta.bar.sub.c2[-1])
beta.bar.sub.c2 
all.equal(c(beta.bar.sub.c0, beta.bar.sub.c), as.numeric(beta.bar.sub.c2)) # this matches the above
l2.sub.c = sqrt(crossprod(beta.bar.sub.c))
dist.t.sub.c = sqrt(sum((beta.bar.sub.c - beta.paper)^2))
l2.sub.c
dist.t.sub.c

# constrained - with intercept and scaling (standardizing, bc need to center for intercept)
Xmat.sub.cs = scale(Xmat.sub, center=x.mean.sub, scale=x.sd.sub)
beta.hat.sub.cs = solve(crossprod(Xmat.sub.cs), crossprod(Xmat.sub.cs, y.c))
beta.bar.sub.cs = beta.hat.sub.cs - solve(crossprod(Xmat.sub.cs), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat.sub.cs), Q)), crossprod(Q, beta.hat.sub.cs))
sum(beta.bar.sub.cs)
beta.bar.sub.cs
beta.bar.sub.cs0 <- y.mean - as.vector(x.mean.sub%*%beta.bar.sub.cs)
beta.bar.sub.cs0
l2.sub.cs = sqrt(crossprod(beta.bar.sub.cs))
dist.t.sub.cs = sqrt(sum((beta.bar.sub.cs - beta.paper)^2))
l2.sub.cs
dist.t.sub.cs
beta.bar.sub.cs2 = beta.bar.sub.cs / x.sd.sub
beta.bar.sub.cs2
sum(beta.bar.sub.cs2) # not satisfied anymore
beta.bar.sub.cs02 <- y.mean - as.vector(x.mean.sub%*%beta.bar.sub.cs2)
beta.bar.sub.cs02
sum(beta.bar.sub.cs02, beta.bar.sub.cs2)
sum(beta.bar.sub.cs2) # not satisfied anymore
l2.sub.cs2 = sqrt(crossprod(beta.bar.sub.cs2))
dist.t.sub.cs2 = sqrt(sum((beta.bar.sub.cs2 - beta.paper)^2))
l2.sub.cs2
dist.t.sub.cs2
