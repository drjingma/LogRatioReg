# constrained least squares fit
# last update: 12/20/20

################################################################################
# sources, libraries, helper functions
################################################################################

# libraries
library(mvtnorm)
library(balance)
library(selbal)
library(microbenchmark)
library(ggplot2)

# constrained least squares
library(limSolve)

# ones from github repos
library(selbal)
library(logratiolasso) # bates & tibshirani 2019

# Dr. Ma sources
source("RCode/func_libs.R")
source("COAT-master/coat.R")

# Kristyn sources
functions_path = "Kristyn/Functions/"
source(paste0(functions_path, "constrainedlm.R"))

# helper functions
MSEyhat <- function(y, yhat){
  sqrt(sum((y - yhat)^2)) / length(y)
}

################################################################################
# data
################################################################################

# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
DataFolder <- "/Data/"
load(paste0("Data/", "BMI.rda"))
log.X.prop = log(X.prop)
n = dim(X)[1]
num.genera = dim(X)[2]

# looking at log.X.prop
all.equal(diag(1/rowSums(X)), solve(diag(rowSums(X))))
all.equal(X.prop, diag(1/rowSums(X)) %*% X)
all.equal(log.X.prop, log(X) - log(diag(rowSums(X))))
V = diag(rowSums(X))
Vinv = diag(1/rowSums(X))
W = X
all.equal(
  crossprod(log.X.prop), 
  crossprod(log(solve(V, W)))
)
all.equal(
  crossprod(log.X.prop), 
  log(t(W) %*% Vinv) %*% log(Vinv %*% W)
)
all.equal(
  crossprod(log.X.prop), 
  log(t(W) %*% Vinv %*% Vinv %*% W)
)

################################################################################
# fitting to genera selected in Lin et. al. 2014
################################################################################

final.selected = c(
  "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes",
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella")
Xmat = log.X.prop[, final.selected]

betahat.paper = c(-0.76, -1.35, 0.61, 1.50)
betahat0.paper = as.numeric(mean(y) - betahat.paper %*% colMeans(Xmat))

################################################################################
# reproduced fitted vs observed plot
################################################################################

yhat.paper = betahat0.paper + Xmat %*% betahat.paper
xyrange = range(y, yhat.paper)
plot(main = "Using original composition", xlab = "Observed BMI", ylab = "Fitted BMI",
     x = y, y = yhat.paper, col = 1, xlim = xyrange, ylim = xyrange)
abline(0, 1, lty = 2)
# yay, looks just like it!

################################################################################
# get the subcomposition
################################################################################

X.sub = X[, final.selected]
X.prop.sub = sweep(X.sub, MARGIN = 1, STATS = rowSums(X.sub), FUN = "/")
Xmat.sub = log(X.prop.sub)
all.equal(Xmat, Xmat.sub) # "Mean relative difference: 0.4333185"
all.equal(crossprod(Xmat), crossprod(Xmat.sub)) # different
all.equal(solve(crossprod(Xmat), crossprod(Xmat, y)), 
          solve(crossprod(Xmat.sub), crossprod(Xmat.sub, y))) # different

################################################################################
# other: divide by row max
################################################################################

X.prop.other = sweep(X.sub, MARGIN = 1, STATS = rowMaxs(X), FUN = "/")
Xmat.other = log(X.prop.other)
all.equal(Xmat, Xmat.other) # "Mean relative difference: 0.4333185"

################################################################################
# another: divide by row max
################################################################################

X.prop.ano = sweep(X.sub, MARGIN = 1, STATS = rep(max(X), n), FUN = "/")
Xmat.ano = log(X.prop.ano)
all.equal(Xmat, Xmat.ano) # "Mean relative difference: 0.4333185"

################################################################################
# compare betahat0 using original and sub- composition (using Lin et. al. est.)
################################################################################

# original composition
betahat0.paper = as.numeric(mean(y) - betahat.paper %*% colMeans(Xmat)) # recall
# subcomposition
betahat0.paper.sub = as.numeric(mean(y) - betahat.paper %*% colMeans(Xmat.sub))
betahat0.paper.other = as.numeric(mean(y) - betahat.paper %*% colMeans(Xmat.other)) # equal
as.numeric(mean(y) - betahat.paper %*% colMeans(log(X[, final.selected])))
# they are equal
betahat0.paper == betahat0.paper.sub # TRUE
print(data.frame(
  original = betahat0.paper, 
  subcomposition = betahat0.paper.sub))
# at what point did they become equal?
all.equal(colMeans(Xmat), colMeans(Xmat.sub)) # "Mean relative difference: 0.4333185"
all.equal(
  betahat.paper %*% colMeans(Xmat), 
  betahat.paper %*% colMeans(Xmat.sub)) # TRUE

################################################################################
# compare betahat using original and sub- composition
################################################################################
stdXY = standardize(Xmat, y, center = TRUE, scale = FALSE)
stdXY.sub = standardize(Xmat.sub, y, center = TRUE, scale = FALSE)
y.cen = y - mean(y)

# get betahat using lm()
# original composition
Xmat.cen = stdXY$Xtilde
betahat = coefficients(lm(y.cen ~ -1 + Xmat.cen))
names(betahat) = colnames(Xmat)
# subcomposition
Xmat.sub.cen = stdXY.sub$Xtilde
betahat.sub = coefficients(lm(y.cen ~ -1 + Xmat.sub.cen))
names(betahat.sub) = colnames(Xmat.sub.cen)
# they are not equal:
all.equal(betahat, betahat.sub) # "Mean relative difference: 0.5880189"
print(data.frame(
  original = betahat, 
  subcomposition = betahat.sub))

# solve for betahat by hand
# original
betahat2 = solve(crossprod(Xmat.cen), crossprod(Xmat.cen, y.cen))
all.equal(as.vector(betahat), as.vector(betahat2)) # by hand matches lm()
# using subcomposition
betahat.sub2 = solve(crossprod(Xmat.sub.cen), crossprod(Xmat.sub.cen, y.cen))
all.equal(as.vector(betahat.sub), as.vector(betahat.sub2)) # by hand matches lm()
# they are not equal
all.equal(betahat2, betahat.sub2) # "Mean relative difference: 0.5880189"
print(data.frame(
  original = betahat2, 
  subcomposition = betahat.sub2))

# betahat using original vs sub- composition are NOT equal

################################################################################
# compare betabar using original and sub- composition
################################################################################

# get betabar using lsei()
# original
Xmat.int = cbind(1, Xmat)
Q = as.matrix(rep(1, dim(Xmat)[2]))
Q.int = rbind(0, Q)
betabar = lsei(A = Xmat.int, B = y, E = t(Q.int), F = 0)$X
# subcomposition
Xmat.int.sub = cbind(1, Xmat.sub)
Q = as.matrix(rep(1, dim(Xmat.sub)[2]))
Q.int = rbind(0, Q)
betabar.sub = lsei(A = Xmat.int.sub, B = y, E = t(Q.int), F = 0)$X
# other
Xmat.int.other = cbind(1, Xmat.other)
Q = as.matrix(rep(1, dim(Xmat.other)[2]))
Q.int = rbind(0, Q)
betabar.other = lsei(A = Xmat.int.other, B = y, E = t(Q.int), F = 0)$X
# another
Xmat.int.ano = cbind(1, Xmat.ano)
Q = as.matrix(rep(1, dim(Xmat.ano)[2]))
Q.int = rbind(0, Q)
betabar.ano = lsei(A = Xmat.int.ano, B = y, E = t(Q.int), F = 0)$X
# they are equal
all.equal(betabar, betabar.sub, betabar.other, betabar.ano) # TRUE
print(data.frame(
  original = betabar, 
  subcomposition = betabar.sub, 
  maxrow = betabar.other, 
  maxall = betabar.ano))

# solve for betabar by hand
# original
XtXinvQ = solve(crossprod(Xmat.cen), Q)
betabar2 = betahat2 - XtXinvQ %*% 
  solve(crossprod(Q, XtXinvQ), crossprod(Q, betahat2))
betahat0.2 = as.numeric(mean(y) - betabar2 %*% colMeans(Xmat))
betabar2 = backStandardize(stdXY, betabar2)
all.equal(as.numeric(unlist(betabar2)), as.vector(betabar)) # matches lsei()
# subcomposition
XtXinvQ.sub = solve(crossprod(Xmat.sub.cen), Q)
betabar.sub2 = betahat.sub2 - XtXinvQ.sub %*% 
  solve(crossprod(Q, XtXinvQ.sub), crossprod(Q, betahat.sub2))
betahat0.2 = as.numeric(mean(y) - betabar.sub2 %*% colMeans(Xmat.sub))
betabar.sub2 = backStandardize(stdXY.sub, betabar.sub2)
all.equal(as.numeric(unlist(betabar.sub2)), as.vector(betabar.sub)) # matches lsei()
# they are equal
all.equal(unlist(betabar2), unlist(betabar.sub2)) # TRUE
print(data.frame(
  original = unlist(betabar2), 
  subcomposition = unlist(betabar.sub2)))

# where did they start being equal?
# the OLS estimates are not equal
all.equal(betahat, betahat.sub) # "Mean relative difference: 0.5880189"
# the intercept parameter estimates are equal
all.equal(
  XtXinvQ %*% 
    solve(crossprod(Q, XtXinvQ), crossprod(Q, betahat2)), 
  XtXinvQ.sub %*% 
    solve(crossprod(Q, XtXinvQ.sub), crossprod(Q, betahat.sub2))
) # "Mean relative difference: 2.138686"
# crossprod(X)
all.equal(crossprod(Xmat.cen), crossprod(Xmat.sub.cen)) # "Mean relative difference: 1.898389"
### Rewrite betabar
# \overline{\beta} = (X^\intercal X)^{-1} \{ X^\intercal y -  
#   1_p [1_p^\intercal (X^\intercal X)^{-1} 1_p]^{-1} 1_p^\intercal \hat{\beta}\}
p = 4
# original
betabar3 = solve(
  crossprod(Xmat.cen), 
  crossprod(Xmat.cen, y.cen) - 
    Q %*% solve(crossprod(Q, XtXinvQ), crossprod(Q, betahat2)))
# check that it's actually betabar
all.equal(betabar2$betahat, betabar3) # TRUE
# subcomposition
betabar.sub3 = solve(
  crossprod(Xmat.sub.cen), 
  crossprod(Xmat.sub.cen, y.cen) - 
    Q %*% solve(crossprod(Q, XtXinvQ.sub), crossprod(Q, betahat.sub2)))
# check that it's actually betabar
all.equal(betabar.sub2$betahat, betabar.sub3) # TRUE
# check 1st product term
all.equal(
  solve(crossprod(Xmat.cen)), 
  solve(crossprod(Xmat.sub.cen))
) # "Mean relative difference: 1.856727"
# check 2nd product term
all.equal(
  crossprod(Xmat.cen, y.cen) - 
    Q %*% solve(crossprod(Q, XtXinvQ), crossprod(Q, betahat2)), 
  crossprod(Xmat.sub.cen, y.cen) - 
    Q %*% solve(crossprod(Q, XtXinvQ.sub), crossprod(Q, betahat.sub2))
) # "Mean relative difference: 1.856727"
all.equal(
  crossprod(Xmat.cen, y.cen), 
  crossprod(Xmat.sub.cen, y.cen)) # "Mean relative difference: 0.2600772"
### Rewrite betabar
# \overline{\beta} = (X^\intercal X)^{-1} \[ \{ I_p -  
#   1_p [1_p^\intercal (X^\intercal X)^{-1} 1_p]^{-1} 1_p^\intercal 
#     (X^\intercal X)^{-1} \} X\intercal y \]
# original
betabar4 = solve(
  crossprod(Xmat.cen), 
  (diag(p) - Q %*% solve(crossprod(Q, XtXinvQ)) %*% 
    crossprod(Q, solve(crossprod(Xmat.cen)))
   ) %*% crossprod(Xmat.cen, y.cen)
)
# check that it's actually betabar
all.equal(betabar2$betahat, betabar4) # TRUE
# subcomposition
betabar.sub4 = solve(
  crossprod(Xmat.sub.cen), 
  (diag(p) - Q %*% solve(crossprod(Q, XtXinvQ.sub)) %*% 
     crossprod(Q, solve(crossprod(Xmat.sub.cen)))
  ) %*% crossprod(Xmat.sub.cen, y.cen)
)
# check that it's actually betabar
all.equal(betabar.sub2$betahat, betabar.sub4) # TRUE
# checking things
all.equal(
  Q %*% solve(crossprod(Q, XtXinvQ)) %*% 
    crossprod(Q, solve(crossprod(Xmat.cen))), 
  Q %*% solve(crossprod(Q, XtXinvQ.sub)) %*% 
    crossprod(Q, solve(crossprod(Xmat.sub.cen)))
) # "Mean relative difference: 1.360309"
all.equal(
  Q %*% solve(crossprod(Q, XtXinvQ)) %*% t(Q), 
  Q %*% solve(crossprod(Q, XtXinvQ.sub)) %*% t(Q)
) # "Mean relative difference: 0.628727"
all.equal(
  solve(crossprod(Q, XtXinvQ)), 
  solve(crossprod(Q, XtXinvQ.sub))
) # "Mean relative difference: 0.628727"

