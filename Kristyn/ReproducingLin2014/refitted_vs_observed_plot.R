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
source(paste0(functions_path, "constrainedlm.R"))

# helper functions
MSEyhat <- function(y, yhat){
  sqrt(sum((y - yhat)^2)) / length(y)
}

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

# testing with centered and/or scaled y
# y = scale(y, center = T, scale = T)

# try removing the 5 obese observations
sort(y)
sort(scale(y))
# they have BMI (y) > 36
# include = y < 36
# # checking if rows still sum to 1
# apply(X.prop, MARGIN = 1, sum)
# X = X[include, ]
# X.prop = X.prop[include, ]
# apply(X.prop, MARGIN = 1, sum)
# log.X.prop = log.X.prop[include, ]
# y = y[include]

final.selected = c(
  "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes",
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella")
Xmat = log.X.prop[, final.selected]

betahat.paper = c(-0.76, -1.35, 0.61, 1.50)
betahat0.paper = as.numeric(mean(y) - betahat.paper %*% colMeans(Xmat))

################################################################################
# REPRODUCED ORIGINAL FITTED vs. OBSERVED PLOT
################################################################################
par(mfrow = c(1, 1))
yhat.paper = betahat0.paper + Xmat %*% betahat.paper
xyrange = range(y, yhat.paper)
plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
     x = y, y = yhat.paper, col = 1,
     xlim = xyrange, ylim = xyrange)
abline(0, 1, lty = 2)
# yay, looks just like it!
# now, to figure out how betahat.paper was obtained...

# check standardize() and backStandardize()
stdXY = standardize(Xmat, y, center = FALSE, scale = FALSE)
backstdXY = backStandardize(stdXY, betahat.paper, scale = FALSE) # matches paper
backstdXY
betahat0.paper
betahat.paper
# sum(betahat)
sum(betahat.paper)
# L2
sqrt(crossprod(betahat.paper))
# Euclidean distance from paper's beta
sqrt(sum((betahat.paper - betahat.paper)^2))
# calculate MSE(yhat)
MSEyhat(y, yhat.paper)

################################################################################
# Constraint Linear Model to Original Data
################################################################################
par(mfrow = c(2, 2))
# constrained - no centering or scaling: 

# check clm()
# calculate by hand
# beta-bar = beta-hat - (X'X)^(-1) 1 [1' (X'X)^(-1) 1]^(-1) (1' beta-hat)
beta.hat = solve(crossprod(Xmat), crossprod(Xmat, y))
# check that it matches lm()
beta.hat.lm = coefficients(lm(y ~ -1 + Xmat))
all.equal(as.numeric(beta.hat), as.numeric(beta.hat.lm)) # they match, as expected
# compute beta.bar, constrained fit
Q = as.matrix(rep(1, dim(Xmat)[2]))
beta.bar = beta.hat - solve(crossprod(Xmat), Q) %*% 
  solve(crossprod(Q, solve(crossprod(Xmat), Q)), crossprod(Q, beta.hat))
# check function output
beta.bar2 = clm(Xmat, y, Q)
all.equal(beta.bar, beta.bar2) # they match :D

# fit to the original data
Q = as.matrix(rep(1, dim(Xmat)[2]))
stdXY = standardize(Xmat, y, center = FALSE, scale = FALSE)
# fit betabar
betabar.orig = clm(Xmat, y, Q)
# we can ignore the warning,  we just need this as input to backStandardize()
# get betabar0
betabar.bstd.orig = backStandardize(stdXY, betabar.orig, scale = FALSE)
# calculat yhat from fit
yhat.orig = betabar.bstd.orig$betahat0 + Xmat %*% betabar.bstd.orig$betahat
# plot observed vs fitted values
plot(main = "No centering/scaling", xlab = "Observed BMI", ylab = "Fitted BMI",
     x = y, y = yhat.paper, col = 2,
     xlim = xyrange, ylim = xyrange)
points(x = y, y = yhat.orig)
abline(0, 1, lty = 2)
# sum(betahat)
sum(betabar.bstd.orig$betahat)
# L2
sqrt(crossprod(betabar.bstd.orig$betahat))
# Euclidean distance from paper's beta
sqrt(sum((betabar.bstd.orig$betahat - betahat.paper)^2))
# calculate MSE(yhat)
MSEyhat(y, yhat.orig)

################################################################################
# Scale
################################################################################

# transform Xmat by doing scale = TRUE
stdXY = standardize(Xmat, y, center = FALSE, scale = TRUE)
all.equal(y, stdXY$Ytilde)
all.equal(Xmat, stdXY$Xtilde)
# fit betabar
betabar = clm(stdXY$Xtilde, stdXY$Ytilde, Q)
# we can ignore the warning,  we just need this as input to backStandardize()

# # # wrong input!
# fitted values with 
#     scale = FALSE
#     scaled Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + stdXY$Xtilde %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

# fitted values with 
#     scale = FALSE
#     original Xmat
# get betabar0
betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
# calculate yhat, 
yhat = betabar.bstd$betahat0 + Xmat %*% betabar.bstd$betahat
# plot observed vs fitted values
plot(main = "Scaling", xlab = "Observed BMI", ylab = "Fitted BMI",
     x = y, y = yhat.paper, col = 2,
     xlim = xyrange, ylim = xyrange)
points(x = y, y = yhat)
abline(0, 1, lty = 2)
# sum(betahat)
sum(betabar.bstd$betahat)
# L2
sqrt(crossprod(betabar.bstd$betahat))
# Euclidean distance from paper's beta
sqrt(sum((betabar.bstd$betahat - betahat.paper)^2))
# calculate MSE(yhat)
MSEyhat(y, yhat)

# # # wrong input!
# fitted values with 
#     scale = TRUE
#     scaled Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = TRUE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + stdXY$Xtilde %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

# # # doesn't satisfy constraint if we back-scale!
# fitted values with 
#     scale = TRUE
#     original Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = TRUE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + Xmat %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

################################################################################
# Center
################################################################################

# transform Xmat by doing center = TRUE
stdXY = standardize(Xmat, y, center = TRUE, scale = FALSE)
all.equal(y, stdXY$Ytilde)
all.equal(Xmat, stdXY$Xtilde)
# fit betabar
betabar = clm(stdXY$Xtilde, stdXY$Ytilde, Q)
# we can ignore the warning,  we just need this as input to backStandardize()

# # # wrong input!
# fitted values with 
#     scale = FALSE
#     centered Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + stdXY$Xtilde %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "Centering", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

# fitted values with 
#     scale = FALSE
#     original Xmat
# get betabar0
betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
# calculate yhat, 
yhat = betabar.bstd$betahat0 + Xmat %*% betabar.bstd$betahat
# plot observed vs fitted values
plot(main = "Centering", xlab = "Observed BMI", ylab = "Fitted BMI",
     x = y, y = yhat.paper, col = 2,
     xlim = xyrange, ylim = xyrange)
points(x = y, y = yhat)
abline(0, 1, lty = 2)
# sum(betahat)
sum(betabar.bstd$betahat)
# L2
sqrt(crossprod(betabar.bstd$betahat))
# Euclidean distance from paper's beta
sqrt(sum((betabar.bstd$betahat - betahat.paper)^2))
# calculate MSE(yhat)
MSEyhat(y, yhat)

# # # wrong input!
# fitted values with
#     scale = TRUE
#     scaled Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = TRUE)
# # calculate yhat,
# yhat = betabar.bstd$betahat0 + stdXY$Xtilde %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

# # # doesn't satisfy constraint if we back-scale!
# # # besides, we didn't even scale, to begin with.
##### the following is strangely close...
# fitted values with
#     scale = TRUE
#     original Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = TRUE)
# # calculate yhat,
# yhat = betabar.bstd$betahat0 + Xmat %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

################################################################################
# Center & Scale
################################################################################

# transform Xmat by doing scale = TRUE
stdXY = standardize(Xmat, y, center = TRUE, scale = TRUE)
all.equal(y, stdXY$Ytilde)
all.equal(Xmat, stdXY$Xtilde)
# fit betabar
betabar = clm(stdXY$Xtilde, stdXY$Ytilde, Q)
# we can ignore the warning,  we just need this as input to backStandardize()

# # # wrong input!
# fitted values with 
#     scale = FALSE
#     centered & scaled Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + stdXY$Xtilde %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

# fitted values with 
#     scale = FALSE
#     original Xmat
# get betabar0
betabar.bstd = backStandardize(stdXY, betabar, scale = FALSE)
# calculate yhat, 
yhat = betabar.bstd$betahat0 + Xmat %*% betabar.bstd$betahat
# plot observed vs fitted values
plot(main = "Centering & Scaling", xlab = "Observed BMI", ylab = "Fitted BMI",
     x = y, y = yhat.paper, col = 2,
     xlim = xyrange, ylim = xyrange)
points(x = y, y = yhat)
abline(0, 1, lty = 2)
# sum(betahat)
sum(betabar.bstd$betahat)
# L2
sqrt(crossprod(betabar.bstd$betahat))
# Euclidean distance from paper's beta
sqrt(sum((betabar.bstd$betahat - betahat.paper)^2))
# calculate MSE(yhat)
MSEyhat(y, yhat)

# # # wrong input!
# fitted values with 
#     scale = TRUE
#     scaled Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = TRUE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + stdXY$Xtilde %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)

# doesn't satisfy constraint if we back-scale!
# fitted values with 
#     scale = TRUE
#     original Xmat
# # get betabar0
# betabar.bstd = backStandardize(stdXY, betabar, scale = TRUE)
# # calculate yhat, 
# yhat = betabar.bstd$betahat0 + Xmat %*% betabar.bstd$betahat
# # plot observed vs fitted values
# plot(main = "", xlab = "Observed BMI", ylab = "Fitted BMI",
#      x = y, y = yhat.paper, col = 2,
#      xlim = xyrange, ylim = xyrange)
# points(x = y, y = yhat)
# abline(0, 1, lty = 2)
# # calculate MSE(yhat)
# MSEyhat(y, yhat)
# # sum(betahat)
# sum(betabar.bstd$betahat)






# 
# 
# 
# 
# 
# 
# 
# 
# 
# # try it on the original data, to see if we get beta.paper
# Xmat = log.X.prop[, final.selected]
# beta.hat = solve(crossprod(Xmat), crossprod(Xmat, y))
# Q = as.matrix(rep(1, dim(Xmat)[2]))
# beta.bar = beta.hat - solve(crossprod(Xmat), Q) %*% 
#   solve(crossprod(Q, solve(crossprod(Xmat), Q)), crossprod(Q, beta.hat))
# stdXY = standardize(Xmat, y, center = FALSE, scale = FALSE)
# backstdXY = backStandardize(stdXY, beta.bar, scale = FALSE); backstdXY
# # try scaling
# stdXY.s = standardize(Xmat, y, center = FALSE, scale = TRUE)
# beta.hat.s = solve(crossprod(stdXY.s$Xtilde), crossprod(stdXY.s$Xtilde, stdXY.s$Ytilde))
# beta.bar.s = beta.hat.s - solve(crossprod(stdXY.s$Xtilde), Q) %*% 
#   solve(crossprod(Q, solve(crossprod(stdXY.s$Xtilde), Q)), crossprod(Q, beta.hat.s))
# backstdXY.s = backStandardize(stdXY.s, beta.bar.s, scale = FALSE); backstdXY.s
# backstdXY.s2 = backStandardize(stdXY.s, beta.bar.s, scale = TRUE); backstdXY.s2
# # try centering
# stdXY.c = standardize(Xmat, y, center = TRUE, scale = FALSE)
# beta.hat.c = solve(crossprod(stdXY.c$Xtilde), crossprod(stdXY.c$Xtilde, stdXY.c$Ytilde))
# beta.bar.c = beta.hat.c - solve(crossprod(stdXY.c$Xtilde), Q) %*% 
#   solve(crossprod(Q, solve(crossprod(stdXY.c$Xtilde), Q)), crossprod(Q, beta.hat.c))
# backstdXY.c = backStandardize(stdXY.c, beta.bar.c, scale = FALSE); backstdXY.c
# backstdXY.c2 = backStandardize(stdXY.c, beta.bar.c, scale = TRUE); backstdXY.c2
# # try centering and scaling
# stdXY.cs = standardize(Xmat, y, center = TRUE, scale = FALSE)
# beta.hat.cs = solve(crossprod(stdXY.cs$Xtilde), crossprod(stdXY.cs$Xtilde, stdXY.cs$Ytilde))
# beta.bar.cs = beta.hat.cs - solve(crossprod(stdXY.cs$Xtilde), Q) %*% 
#   solve(crossprod(Q, solve(crossprod(stdXY.cs$Xtilde), Q)), crossprod(Q, beta.hat.cs))
# backstdXY.cs = backStandardize(stdXY.cs, beta.bar.cs, scale = FALSE); backstdXY.cs
# backstdXY.cs2 = backStandardize(stdXY.cs, beta.bar.cs, scale = TRUE); backstdXY.cs2
# 
