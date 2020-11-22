# workdir = "/home/kristyn/Documents/research/supervisedlogratios/LogRatioReg"
# setwd(workdir)

# libraries
library(mvtnorm)
library(balance)
library(selbal)
library(propr)
library(microbenchmark)
library(ggplot2)
library(logratiolasso) # bates & tibshirani 2019
image_path = "/home/kristyn/Pictures"

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
source(paste0(functions_path, "propr.R"))
source(paste0(functions_path, "selbal.R"))

# data
# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
DataFolder <- "/Data/"
load(paste0(workdir, DataFolder, "BMI.rda"))
log.X.prop = log(X.prop)
# dim(raw_data) # 98 x 89
# dim(X) # 98 x 87
# dim(X.prop) # 98 x 87
n = dim(X)[1]
num.genera = dim(X)[2]

# testing with centered and/or scaled y
y = scale(y, center = F, scale = F)

final.selected = c(
  "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes",
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella")

################################################################################
# fit lm with intercept : lm(y ~ log.X.prop)
################################################################################
final.data = data.frame(log.X.prop[, final.selected], y)
final.lm = lm(y ~ ., final.data)
coefficients(final.lm)
beta.lm = coefficients(final.lm)[-1]

################################################################################
# fit lm withOUT intercept : lm(y ~ -1 + log.X.prop)
################################################################################
final.lm.noint = lm(y ~ -1 + ., final.data)
coefficients(final.lm.noint)
beta.lm.noint = coefficients(final.lm.noint)

################################################################################
# fit lm with intercept, scaling input : lm(y ~ -1 + log.X.prop.scaled)
################################################################################
log.X.prop.scaled = scale(log.X.prop, 
                          center = FALSE, 
                          scale = apply(log.X.prop, 2, sd))
final.data.scaledX = data.frame(log.X.prop.scaled[, final.selected], y)
final.lm.scaledX = lm(y ~ ., final.data.scaledX)
coefficients(final.lm.scaledX)
beta.lm.scaledX = coefficients(final.lm.scaledX)[-1]

################################################################################
# fit lm withOUT intercept, scaling input : lm(y ~ -1 + log.X.prop.scaled)
################################################################################
final.lm.noint.scaledX = lm(y ~ -1 + ., final.data.scaledX)
coefficients(final.lm.noint.scaledX)
beta.lm.noint.scaledX = coefficients(final.lm.noint.scaledX)

################################################################################
# fit lm with intercept, standardizing input : lm(y.scaled ~ -1 + log.X.prop.std)
################################################################################
log.X.prop.std = scale(log.X.prop, 
                          center = apply(log.X.prop, 2, mean), 
                          scale = apply(log.X.prop, 2, sd))
final.data.stdX = data.frame(log.X.prop.std[, final.selected], y = y)
final.lm.stdX = lm(y ~ ., final.data.stdX)
coefficients(final.lm.stdX)
beta.lm.stdX = coefficients(final.lm.stdX)[-1]

################################################################################
# fit lm withOUT intercept, standardizing input & output : lm(y.scaled ~ -1 + log.X.prop.scaled)
################################################################################
final.lm.noint.stdX = lm(y ~ -1 + ., final.data.stdX)
coefficients(final.lm.noint.stdX)
beta.lm.noint.stdX = coefficients(final.lm.noint.stdX)

################################################################################
################################################################################
################################################################################
# " We are interested in identifying a subset of important genera whose 
#   *subcomposition* is associated with BMI."
## trying the composition that includes only these genera...
##  didn't work. Property 3 doesn't refer to a subcomposition in this way, anyway.
# X.sel = X[, final.selected]
# X.prop.sel = sweep(X.sel, MARGIN = 1, STATS = rowSums(X.sel), FUN = "/")
# log.X.prop.sel = log(X.prop.sel)
# final.data.sel = data.frame(log.X.prop.sel, y)
# #
# final.lm.sel = lm(y ~ ., final.data.sel)
# coefficients(final.lm.sel)
# #
# log.X.prop.sel.scaled = scale(log.X.prop.sel,
#                               center = FALSE, #apply(log.X.prop.sel, 2, mean),
#                               scale = apply(log.X.prop.sel, 2, sd))
# final.data.sel.scaledX = data.frame(log.X.prop.sel.scaled, y)
# #
# final.lm.sel.scaledX = lm(y ~ ., final.data.sel.scaledX)
# coefficients(final.lm.sel.scaledX)

################################################################################
################################################################################
################################################################################
beta.paper = c(-0.76, -1.35, 0.61, 1.50)
betas = data.frame(
  paper = beta.paper, 
  lm = beta.lm, 
  lm.noint = beta.lm.noint, 
  lm.scaledX = beta.lm.scaledX, 
  lm.noint.scaledX = beta.lm.noint.scaledX, 
  lm.stdX = beta.lm.stdX, 
  lm.noint.stdX = beta.lm.noint.stdX
)
round(betas, 2)
l2norm = rep(NA, ncol(betas))
names(l2norm) = colnames(betas)
dist.paper = rep(NA, ncol(betas))
names(dist.paper) = colnames(betas)
for(i in 1:ncol(betas)){
  beta.i = betas[, i]
  l2norm[i] = sqrt(crossprod(beta.i))
  dist.paper[i] = sqrt(sum((beta.i - beta.paper)^2))
}
l2norm
dist.paper
