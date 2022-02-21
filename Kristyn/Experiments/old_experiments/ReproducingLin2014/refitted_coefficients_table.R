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

# testing with centered and/or scaled y
# y = scale(y, center = F, scale = F)

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
# fit lm with intercept, standardizing input : lm(y ~ -1 + log.X.prop.std)
################################################################################
log.X.prop.std = scale(log.X.prop, 
                          center = apply(log.X.prop, 2, mean), 
                          scale = apply(log.X.prop, 2, sd))
final.data.stdX = data.frame(log.X.prop.std[, final.selected], y = y)
final.lm.stdX = lm(y ~ ., final.data.stdX)
coefficients(final.lm.stdX)
beta.lm.stdX = coefficients(final.lm.stdX)[-1]

################################################################################
# fit lm withOUT intercept, standardizing input : lm(y ~ -1 + log.X.prop.scaled)
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


################################################################################
# Trying to take out genera that are too sparse #
################################################################################
# we don't want the proportion of zeroes for an otu to be too big
#   (i.e. otu is too sparse, not in many samples)
length.too.sparse = 25
too.sparse = seq(from = 0.866, to = 1, length = length.too.sparse)
prop_zero_per_sample = as.vector(apply(X, 1, function(a) sum(a == 0.5) / ncol(X)))
prop_zero_per_otu = as.vector(apply(X, 2, function(a) sum(a == 0.5) / nrow(X)))
# lm
lm.l2norms = rep(NA, length.too.sparse)
lm.dists = rep(NA, length.too.sparse)
lm.noint.l2norms = rep(NA, length.too.sparse)
lm.noint.dists = rep(NA, length.too.sparse)
# lm, scale
lm.scale.l2norms = rep(NA, length.too.sparse)
lm.scale.dists = rep(NA, length.too.sparse)
lm.noint.scale.l2norms = rep(NA, length.too.sparse)
lm.noint.scale.dists = rep(NA, length.too.sparse)
# lm, std
lm.std.l2norms = rep(NA, length.too.sparse)
lm.std.dists = rep(NA, length.too.sparse)
lm.noint.std.l2norms = rep(NA, length.too.sparse)
lm.noint.std.dists = rep(NA, length.too.sparse)
# start calculations
beta.paper = c(-0.76, -1.35, 0.61, 1.50)
for(i in 1:length(too.sparse)){
  prop.zero = too.sparse[i]
  print(paste0("### for genera with less than ", prop.zero, " percent zeroes, ... ###"))
  X.prune = X[, prop_zero_per_otu < prop.zero] # take out OTUs with less than (1 - prop.zero)% nonzero
  contains.final.selected = all(final.selected %in% colnames(X.prune))
  if(!contains.final.selected){
    print("skipping, since doesn't contain all 4 genera")
    next
  }
  X.prop.prune <- sweep(X.prune, MARGIN = 1, STATS = rowSums(X.prune), FUN = "/")
  log.X.prop.prune = log(X.prop.prune)
  ################################################################################
  # fit lm with intercept : lm(y ~ log.X.prop)
  ################################################################################
  data.i = data.frame(log.X.prop.prune[, final.selected], y)
  lm.i = lm(y ~ ., data.i)
  beta.i = coefficients(lm.i)[-1]
  lm.l2norms[i] = sqrt(crossprod(beta.i))
  lm.dists[i] = sqrt(sum((beta.i - beta.paper)^2))
  ################################################################################
  # fit lm withOUT intercept : lm(y ~ -1 + log.X.prop)
  ################################################################################
  lm.i = lm(y ~ -1 + ., data.i)
  beta.i = coefficients(lm.i)
  lm.noint.l2norms[i] = sqrt(crossprod(beta.i))
  lm.noint.dists[i] = sqrt(sum((beta.i - beta.paper)^2))
  ################################################################################
  # fit lm with intercept, scaling input : lm(y ~ -1 + log.X.prop.scaled)
  ################################################################################
  log.X.prop.scaled = scale(log.X.prop.prune, 
                            center = FALSE, 
                            scale = apply(log.X.prop.prune, 2, sd))
  data.i = data.frame(log.X.prop.scaled[, final.selected], y)
  lm.i = lm(y ~ ., data.i)
  beta.i = coefficients(lm.i)[-1]
  lm.scale.l2norms[i] = sqrt(crossprod(beta.i))
  lm.scale.dists[i] = sqrt(sum((beta.i - beta.paper)^2))
  ################################################################################
  # fit lm withOUT intercept, scaling input : lm(y ~ -1 + log.X.prop.scaled)
  ################################################################################
  lm.i = lm(y ~ -1 + ., data.i)
  beta.i = coefficients(lm.i)
  lm.noint.scale.l2norms[i] = sqrt(crossprod(beta.i))
  lm.noint.scale.dists[i] = sqrt(sum((beta.i - beta.paper)^2))
  ################################################################################
  # fit lm with intercept, standardizing input : lm(y.scaled ~ -1 + log.X.prop.std)
  ################################################################################
  log.X.prop.std = scale(log.X.prop.prune, 
                         center = apply(log.X.prop.prune, 2, mean), 
                         scale = apply(log.X.prop.prune, 2, sd))
  data.i = data.frame(log.X.prop.std[, final.selected], y = y)
  lm.i = lm(y ~ ., data.i)
  beta.i = coefficients(lm.i)[-1]
  lm.std.l2norms[i] = sqrt(crossprod(beta.i))
  lm.std.dists[i] = sqrt(sum((beta.i - beta.paper)^2))
  ################################################################################
  # fit lm withOUT intercept, standardizing input & output : lm(y.scaled ~ -1 + log.X.prop.scaled)
  ################################################################################
  lm.i = lm(y ~ -1 + ., data.i)
  beta.i = coefficients(lm.i)
  lm.noint.std.l2norms[i] = sqrt(crossprod(beta.i))
  lm.noint.std.dists[i] = sqrt(sum((beta.i - beta.paper)^2))
}
# plot results
library(reshape2)
l2data = data.frame(
  x = too.sparse,
  lm = lm.l2norms, 
  noint = lm.noint.l2norms, 
  scale = lm.scale.l2norms, 
  noint.scale = lm.noint.scale.l2norms, 
  std = lm.std.l2norms, 
  noint.std = lm.noint.std.l2norms
)
l2data.melt = melt(
  l2data, 
  measure.vars = c("lm", "noint", "scale", "noint.scale", "std", "noint.std"), 
  variable.name = "fit")
ggplot(data = l2data.melt, aes(x = x, y = value)) + 
  geom_line() + 
  facet_wrap(vars(fit), scales = "free")
ggsave("l2norms_112320.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in")
)
dist.data = data.frame(
  x = too.sparse,
  lm = lm.dists, 
  noint = lm.noint.dists, 
  scale = lm.scale.dists, 
  noint.scale = lm.noint.scale.dists, 
  std = lm.std.dists, 
  noint.std = lm.noint.std.dists
)
dist.data.melt = melt(
  dist.data, 
  measure.vars = c("lm", "noint", "scale", "noint.scale", "std", "noint.std"), 
  variable.name = "fit")
ggplot(data = dist.data.melt, aes(x = x, y = value)) + 
  geom_line() + 
  facet_wrap(vars(fit), scales = "free")
ggsave("dists_112320.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       
       scale = 1,
       width = 6,
       height = 4,
       units = c("in")
)
