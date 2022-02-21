################################################################################
# no refitting, parallelized
################################################################################

pe0 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/prediction_performance_parallelized.rds"))
dim(pe0)
mean.pe0 = colMeans(pe0)
sqrt.mean.pe0 = sqrt(mean.pe0)
mean(mean.pe0) # 28.36351
sd(mean.pe0) # 10.2603
# standard error?
sd.pe0 = apply(pe0, 2, sd)
mean(sd.pe0) / sqrt(dim(pe0)[1]) # 10.57

################################################################################
# refitting, parallelized
################################################################################

pe1 = readRDS(paste0(
  "Kristyn/ReproducingLin2014/bootstrap_simulations_121420",
  "/prediction_performance_refitted_parallelized.rds"))
dim(pe1)
mean.pe1 = colMeans(pe1)
mean(mean.pe1) # 33.85995
sd(mean.pe1) # 13.73733
# standard error?
sd.pe1 = apply(pe1, 2, sd)
mean(sd.pe1) / sqrt(dim(pe1)[1]) # 12.04
