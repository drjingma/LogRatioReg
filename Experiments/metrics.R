# Purpose: compare methods on simulated data sets
# Date: 03/20/2025

# libraries and settings

output_dir = "../Outputs/data_final/"
source("Functions/util.R")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

# detach(package:plyr)

# tuning parameter settings ####
hparam = "1se"
replicates <- 1:50
run_selbal <- T
dataset <- 'Crohn_UC'
dim <- 50
PC <- 1

if (PC == 1){
  condition <- paste0(dataset, "_highCor_sparsity1_ratio1_dim", dim)
} else condition <- paste0(dataset, "_highCor_PC", PC, "_sparsity1_ratio1_dim", dim)

data.list <- readRDS(paste0('../Data/',condition,'.rds'))
truth <- rownames(data.list$sbp)[data.list$sbp!=0]
p <- ncol(data.list$X)

file.end = paste0(
  "/",condition,
  "_hparam", hparam)

cat(file.end,'\n')

metrics <- list()

# summary metrics ----

## codacore - 1 balance #########################################################
type <- 'codacore'
list.metrics = lapply(replicates, function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,"_", type, 
    "1_metrics",
    ".rds")))
type <- 'CoDaCoRe'
metrics[[type]] <- data.frame(Reduce(rbind, list.metrics))
rownames(metrics[[type]]) <- paste0('run',replicates)

## codaLasso #######################################################################
type <- 'classo'
list.metrics = lapply(replicates, function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid, "_", type, 
    "_metrics",
    ".rds")))
type <- 'codaLasso'
metrics[[type]] <- data.frame(Reduce(rbind,list.metrics))
rownames(metrics[[type]]) <- paste0('run',replicates)

## selbal #######################################################################
if (run_selbal){
  type <- 'selbal'
  list.metrics = lapply(replicates, function(jid) readRDS(
    paste0(
      output_dir, file.end, "_", jid, "_", type, 
      "_metrics",
      ".rds")))
  metrics[[type]] <- data.frame(Reduce(rbind,list.metrics))
  rownames(metrics[[type]]) <- paste0('run',replicates)
}

## slr ----
type <- 'slr_constrainedPC'
list.metrics = lapply(replicates, function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,"_", type, 
    "_metrics",
    ".rds")))
type = 'SLR'
metrics[[type]] <- data.frame(Reduce(rbind, list.metrics))
rownames(metrics[[type]]) <- paste0('run',replicates)

# combine all metrics ----
test <- lapply(1:length(metrics), function(a) data.frame(Method=names(metrics)[a],metrics[[a]]))
df <- Reduce(rbind,test) %>%
  dplyr::filter(Method %in% c('codaLasso','CoDaCoRe','selbal','SLR')) %>%
  mutate(Method = factor(Method))
head(df)


# Variable selection heat map ----
M.list <- list()
heatmap.list <- list()

## codacore - 1 balance #########################################################
type <- "codacore" 
metrics <- lapply(replicates, function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,
    "_", type, 
    "1.rds")))
M <- Reduce(cbind,lapply(metrics, function(jid)
  jid$ensemble[[1]]$hard$numerator - jid$ensemble[[1]]$hard$denominator
))
colnames(M) <- replicates
rownames(M) <- colnames(metrics[[1]]$x)
M <- M[rowSums(abs(M))>0,]

M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
type <- "CoDaCoRe"
p.temp <- pheatmap::pheatmap(M,
                             main=type,
                             cluster_rows = F,cluster_cols = F, silent = FALSE)
grid.newpage()
heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)
M.list[[type]] <- M


## codaLasso #######################################################################
type <- 'classo'
metrics = lapply(replicates, function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,
    "_", type,
    ".rds")))
M <- sapply(metrics, function(a) sign(a$cll$betas[-1]))
colnames(M) <- replicates
rownames(M) <- colnames(metrics[[1]]$x)
M <- M[rowSums(abs(M))>0,]
M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
type <- 'codaLasso'
p.temp <- pheatmap::pheatmap(M,
                             main=type,
                             cluster_rows = F,cluster_cols = F, silent = FALSE)


grid.newpage()
heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)
M.list[[type]] <- M

## selbal #######################################################################
if (run_selbal){
  type <- 'selbal'
  metrics = lapply(replicates, function(jid) readRDS(
    paste0(
      output_dir, file.end, "_", jid,
      "_selbal",
      ".rds")))
  All_features <- unique(Reduce(f = 'c',lapply(metrics, function(a) a$balance)))
  M <- sapply(metrics, function(m) {
    M <- rep(0, length(All_features))
    M[match(m$numerator,All_features)] <- +1
    M[match(m$denominator,All_features)] <- -1
    M
  })
  rownames(M) <- All_features
  colnames(M) <- replicates
  M <- M[rowSums(abs(M))>0,]
  
  M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
  p.temp <- pheatmap::pheatmap(M,
                               main=type,
                               cluster_rows = F,cluster_cols = F, silent = FALSE)
  
  grid.newpage()
  heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)
  M.list[[type]] <- M
  
}

## slr ----
type <- 'slr_constrainedPC'
metrics = lapply(replicates, function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid, '_', type,
    ".rds")))
All_features <- names(metrics[[1]]$`log-contrast coefficients`)
M <- sapply(metrics, function(m) {
  M <- rep(0, length(All_features))
  M[match(rownames(m$sbp),All_features)] <- m$sbp
  M
})
rownames(M) <- All_features
colnames(M) <- replicates
M <- M[rowSums(abs(M))>0,]

M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
type = 'SLR'
p.temp <- pheatmap::pheatmap(M,
                             main=type,
                             cluster_rows = F,cluster_cols = F, silent = FALSE)

grid.newpage()
heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)

M.list[[type]] <- M

getErrorRate(rownames(M.list[['SLR']])[abs(rowMeans(M.list[['SLR']]))>0.75], truth, colnames(data.list$X))


# plot ####
plt.heatmap <- plot_grid(heatmap.list[['CoDaCoRe']],
                         heatmap.list[['codaLasso']],
                         heatmap.list[['selbal']],
                         heatmap.list[['SLR']])
saveRDS(plt.heatmap, 
        file=paste0('../Outputs/figures/input4correlation/',file.end,'_heatmaps.rds'))

save_plot(plt.heatmap,
          file=paste0('../Outputs/figures/',file.end,'_heatmap.png'),base_asp = 1.25,
          base_height = 8)


plt1 <- ggplot(df, aes(x=Method, y=sensitivity, color=Method)) +
  geom_boxplot() + theme_cowplot() +
  ylab('sensitivity') + theme(axis.title.x = element_blank(),
                              axis.text.x = element_blank()) +
  geom_point(shape=16)

plt1.1 <- ggplot(df, aes(x=Method, y=specificity, color=Method)) +
  geom_boxplot() + theme_cowplot() +
  ylab('specificity') + theme(axis.title.x = element_blank(),
                              axis.text.x = element_blank()) +
  geom_point(shape=16)

plt2 <- ggplot(df, aes(x=Method, y=auc, color=Method)) +
  geom_boxplot() + theme_cowplot() +
  ylab('AUC') + theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      legend.position = 'right') +
  geom_point(shape=16)


plt3 <- ggplot(df, aes(x=Method, y=numselected, color=Method)) +
  geom_boxplot() + theme_cowplot() +
  ylab('number selected') + theme(axis.title.x = element_blank(),
                                  axis.text.x = element_blank()) +
  geom_point(shape=16)

 
plt <- plot_grid(plot_grid(plt1, plt1.1, plt2, plt3, nrow=2))

saveRDS(list(sensitivity = plt1,
             specifity = plt1.1,
             auc = plt2,
             perc = plt3),
        file=paste0('../Outputs/figures/input4correlation/',file.end,'_plots.rds'))

save_plot(plt,file=paste0('../Outputs/figures/',file.end,'_plots.png'),
          base_asp = 2,
          base_height = 6)

# for paper ----
cat('stability selection ...\n')
print(sapply(M.list, function(m) 
  getErrorRate(rownames(m)[abs(rowMeans(m))>0.75], 
               truth, colnames(data.list$X))))


