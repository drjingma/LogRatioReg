# Purpose: compare methods on data sets
# Date: 03/20/2025

rm(list=ls())

# libraries and settings

output_dir = "../Outputs/data_metrics/"

library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(grid)
library(gridExtra)

# tuning parameter settings ####
hparam = "1se"
ncomp <- 1
nreps <- 20
run_selbal <- T

dataname <- 'CRC' # IBD or CRC
cat(dataname,'\n')
if (dataname == "CRC"){
  file.end = paste0(
    "/",dataname,
    "_split", 0.7,
    "_hparam", hparam
    )
} else {
  file.end = paste0(
    "/",dataname,
    "_hparam", hparam
  )
}


metrics <- list()


# codacore - 1 balance #########################################################
type <- 'codacore'
list.metrics = lapply(seq(nreps), function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,"_", type, 
    "1_metrics",
    ".rds")))
type <- 'CoDaCoRe'
metrics[[type]] <- data.frame(Reduce(rbind, list.metrics))
rownames(metrics[[type]]) <- paste0('run',seq(nreps))

# codaLasso #######################################################################
type <- 'classo'
list.metrics = lapply(seq(nreps), function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid, "_", type, 
    "_metrics",
    ".rds")))
type <- 'codaLasso'
metrics[[type]] <- data.frame(Reduce(rbind,list.metrics))
rownames(metrics[[type]]) <- paste0('run',seq(nreps))


# selbal #######################################################################

if (run_selbal){
  type <- 'selbal'
  list.metrics = lapply(seq(nreps), function(jid) readRDS(
    paste0(
      output_dir, file.end, "_", jid, "_", type, 
      "_metrics",
      ".rds")))
  metrics[[type]] <- data.frame(Reduce(rbind,list.metrics))
  rownames(metrics[[type]]) <- paste0('run',seq(nreps))
}


# slr ----
type <- 'slr_constrainedPC'
list.metrics = lapply(seq(nreps), function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,"_", type, 
    "_metrics",
    ".rds")))
type = 'SLR'
metrics[[type]] <- data.frame(Reduce(rbind, list.metrics))
rownames(metrics[[type]]) <- paste0('run',seq(nreps))


# combine all metrics ----

test <- lapply(1:length(metrics), function(a) 
  data.frame(Method=names(metrics)[a],metrics[[a]]))
df <- Reduce(rbind,test) %>%
  dplyr::filter(Method %in% c('CoDaCoRe','codaLasso','SLR','selbal'))
head(df)

# Variable selection heat maps ----
M.list <- list()
heatmap.list <- list()

## codacore - 1 balance #########################################################
type <- "codacore" 
metrics <- lapply(seq(nreps), function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,
    "_", type, 
    "1.rds")))
M <- Reduce(cbind,lapply(seq(nreps), function(jid)
  metrics[[jid]]$ensemble[[1]]$hard$numerator - metrics[[jid]]$ensemble[[1]]$hard$denominator
))
colnames(M) <- seq(nreps)
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

## classo ----
type <- 'classo'
metrics = lapply(seq(nreps), function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid,
    "_", type,
    ".rds")))
M <- matrix(0, nrow=ncol(metrics[[1]]$x), ncol=nreps)
rownames(M) <- colnames(metrics[[1]]$x)
colnames(M) <- seq(nreps)
for (a in seq(nreps)){
  M[,a] <- sign(metrics[[a]]$cll$betas[-1])
}
M <- M[rowSums(abs(M))>0,]
M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
type <- 'codaLasso'
p.temp <- pheatmap::pheatmap(M,
                             main=type,
                             cluster_rows = F,cluster_cols = F, silent = FALSE)

grid.newpage()
heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)
M.list[[type]] <- M

if (run_selbal){
  ## selbal #######################################################################
  type <- 'selbal'
  metrics = lapply(seq(nreps), function(jid) readRDS(
    paste0(
      output_dir, file.end, "_", jid,
      "_selbal",
      ".rds")))
  All_features <- unique(Reduce(f = 'c',lapply(metrics, function(a) a$global.balance$Taxa)))
  M <- matrix(0, nrow=length(All_features), ncol=nreps)
  rownames(M) <- All_features
  colnames(M) <- seq(nreps)
  for (a in seq(nreps)){
    M[match(metrics[[a]]$global.balance$Taxa,All_features),a] <- ifelse(metrics[[a]]$global.balance$Group == 'NUM', +1, -1)
  }
  M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
  p.temp <- pheatmap::pheatmap(M,
                               main=type,
                               cluster_rows = F,cluster_cols = F, silent = FALSE)
  
  grid.newpage()
  heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)
  M.list[[type]] <- M
  
}

## slr - PC ----
type <- 'slr_constrainedPC'
metrics = lapply(seq(nreps), function(jid) readRDS(
  paste0(
    output_dir, file.end, "_", jid, '_', type,
    ".rds")))
All_features <- unique(Reduce(f = 'c',lapply(metrics, function(a) rownames(a$sbp))))
M <- matrix(0, nrow=length(All_features), ncol=nreps)
rownames(M) <- All_features
colnames(M) <- seq(nreps)
for (a in seq(nreps)){
  M[match(rownames(metrics[[a]]$sbp),All_features),a] <- metrics[[a]]$sbp
}
M <- M[rowSums(abs(M))>0,]
M <- M[order(rowMeans(M==-1),rowMeans(M==0)),]
type = 'SLR'
p.temp <- pheatmap::pheatmap(M,
                             main=type,
                             cluster_rows = F,cluster_cols = F, silent = FALSE)
grid.newpage()
heatmap.list[[type]] <- arrangeGrob(p.temp$gtable)

M.list[[type]] <- M

saveRDS(M.list, file=paste0('../Outputs/figures/',file.end,'_sbp.rds'))

# summary ----
heatmap.grid <- plot_grid(heatmap.list[['CoDaCoRe']],
                          heatmap.list[['codaLasso']],
                          heatmap.list[['selbal']],
                          heatmap.list[['SLR']], nrow=2)
save_plot(heatmap.grid,
          file=paste0('../Outputs/figures/',file.end,'_heatmap.png'),base_asp = 1.25,
          base_height = 8)

# plot ####
df$numselected <- c(sapply(M.list,function(a) colSums(abs(a))))

plt2 <- ggplot(df, aes(x=Method, y=auc, color=Method)) + 
  geom_boxplot() + theme_cowplot() +
  ylab('AUC') + theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank()) +
  geom_point(shape=16)


plt3 <- ggplot(df, aes(x=Method, y=numselected, color=Method)) + 
  geom_boxplot() + theme_cowplot() +
  ylab('number selected') + theme(axis.title.x = element_blank(),
                                  axis.text.x = element_blank()) +
  geom_point(shape=16)

legend = cowplot::get_plot_component(plt2, 'guide-box-right', return_all = TRUE)

plt2 <- plt2 + theme(legend.position="none")
plt3 <- plt3 + theme(legend.position="none")

plt <- plot_grid(heatmap.grid,
                 plot_grid(plot_grid(plt2, plt3, nrow=2,
                                     labels = c('B','C'), scale = 0.9),
                           legend,
                           nrow = 1,
                           rel_widths = c(0.65,.35)),
                 labels = c('A',''), scale = 0.95,
                 rel_widths = c(0.6,0.4)
)
plt
save_plot(plt,file=paste0('../Outputs/figures/',file.end,'.png'),
          base_asp = 2,
          base_height = 5)

## For Venndiagram

stable.variables <- lapply(M.list, function(m)
  names(which(rowMeans(abs(m)) >= 0.75))
)

saveRDS(stable.variables,  paste0('../Figures/', file.end, '_variables.rds'))

