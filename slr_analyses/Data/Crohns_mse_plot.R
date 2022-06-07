rm(list=ls())
# Purpose: compare methods on prediction accuracy & timing using Crohns data
# Date: 5/26/2022

################################################################################
# libraries and settings

output_dir = "Kristyn/Data/outputs_mse"

source("Kristyn/Functions/util.R")

library(tidyverse)
library(reshape2)

numSplits = 20

# tuning parameter settings
K = 10
scaling = TRUE

file.end0 = paste0(
  # "_sim", b,
  "_Crohns", 
  "_gbm")

################################################################################
# plot metrics

# import metrics
classo_list = list()
slr_0.05_list = list()
slr_0.01_list = list()
slrscreen_list = list()
selbal_list = list()
codacore_list = list()
for(i in 1:numSplits){
  print(i)
  
  # compositional lasso
  cl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(cl_tmp) = NULL
  classo_list[[i]] = data.table::data.table(cl_tmp)
  
  ###
  
  # slr - alpha = 0.05
  slr_alpha0.05_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_alpha0.05_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(slr_alpha0.05_tmp) = NULL
  slr_0.05_list[[i]] = data.table::data.table(slr_alpha0.05_tmp)
  
  # slr - alpha = 0.01
  slr_alpha0.01_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_alpha0.01_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(slr_alpha0.01_tmp) = NULL
  slr_0.01_list[[i]] = data.table::data.table(slr_alpha0.01_tmp)
  
  # slr - screen
  slrscreen_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slrscreen_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(slrscreen_tmp) = NULL
  slrscreen_list[[i]] = data.table::data.table(slrscreen_tmp)
  ###
  
  # selbal
  slbl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/selbal_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(slbl_tmp) = NULL
  selbal_list[[i]] = data.table::data.table(slbl_tmp)
  
  # codacore
  cdcr_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/codacore_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(cdcr_tmp) = NULL
  codacore_list[[i]] = data.table::data.table(cdcr_tmp)
}

# metrics boxplots
classo.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(classo_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo")
###
slr_0.05.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(slr_0.05_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-0.05")
slr_0.01.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(slr_0.01_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-0.01")
slrscreen.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slrscreen_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr")
###
selbal.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(selbal_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "selbal")
codacore.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(codacore_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "codacore")
###
data.gg = rbind(
  classo.gg,
  # slr_0.05.gg,
  # slr_0.01.gg,
  slrscreen.gg,
  selbal.gg, 
  codacore.gg
) %>%
  mutate(
    Metric = factor(
      Metric, 
      levels = c(
        "acc", "auc", "time"
      ), 
      labels = c(
        "Accuracy", "AUC", "Timing"
      ))
  )

data.gg_main = data.gg %>%
  dplyr::filter(
    Metric %in% c(
      "Accuracy", "AUC", "Timing"
    )
  )
plt_main = ggplot(
  data.gg_main, 
  aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(Metric), scales = "free_y") +
  geom_boxplot() +
  stat_summary(
    fun = mean, geom = "point", shape = 4, size = 1.5,
    color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank())
plt_main
ggsave(
  filename = paste0(
    "20220530",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 6, height = 2.5, units = c("in")
)
