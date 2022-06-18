# Purpose: plot results from HIV_mse.R
# Date: 6/17/2022
rm(list=ls())

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs_mse"

source("slr_analyses/Functions/util.R")

library(tidyverse)
library(reshape2)

numSplits = 20

# tuning parameter settings
K = 10
scaling = TRUE

file.end0 = paste0(
  # "_sim", b,
  "_HIV", 
  "_gbm")

################################################################################
# plot metrics

# import metrics
classo_list = list()
slr_list = list()
selbal_list = list()
codacore_list = list()
lrlasso_list = list()
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
  
  # slr
  slrscreen_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(slrscreen_tmp) = NULL
  slr_list[[i]] = data.table::data.table(slrscreen_tmp)
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
  
  # log-ratio lasso
  lrl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/lrlasso_metrics",
    "_sim", i, file.end0, ".rds"
  ))))
  rownames(lrl_tmp) = NULL
  lrlasso_list[[i]] = data.table::data.table(lrl_tmp)
}

# metrics boxplots
classo.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(classo_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo")
###
slrscreen.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr_list)), 
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
lrlasso.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(lrlasso_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "lrlasso")
###
data.gg = rbind(
  classo.gg,
  slrscreen.gg,
  selbal.gg, 
  codacore.gg, 
  lrlasso.gg
) %>%
  mutate(
    Metric = factor(
      Metric, 
      levels = c(
        "acc", "auc", "f1", "percselected",  "time"
      ), 
      labels = c(
        "Accuracy", "AUC", "F1", "% Selected", "Timing"
      ))
  )

data.gg_main = data.gg
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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank())
plt_main
ggsave(
  filename = paste0(
    "20220617",
    file.end0,
    "_", "metrics", ".pdf"),
  plot = plt_main,
  width = 6, height = 3.5, units = c("in")
)
