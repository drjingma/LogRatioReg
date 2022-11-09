# Purpose: plot results from HIV_mse.R
# Date: 8/25/2022
rm(list=ls())

data_set = "Crohn" # "HIV", "Crohn"
date = "20221108"

response_type = "binary"

label_means = TRUE

################################################################################
# libraries and settings

output_dir = "slr_analyses/Data/outputs_metrics"

source("slr_analyses/Functions/util.R")

library(tidyverse)
library(reshape2)
library(ggrepel)

numSplits = 20

# tuning parameter settings
hparam = "1se"
filter.perc = 0.8 # 0.8, 1
split.perc = 0.7  # 0.7, 0.8

file.end0 = paste0(
  "_", data_set,
  "_split", split.perc, 
  "_filter", filter.perc,
  "_hparam", hparam,
  "_gbm"
)

################################################################################
# plot metrics

# import metrics
classo_list = list()
classo0_list = list()
for(i in (1:numSplits)){
  print(i)
  
  # compositional lasso (MSE for continuous response, AUC for binary response)
  cl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo1_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(cl_tmp) = NULL
  classo_list[[i]] = data.table::data.table(cl_tmp)
  
  # compositional lasso with Gordon-Rodriguez' validation metric
  cl0_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo0_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(cl0_tmp) = NULL
  classo0_list[[i]] = data.table::data.table(cl0_tmp)
  
}

# metrics boxplots
classo.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(classo_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo-auc")
classo0.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(classo0_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "classo-orig")
###
data.gg = rbind(
  classo.gg,
  classo0.gg
)
data.gg = data.gg %>% 
  dplyr::filter(
    Metric %in% c("auc", "percselected", "time")
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

data.gg = data.gg %>% 
  mutate(
    Method = factor(
      Method, 
      levels = c(
        "classo-auc", "classo-orig"
      )
    )
  ) 

data.gg_main = data.gg 
means.gg = data.gg_main %>% 
  group_by(Metric) %>%
  mutate(
    yrange = abs(max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) %>%
  group_by(Metric, Method) %>% 
  dplyr::summarize(mean = signif(mean(value, na.rm = TRUE), 2), yrange = first(yrange))
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

if(label_means){
  plt_main = plt_main  +
    geom_text_repel(
      data = means.gg, aes(label = mean, y = mean), # + 0.05 * yrange), 
      size = 2, color = "black")
  filename = paste0(
    date,
    file.end0,
    "_", "metrics_codalasso_labeledmeans", ".png")
} else{
  filename = paste0(
    date,
    file.end0,
    "_", "metrics_codalasso", ".png")
}


plt_main
ggsave(
  filename = filename,
  plot = plt_main,
  width = 6, height = 2.75, units = c("in") # height is 2.5 usually
)
