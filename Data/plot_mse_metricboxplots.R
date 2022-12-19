# Purpose: plot results from HIV_mse.R
# Date: 8/25/2022
rm(list=ls())

data_set = "sCD14" # "HIV", "sCD14", "Crohn"
date = "20221218"

response_type = NA
if(data_set %in% c("sCD14")){
  response_type = "continuous"
} else{
  response_type = "binary"
}

logtime = TRUE
label_means = TRUE

remove_xaxis = TRUE

################################################################################
# libraries and settings

output_dir = "Data/outputs_metrics"

source("Functions/util.R")

library(tidyverse)
library(reshape2)
library(ggrepel)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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
slr_spec_list = list()
slr_hier_list = list()
selbal_list = list()
# selbal_covar_list = list() # only applicable to HIV data set
codacore_list = list()
lrlasso_list = list()
for(i in (1:numSplits)){
  print(i)
  
  # compositional lasso (MSE for continuous response, AUC for binary response)
  cl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/classo_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(cl_tmp) = NULL
  classo_list[[i]] = data.table::data.table(cl_tmp)
  
  ###
  
  # slr - spectral clustering
  slr_spec_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_spectral_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(slr_spec_tmp) = NULL
  slr_spec_list[[i]] = data.table::data.table(slr_spec_tmp)
  
  # slr - hierarchical clustering
  slr_hier_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/slr_hierarchical_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(slr_hier_tmp) = NULL
  slr_hier_list[[i]] = data.table::data.table(slr_hier_tmp)
  
  ###
  
  # selbal
  slbl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/selbal_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(slbl_tmp) = NULL
  selbal_list[[i]] = data.table::data.table(slbl_tmp)
  
  # # selbal - covar
  # slbl_covar_tmp = t(data.frame(readRDS(paste0(
  #   output_dir, "/selbal_covar_metrics",
  #   file.end0, "_sim", i, ".rds"
  # ))))
  # rownames(slbl_covar_tmp) = NULL
  # selbal_covar_list[[i]] = data.table::data.table(slbl_covar_tmp)
  
  #
  
  # codacore
  cdcr_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/codacore1_metrics",
    file.end0, "_sim", i, ".rds"
  ))))
  rownames(cdcr_tmp) = NULL
  codacore_list[[i]] = data.table::data.table(cdcr_tmp)
  
  # log-ratio lasso
  lrl_tmp = t(data.frame(readRDS(paste0(
    output_dir, "/lrlasso_metrics",
    file.end0, "_sim", i, ".rds"
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
slr_spec.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(slr_spec_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-spec")
slr_hier.gg = 
  pivot_longer(as.data.frame(data.table::rbindlist(slr_hier_list)), 
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "slr-hier")
###
selbal.gg =
  pivot_longer(as.data.frame(data.table::rbindlist(selbal_list)),
               cols = everything(),
               names_to = "Metric") %>%
  mutate("Method" = "selbal")

# selbal_covar.gg = 
#   pivot_longer(as.data.frame(data.table::rbindlist(selbal_covar_list)), 
#                cols = everything(),
#                names_to = "Metric") %>%
#   mutate("Method" = "selbal-MSM")
###
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
  slr_spec.gg,
  # slr_hier.gg,
  selbal.gg,
  # selbal_covar.gg,
  codacore.gg, 
  lrlasso.gg
)

if(response_type == "binary"){
  prediction_metric = "auc"
  prediction_metric_label = "AUC"
} else if(response_type == "continuous"){
  prediction_metric = "mse"
  prediction_metric_label = "MSE"
}

data.gg = data.gg %>% 
  dplyr::filter(
    Metric %in% c(prediction_metric, "percselected", "time")
  ) %>%
  mutate(
    Metric = factor(
      Metric, 
      levels = c(
        prediction_metric, "percselected",  "time"
      ), 
      labels = c(
        prediction_metric_label, "% Selected", "Time (s)"
      ))
  ) 

data.gg = data.gg %>% 
  mutate(
    Method = factor(
      Method, 
      levels = c(
        "selbal", "classo", "codacore", "lrlasso", 
        "slr-spec", "slr-hier"
      )
    )
  ) 

if(logtime){
  data.gg = data.gg %>% 
    mutate(
      value = ifelse(Metric == "Time (s)", log(value), value)
    ) %>%
    mutate(
      Metric = factor(
        Metric, 
        levels = c(
          prediction_metric_label, "% Selected", "Time (s)"
        ), 
        labels = c(
          prediction_metric_label, "% Selected", "log(Time (s))"
        ))
    ) 
}

data.gg_main = data.gg 
means.gg = data.gg_main %>% 
  group_by(Metric) %>%
  mutate(
    yrange = abs(max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
  ) %>%
  group_by(Metric, Method) %>% 
  dplyr::summarize(mean = signif(mean(value, na.rm = TRUE), 2), yrange = first(yrange))

plt_colors = gg_color_hue(6)[1:5]
names(plt_colors) = c("selbal", "classo", "codacore", "lrlasso", "slr-spec")
plt_main = ggplot(
  data.gg_main, 
  aes(x = Method, y = value, color = Method)) +
  facet_wrap(vars(Metric), scales = "free_y") +
  geom_boxplot() +
  stat_summary(
    fun = mean, geom = "point", shape = 4, size = 1.5,
    color = "red")
if(label_means){
  plt_main = plt_main  +
    geom_text_repel(
      data = means.gg, aes(label = mean, y = mean), # + 0.05 * yrange), 
      size = 2, color = "black")
  filename = paste0(
    date,
    file.end0,
    "_", "metrics_labeledmeans")
} else{
  filename = paste0(
    date,
    file.end0,
    "_", "metrics")
}
if(logtime){
  filename = paste0(
    filename, 
    "_logtime", 
    ".png"
  )
} else{
  filename = paste0(
    filename, 
    ".png"
  )
}
plt_main = plt_main +
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) + 
  scale_color_manual(values = plt_colors)
plt_main

ggsave(
  filename = filename,
  plot = plt_main,
  width = 6, height = 2.75, units = c("in") # height is 2.5 usually
)
