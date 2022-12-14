# Purpose: plot results from HIV_mse.R
# Date: 8/10/2022
rm(list=ls())

data_set = "HIV" # "HIV", "sCD14", "Crohn"
date = "20221214"

################################################################################
# libraries and settings

output_dir = "Data/outputs_metrics"

source("Functions/util.R")

library(tidyverse)
library(reshape2)

numSplits = 20

# tuning parameter settings
hparam = "1se"
filter.perc = 0.8 # 0.8, 1
split.perc = 0.7  # 0.7, 0.8

# tuning parameter settings
K = 10
scaling = TRUE

file.end0 = paste0(
  "_", data_set,
  "_split", split.perc, 
  "_filter", filter.perc,
  "_hparam", hparam,
  "_gbm")

################################################################################
# plot metrics
if(data_set == "Crohn"){
  W = readRDS(paste0(
    output_dir, "/data",
    file.end0, "_sim", 1, ".rds"
  ))$XTr
} else if(data_set == "HIV"){
  W = readRDS(paste0(
    output_dir, "/data",
    file.end0, "_sim", 1, ".rds"
  ))$XTr
} else if(data_set == "sCD14"){
  W = readRDS(paste0(
    output_dir, "/data",
    file.end0, "_sim", 1, ".rds"
  ))$XTr
}
p = ncol(W)

# import metrics
slr_spec_sbps = matrix(NA, nrow = p, ncol = numSplits)
slr_hier_sbps = matrix(NA, nrow = p, ncol = numSplits)
selbal_sbps = matrix(NA, nrow = p, ncol = numSplits)
codacore1_sbps = matrix(NA, nrow = p, ncol = numSplits) # codacore has up to 5 balances selected...
# codacore2_sbps = matrix(NA, nrow = p, ncol = numSplits)
rownames(slr_spec_sbps) <- rownames(slr_hier_sbps) <- rownames(selbal_sbps) <- 
  rownames(codacore1_sbps) <- 
  # rownames(codacore2_sbps) <- 
  colnames(W)
for(i in 1:numSplits){
  print(i)
  
  # slr - spectral clustering
  slr_spec_sbps[, i] = readRDS(paste0(
    output_dir, "/slr_spectral_sbp",
    file.end0, "_sim", i, ".rds"
  ))
  
  # slr - hierarchical clustering
  slr_hier_sbps[, i] = readRDS(paste0(
    output_dir, "/slr_hierarchical_sbp",
    file.end0, "_sim", i, ".rds"
  ))
  
  ###
  
  # selbal
  selbal_sbps[, i] = readRDS(paste0(
    output_dir, "/selbal_sbp",
    file.end0, "_sim", i, ".rds"
  ))
  
  # codacore
  codacore_sbp_tmp = readRDS(paste0(
    output_dir, "/codacore1_sbp",
    file.end0, "_sim", i, ".rds"
  ))
  print(paste0("ncol(codacore_sbp) = ", ncol(codacore_sbp_tmp)))
  
  codacore1_sbps[, i] = codacore_sbp_tmp[, 1]
  # if(ncol(codacore_sbp_tmp) >= 2){
  #   codacore2_sbps[, i] = codacore_sbp_tmp[, 2]
  # }
  
}
################################################################################
# get active sets and selected balances (if applicable)
################################################################################

# slr - spectral ###############################################################
slr_spec_props0 = data.frame(
  taxa = rownames(slr_spec_sbps),
  active = apply(slr_spec_sbps != 0, 1, function(row) mean(row)), 
  numerator = apply(slr_spec_sbps == 1, 1, function(row) mean(row)), 
  denominator = apply(slr_spec_sbps == -1, 1, function(row) mean(row))
)
slr_spec_props = slr_spec_props0 %>% filter(active != 0) %>%
  arrange(active) %>% 
  dplyr::select(!(active)) %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
slr_spec_props = slr_spec_props %>% 
  mutate(taxa = factor(taxa, levels = unique(slr_spec_props$taxa)))

# slr - hierarchical ###########################################################
slr_hier_props0 = data.frame(
  taxa = rownames(slr_hier_sbps),
  active = apply(slr_hier_sbps != 0, 1, function(row) mean(row)), 
  numerator = apply(slr_hier_sbps == 1, 1, function(row) mean(row)), 
  denominator = apply(slr_hier_sbps == -1, 1, function(row) mean(row))
)
slr_hier_props = slr_hier_props0 %>% filter(active != 0) %>%
  arrange(active) %>% 
  dplyr::select(!(active)) %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
slr_hier_props = slr_hier_props %>% 
  mutate(taxa = factor(taxa, levels = unique(slr_hier_props$taxa)))

# selbal #######################################################################
selbal_props0 = data.frame(
  taxa = rownames(selbal_sbps),
  active = apply(selbal_sbps != 0, 1, function(row) mean(row)), 
  numerator = apply(selbal_sbps == 1, 1, function(row) mean(row)), 
  denominator = apply(selbal_sbps == -1, 1, function(row) mean(row))
)
selbal_props = selbal_props0 %>% filter(active != 0) %>%
  arrange(active) %>% 
  dplyr::select(!(active)) %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
selbal_props = selbal_props %>% 
  mutate(taxa = factor(taxa, levels = unique(selbal_props$taxa)))

# codacore1 ####################################################################
codacore1_props0 = data.frame(
  taxa = rownames(codacore1_sbps),
  active = apply(codacore1_sbps != 0, 1, function(row) mean(row)), 
  numerator = apply(codacore1_sbps == 1, 1, function(row) mean(row)), 
  denominator = apply(codacore1_sbps == -1, 1, function(row) mean(row))
)
codacore1_props = codacore1_props0 %>% filter(active != 0) %>%
  arrange(active) %>% 
  dplyr::select(!(active)) %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
codacore1_props = codacore1_props %>% 
  mutate(taxa = factor(taxa, levels = unique(codacore1_props$taxa)))

# combine ######################################################################

method_props = rbind(
  cbind(slr_spec_props, method = "slr-spec"), 
  # cbind(slr_spec_props, method = "slr-hier"), 
  cbind(selbal_props, method = "selbal"), 
  cbind(codacore1_props, method = "codacore") 
) %>% 
  arrange(side, desc(proportion)) %>% 
  mutate(method = factor(
    method, labels = c("slr-spec", "selbal", "codacore"), 
    levels = rev(c("selbal", "codacore", "slr-spec"))))

barplt = ggplot(
  method_props, aes(x = reorder(taxa, -proportion, function(vec) sum(vec, na.rm = TRUE)), y = proportion, fill = side)) + 
  facet_wrap(vars(method), ncol = 1) +
  geom_bar(stat = "identity") + 
  # coord_flip() + 
  theme_bw() +
  theme(
    text = element_text(size=10),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, 1))
barplt
ggsave(
  filename = paste0(
    date,
    file.end0,
    "_", "selectionbars", ".png"),
  plot = barplt,
  width = 6, height = 9, units = c("in")
)
