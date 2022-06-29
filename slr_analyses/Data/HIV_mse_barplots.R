# Purpose: plot results from HIV_mse.R
# Date: 6/28/2022
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
W = selbal::HIV[, 1:60]
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
    "_sim", i, file.end0, ".rds"
  ))
  
  # slr - hierarchical clustering
  slr_hier_sbps[, i] = readRDS(paste0(
    output_dir, "/slr_hierarchical_sbp",
    "_sim", i, file.end0, ".rds"
  ))
  
  ###
  
  # selbal
  selbal_sbps[, i] = readRDS(paste0(
    output_dir, "/selbal_sbp",
    "_sim", i, file.end0, ".rds"
  ))
  
  # codacore
  codacore_sbp_tmp = readRDS(paste0(
    output_dir, "/codacore_sbp",
    "_sim", i, file.end0, ".rds"
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
  select(-"active") %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
slr_spec_props = slr_spec_props %>% 
  mutate(taxa = factor(taxa, levels = slr_spec_props$taxa))

slr_spec_bar = ggplot(
  slr_spec_props, aes(x = taxa, y = proportion, fill = side)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) +
  ggtitle("slr-spec")
ggsave(
  filename = paste0(
    "20220622",
    file.end0,
    "_", "slr_spectral", ".pdf"),
  plot = slr_spec_bar,
  width = 6, height = 3.5, units = c("in")
)

# slr - hierarchical ###########################################################
slr_hier_props0 = data.frame(
  taxa = rownames(slr_hier_sbps),
  active = apply(slr_hier_sbps != 0, 1, function(row) mean(row)), 
  numerator = apply(slr_hier_sbps == 1, 1, function(row) mean(row)), 
  denominator = apply(slr_hier_sbps == -1, 1, function(row) mean(row))
)
slr_hier_props = slr_hier_props0 %>% filter(active != 0) %>%
  arrange(active) %>% 
  select(-"active") %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
slr_hier_props = slr_hier_props %>% 
  mutate(taxa = factor(taxa, levels = slr_hier_props$taxa))

slr_hier_bar = ggplot(
  slr_hier_props, aes(x = taxa, y = proportion, fill = side)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) +
  ggtitle("slr-hier")
ggsave(
  filename = paste0(
    "20220622",
    file.end0,
    "_", "slr_hierarchical", ".pdf"),
  plot = slr_hier_bar,
  width = 6, height = 3.5, units = c("in")
)

# selbal #######################################################################
selbal_props0 = data.frame(
  taxa = rownames(selbal_sbps),
  active = apply(selbal_sbps != 0, 1, function(row) mean(row)), 
  numerator = apply(selbal_sbps == 1, 1, function(row) mean(row)), 
  denominator = apply(selbal_sbps == -1, 1, function(row) mean(row))
)
selbal_props = selbal_props0 %>% filter(active != 0) %>%
  arrange(active) %>% 
  select(-"active") %>% 
  pivot_longer(
    cols = numerator:denominator, 
    names_to = "side", values_to = "proportion") %>%
  filter(proportion != 0)
selbal_props = selbal_props %>% 
  mutate(taxa = factor(taxa, levels = selbal_props$taxa))

selbal_bar = ggplot(
  selbal_props, aes(x = taxa, y = proportion, fill = side)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_bw() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_blank()) +
  ggtitle("selbal")
ggsave(
  filename = paste0(
    "20220622",
    file.end0,
    "_", "selbal", ".pdf"),
  plot = selbal_bar,
  width = 6, height = 3.5, units = c("in")
)





# codacore #####################################################################










