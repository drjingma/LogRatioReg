# This data set is from Franzosa et al. (19') on IBD and includes data from two cohorts.  
# The number of genera is very large. 
# To reduce the dimensionality, we filter taxa at mean abundance 0.005\% and prevalence of 0.5.
library(phyloseq)
library(tidyverse)

rm(list=ls())
#' @param f a data frame 
#' @param meta metadata
create_phyloseq <- function(f,meta){
  genera <- colnames(f)[-1]
  TAX <- Reduce(rbind,lapply(genera, function(s) unlist(strsplit(s,';'))))
  colnames(TAX) <- c("Kingdom","Phylum","Class","Order","Family","Genus")
  rownames(TAX) <- NULL
  rownames(TAX) <- paste0("g",1:nrow(TAX))
  OTU <- apply(f[,-1], 2, as.numeric)
  colnames(OTU) <- paste0("g",1:nrow(TAX))
  rownames(OTU) <- paste0("sa",1:nrow(OTU))
  
  rownames(meta) <- paste0("sa",1:nrow(OTU))
  genera <- phyloseq(otu_table(OTU,taxa_are_rows = FALSE),tax_table(TAX),sample_data(meta))
  genera
}

# load data
filename <- 'Yachida_CRC_2019'
data_dir <- paste0("../Data/",filename,'/')
genera.counts <- read_delim(paste0(data_dir, "genera.counts.tsv"),delim = '\t')
meta <- read_delim(paste0(data_dir, "metadata.tsv"),delim = "\t",)
meta <- as.data.frame(meta)

identical(genera.counts$Sample,meta$Sample) # FALSE
meta <- meta[match(genera.counts$Sample,meta$Sample),]

# create a phyloseq table
genera <- create_phyloseq(genera.counts,meta)
genera
count_table <- otu_table(genera)
count_table <- as(count_table,'matrix')
summary(rowSums(count_table))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 2892954 28574264 33823006 34744654 40426291 64470852 

saveRDS(tax_table(genera), file='../Data/Yachida_CRC_2019/CRC_taxonomy.rds')

# feature filtering 
genera1 = transform_sample_counts(genera, function(x) x / sum(x) )
genera2 <- phyloseq::filter_taxa(genera1,  function(x){ (sum(x > 0) > nsamples(genera1)*0.3) && (mean(x) > 5*1e-5)}, prune = T)
genera3 <- phyloseq::prune_taxa(rownames(tax_table(genera)) %in% rownames(tax_table(genera2)),genera)
genera3

# obtain otu count table
count_table <- otu_table(genera3)
count_table <- as(count_table,'matrix')
dim(count_table)
mean(count_table==0)
summary(colMeans(count_table==0))

# save data for regression 
data.list <- list()
idx <- meta$Study.Group %in% c('Healthy', 'Stage_I_II', 'Stage_III_IV')
data.list$X <- count_table[idx,]
summary(colMeans(data.list$X==0))
mean(data.list$X==0)
data.list$y <- ifelse(meta$Study.Group[idx] %in% 'Healthy', 'Control', 'CRC')

save(data.list,file=paste0(data_dir,filename,".rda"))


