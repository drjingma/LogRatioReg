# This data set is from Franzosa et al. (19') on IBD and includes data from two cohorts.  
# The number of genera is very large. 
# To reduce the dimensionality, we filter taxa at mean abundance 0.005\% and prevalence of 0.5.
library(phyloseq)
library(tidyverse)

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
  
  genera <- phyloseq(otu_table(OTU,taxa_are_rows = FALSE),tax_table(TAX),sample_data(meta))
  genera
}

# load data
data_dir <- "../Data/Franzosa_IBD_2019/"
genera.counts <- read_delim(paste0(data_dir, "genera.counts.tsv"),delim = '\t')
meta <- read_delim(paste0(data_dir, "metadata.tsv"),delim = "\t",)
meta <- as.data.frame(meta)
meta$Cohort <- sapply(meta$Sample,function(s) unlist(strsplit(s,".",fixed = TRUE))[1])

identical(genera.counts$Sample,meta$Sample) # TRUE

# create a phyloseq table
genera <- create_phyloseq(genera.counts,meta)
genera

# feature filtering 
genera1 = transform_sample_counts(genera, function(x) x / sum(x) )
genera2 <- phyloseq::filter_taxa(genera1,  function(x){ (sum(x > 0) > nsamples(genera1)*0.5) && (mean(x) > 5*1e-5)  }, prune = T)
genera3 <- phyloseq::prune_taxa(rownames(tax_table(genera)) %in% rownames(tax_table(genera2)),genera)
genera3

saveRDS(tax_table(genera), file='../Data/Franzosa_IBD_2019/taxonomy.rds')

# obtain otu count table
count_table <- otu_table(genera3)
count_table <- as(count_table,'matrix')
dim(count_table)
summary(rowSums(count_table))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 126087  16777852  23911738  31569255  28212911 109063353 
summary(colMeans(count_table==0))

# save data for regression 
Franzosa_PRISM <- list()
Franzosa_Validation <- list()
PRISM.idx <- (meta$Cohort=='PRISM')
Validation.idx <- (meta$Cohort=='Validation')

CD.idx <- (meta$Study.Group %in% c('Control','CD'))

Franzosa_PRISM$X <- count_table[PRISM.idx*CD.idx==1,]
mean(Franzosa_PRISM$X==0) # 4\%
summary(colMeans(Franzosa_PRISM$X==0))
Franzosa_PRISM$y <- meta$Study.Group[PRISM.idx*CD.idx==1]
save(Franzosa_PRISM,file=paste0(data_dir,"Franzosa_PRISM_CD.rda"))

Franzosa_Validation$X <- count_table[Validation.idx*CD.idx==1,]
mean(Franzosa_Validation$X==0) # 4\%
summary(colMeans(Franzosa_Validation$X==0))
Franzosa_Validation$y <- meta$Study.Group[Validation.idx*CD.idx==1]
save(Franzosa_Validation,file=paste0(data_dir, "Franzosa_Validation_CD.rda"))


UC.idx <- (meta$Study.Group %in% c('Control','UC'))
Franzosa_PRISM$X <- count_table[PRISM.idx*UC.idx==1,]
mean(Franzosa_PRISM$X==0) # 3\%
summary(colMeans(Franzosa_PRISM$X==0))
Franzosa_PRISM$y <- meta$Study.Group[PRISM.idx*UC.idx==1]
save(Franzosa_PRISM,file=paste0(data_dir, "Franzosa_PRISM_UC.rda"))
Franzosa_Validation$X <- count_table[Validation.idx*UC.idx==1,]
mean(Franzosa_Validation$X==0) # 3\%
summary(colMeans(Franzosa_Validation$X==0))
Franzosa_Validation$y <- meta$Study.Group[Validation.idx*UC.idx==1]

save(Franzosa_Validation,file=paste0(data_dir,"Franzosa_Validation_UC.rda"))

