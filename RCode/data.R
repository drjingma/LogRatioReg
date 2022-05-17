library(curatedMetagenomicData)
LC <- curatedMetagenomicData("LeChatelierE_2013.metaphlan_bugs_list.stool",dryrun = F)
head(pData(LC[[1]]))
# First run PCoA to visualize the association between data and BMI.


library(tidyverse)
abundance <- read.csv("Data/abundance/abundance_Obesity.txt",header = TRUE,sep = '\t',row.names = 1,skip = 1)
marker <- read.csv("Data/marker/marker_Obesity.txt",sep='\t')


rm(list=ls())
library(phyloseq)
library(ggplot2)
library(glmnet)
library(gdata)
library(igraph)
library(evd)
library(RColorBrewer)
library(fields)
library(gridExtra)
library(cowplot)
theme_set(theme_cowplot())
library(vcd)
library(pheatmap)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

fileFolder <- "../../TestBMN/0_Data/"

load(paste0(fileFolder,"twins_rna_seq_genus_ready.rda"))
sample_info <- read.csv(paste0(fileFolder,"CHM_2016_upload_for_EBI.txt"), header = T, sep="\t")
sample_info$sample_alias <- as.character(sample_info$sample_alias)
sample_info$family_id <- as.character(sample_info$family_id)

dataset$sample_info <- sample_info[match(dataset$sample_info$sample_alias, sample_info$sample_alias),]

## Compile taxa info
dataset$taxa_info$Genus <- sapply(dataset$taxa_info$Genus, function(s) unlist(strsplit(s, "g__"))[2])
dataset$taxa_info$Kingdom <- sapply(dataset$taxa_info$Details, function(s) unlist(strsplit(s, "; "))[1])
dataset$taxa_info$Phylum <- sapply(dataset$taxa_info$Details, function(s) unlist(strsplit(s, "; "))[2])
dataset$taxa_info$Class <- sapply(dataset$taxa_info$Details, function(s) unlist(strsplit(s, "; "))[3])
dataset$taxa_info$Order <- sapply(dataset$taxa_info$Details, function(s) unlist(strsplit(s, "; "))[4])
dataset$taxa_info$Family <- sapply(dataset$taxa_info$Details, function(s) unlist(strsplit(s, "; "))[5])
dataset$taxa_info$Kingdom <- sapply(dataset$taxa_info$Kingdom, function(s) unlist(strsplit(s, "k__"))[2])
dataset$taxa_info$Phylum <- sapply(dataset$taxa_info$Phylum, function(s) unlist(strsplit(s, "p__"))[2])
dataset$taxa_info$Class <- sapply(dataset$taxa_info$Class, function(s) unlist(strsplit(s, "c__"))[2])
dataset$taxa_info$Order <- sapply(dataset$taxa_info$Order, function(s) unlist(strsplit(s, "o__"))[2])
dataset$taxa_info$Family <- sapply(dataset$taxa_info$Family, function(s) unlist(strsplit(s, "f__"))[2])
dataset$taxa_info$Details <- NULL
dataset$taxa_info <- dataset$taxa_info[,c(2,3,4,5,6,1)]
colnames(dataset$dat) <- dataset$taxa_info$Genus

# Build a phyloseq object
OTU <- otu_table(t(dataset$dat), taxa_are_rows = T)
taxmat <- as.matrix(dataset$taxa_info)
rownames(taxmat) <- dataset$taxa_info$Genus
taxa <- tax_table(taxmat)

rownames(dataset$sample_info) <- rownames(dataset$dat)
sampledata <- sample_data(dataset$sample_info)

twins2Physeq <- phyloseq(OTU, taxa, sampledata)
twins2Physeq1 <- prune_samples(sample_sums(twins2Physeq)>30, twins2Physeq)
twins2Physeq1 <- subset_samples(twins2Physeq1,host.age<50)

# Total reads per sample are sometimes modeled using Poisson or negative binomial. 
qplot(colSums(otu_table(twins2Physeq1)),binwidth=1000) +
  xlab("Total counts per sample")

# fit.sd <- goodfit(log(rowSums(raw.counts)+1), type='nbinomial')

# Show depth vs percent of zeros
raw.counts <- as(otu_table(twins2Physeq1), "matrix")
percent.zero <- apply(raw.counts, 2, function(a) mean(a==0)) 
df <- data.frame('Depth'=colSums(raw.counts),'Percent.Zero'=percent.zero) %>%
  ggplot(aes(x=Depth, y=Percent.Zero)) + geom_point() + geom_smooth() + xlab("Sequencing depth") + ylab("Percentage of zeros")

# Visualize prevalence
rank_names(twins2Physeq1)
# Create table, number of features for each phyla
table(tax_table(twins2Physeq1)[, "Phylum"], exclude = NULL)
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(twins2Physeq1),
               MARGIN = ifelse(taxa_are_rows(twins2Physeq1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(twins2Physeq1),
                    tax_table(twins2Physeq1))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(twins2Physeq1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Prune taxa based on prevalence
prevalenceThreshold = 0.1 * nsamples(twins2Physeq1)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
twins2Physeq2 = prune_taxa(keepTaxa, twins2Physeq1)
abundance_counts <- as(otu_table(twins2Physeq2), "matrix") + 1

## Check the distribution of the balance X_i over X_{A_i,j}
## Note y has heavy tails
# p <- nrow(abundance_counts)
# for (i in 1:(p-1)){
#   for (j in (i+1):p)
#   setA <- setdiff(1:p,c(i,j))
#   y <- log(abundance_counts[i,]) - colMeans(log(abundance_counts[setA,]))
#   plt <- qplot(y)
#   save_plot(paste0('figures/qplot_',i,'_',j,'.png'),plt)
# }
# X <- t(log(abundance_counts[setA,]))
# fit_classo <- cv.func('ConstrLasso',y=y,x=X,C=matrix(1,p-2,1), nlam = 100, nfolds=10)
# 
# residuals <- y - X %*% fit_classo$bet[,which.min(fit_classo$cvm)]
# qplot(residuals)

twins2Physeq3 <- transform_sample_counts(twins2Physeq2, function(OTU) OTU/sum(OTU) )
abundance_relative <- 100*as(otu_table(twins2Physeq3), "matrix") 
sample.info <- as(sample_data(twins2Physeq3),'data.frame')


# Use age as an outcome for prediction
# library(propr)
# W <- t(as(otu_table(twins2Physeq2), "matrix") )
# y <- sample.info$host.age
# p <- ncol(W)
# n <- nrow(W)
# "supervised" = { 
#   S <- matrix(0,p,p)
#   rownames(S) <- colnames(S) <- colnames(W)
#   for (j in 1:(p-1)){
#     for (k in (j+1):p){
#       index <- which(W[,j] * W[,k]>0)
#       newx <- W[,j]
#       newx[-index] <- 0
#       newx[index] <- log(W[index,j]) - log(W[index,k])
#       newx <- newx - mean(newx)
#       newy <- y - mean(y)
#       S[j,k] <- S[k,j] <- abs(cor(newx,y))
#       # S[j,k] <- S[k,j] <- abs(crossprod(newx,newy)/(sqrt(crossprod(newx)) * sqrt(crossprod(newy))))
#     }
#   }
#   h_supervised <- hclust(as.dist(1-S),method = linkage)
#   # plot(h_supervised)
#   # ggtree(h_supervised) + geom_point(aes(shape=isTip, color=isTip), size=3)
#   sbp_supervised <- sbp.fromHclust(h_supervised)
#   ba_supervised <- balance.fromSBP(WC,sbp_supervised)
#   fit_supervised <- run.glmnet(x=ba_supervised[1:n,],y,xt=ba_supervised[-(1:n),],y_test)
# }

# heatmap
# real.metaphlan.dist.bray <- vegdist(t(abundance_relative), method="bray") 
# 
# bk <- unique(c(
#   0,
#   min(min(abundance_relative[abundance_relative>0])/2,0.1),
#   seq(min(min(abundance_relative[abundance_relative>0])/2,1),0.1,length=10),
#   seq(0.1,1,length=10),
#   seq(1,5,length=10),
#   seq(5,20,length=20),
#   seq(20,50,length=20),
#   seq(50,100,length=20)
# ))
# hmcols<- colorRampPalette(c('white','#1f78b4','#ff7f00','#e41a1c'))(length(bk)-1)
# 
# ##########
# ### the color lables must be columns,so I have to transform the data
# pheatmap(abundance_relative,border_color="grey60",
#          clustering_distance_rows = 'correlation',
#          clustering_distance_cols = real.metaphlan.dist.bray,
#          clustering_method = "complete",
#          annotation=sample.info[,c("host.sex","zygosity")],
#          color = hmcols, breaks=bk,
#          show_rownames = TRUE,
#          scale="none"
# )
# 
# ## Apply PCA
# d <- t(apply(abundance_relative,2,clr))
# d_svd <- svd(d)
# 
# ds <- scale(d)
# pcar <- princomp(ds)
# fviz_eig(pcar, geom = "bar", bar_width = 0.4) + ggtitle("")
# 
# # eigen decomposition of a distance
# real.metaphlan.dist.bray <- vegdist(t(abundance_relative), method="bray") 
# # this distance matrix is not guaranteed to be PD
# db <- eigen(as.matrix(real.metaphlan.dist.bray))
# 

library(selbal)
x <- HIV[,1:60] + 1
y <- as.numeric(HIV[,62])

x <- sCD14[,1:60] + 1
y <- sCD14[,61]

p <- ncol(x)
rhoMat <- matrix(0,p,p)
rownames(rhoMat) <- colnames(rhoMat) <- colnames(x)
rhoMat2 <- rhoMat
testmat <- rhoMat
for (j in 1:p){
  for (k in 1:p){
    if (k==j){next}
    else {
      testmat[j,k] <- cor.test(log(x[,j])-log(x[,k]),y)$p.value
      rhoMat[j,k] <- abs(stats::cor(log(x[,j])-log(x[,k]),y))
    }
  }
}
diag(rhoMat) <- 0
library(pheatmap)

pheatmap(rhoMat,show_colnames = F, show_rownames = F,
         cluster_rows = T, cluster_cols = T,
         main = 'Heatmap of the original correlation')


b = spectral.clustering(rhoMat)
table(b)
index <- which(b==1)
which(b==-1)
source("RCode/slr-main.R")
