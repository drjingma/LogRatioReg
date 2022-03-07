# looking at the genera, making plots and stuff

# data
# 98 samples, 87 genera
# replace zero counts with 0.5 (maximum rounding error)
load(paste0("Data/", "BMI.rda"))
# dim(raw_data) # 98 x 89
# dim(X) # 98 x 87
# dim(X.prop) # 98 x 87
log.X.prop = log(X.prop)
n = dim(X)[1]
num.genera = dim(X)[2]

prop_zero_per_sample <- apply(X, 1, function(a) sum(a<1)/ncol(X))
prop_zero_per_otu <- apply(X, 2, function(a) sum(a<1)/nrow(X))

dim(X.prop)
top = c(
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Acidaminococcus", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Allisonella", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Clostridiaceae.Clostridium", 
  "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Rikenellaceae.Alistipes", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Oscillibacter", 
  "Bacteria.Proteobacteria.Gammaproteobacteria.Xanthomonadales.Xanthomonadaceae.Stenotrophomonas",
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Veillonellaceae.Megasphaera", 
  "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea"
  )
X.prop.top = X.prop[, top]
colnames(X.prop.top) = c("A", "B", "C", "D", "E", "F", "G", "H")
library(reshape2)
X.prop.top.gg = melt(as.data.frame(X.prop.top))
ggplot(data = X.prop.top.gg, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = variable))
ggplot(data = X.prop.top.gg, aes(x = variable, y = log(value))) + 
  geom_boxplot(aes(fill = variable))
