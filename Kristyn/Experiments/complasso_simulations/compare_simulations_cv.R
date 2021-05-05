library(ggplot2)
library(reshape2)
library(ggpubr)

output_dir = "Kristyn/Experiments/complasso_simulations/output"
rng.seed = 123
n = 100
p = 200
rho = 0.2 # 0.2, 0.5
intercept = TRUE
K = 10
beta.settings = "old"

# other stuff
metrics0 = c("PEtr", "PEte", "EA1", "EA2", "EAInfty", 
             "FP", "FN", "TPR", "betaSparsity", 
             "FPold", "FNold", "TPRold", "betaSparsityOld")
metrics = c("PEtr", "PEte", "EA1", "EA2", "EAInfty", 
            "FP", "FN", "TPR", "betaSparsity")
file.end = paste0(
  "_dim", n, "x", p, 
  "_rho", rho, 
  "_int", intercept,
  "_K", K,
  "_seed", rng.seed,
  ".rds")

################################################################################
# Compositional Lasso #
################################################################################
if(beta.settings == "old" | beta.settings == "linetal2014"){
  complasso.sims = readRDS(paste0(
    output_dir, "/complasso_cv_sims_old", file.end
  ))
  complasso.summaries = readRDS(paste0(
    output_dir, "/complasso_cv_summaries_old", file.end
  ))
} else{
  complasso.sims = readRDS(paste0(
    output_dir, "/complasso_cv_sims", file.end
  ))
  complasso.summaries = readRDS(paste0(
    output_dir, "/complasso_cv_summaries", file.end
  ))
}
print(complasso.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios #
################################################################################
if(beta.settings == "old" | beta.settings == "linetal2014"){
  slr.sims = readRDS(paste0(
    output_dir, "/slr_cv_sims_old", file.end
  ))
  slr.summaries = readRDS(paste0(
    output_dir, "/slr_cv_summaries_old", file.end
  ))
} else{
  slr.sims = readRDS(paste0(
    output_dir, "/slr_cv_sims", file.end
  ))
  slr.summaries = readRDS(paste0(
    output_dir, "/slr_cv_summaries", file.end
  ))
}
print(slr.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios 2 - alpha = 1 #
################################################################################
if(beta.settings == "old" | beta.settings == "linetal2014"){
  slr2a1.sims = readRDS(paste0(
    output_dir, "/slr2_alpha1_cv_sims_old", file.end
  ))
  slr2a1.summaries = readRDS(paste0(
    output_dir, "/slr2_alpha1_cv_summaries_old", file.end
  ))
} else{
  slr2a1.sims = readRDS(paste0(
    output_dir, "/slr2_alpha1_cv_sims", file.end
  ))
  slr2a1.summaries = readRDS(paste0(
    output_dir, "/slr2_alpha1_cv_summaries", file.end
  ))
}
print(slr2a1.summaries[metrics, c("mean", "se")])

################################################################################
# Supervised Log Ratios 2 - alpha = 1 #
################################################################################
if(beta.settings == "old" | beta.settings == "linetal2014"){
  slr2a0.5.sims = readRDS(paste0(
    output_dir, "/slr2_alpha0.5_cv_sims_old", file.end
  ))
  slr2a0.5.summaries = readRDS(paste0(
    output_dir, "/slr2_alpha0.5_cv_summaries_old", file.end
  ))
} else{
  slr2a0.5.sims = readRDS(paste0(
    output_dir, "/slr2_alpha0.5_cv_sims", file.end
  ))
  slr2a0.5.summaries = readRDS(paste0(
    output_dir, "/slr2_alpha0.5_cv_summaries", file.end
  ))
}
print(slr2a0.5.summaries[metrics, c("mean", "se")])


# plot

cl.sims.gg = melt(data.frame(t(complasso.sims)))
cl.sims.gg$type = "CompLasso"
slr.sims.gg = melt(data.frame(t(slr.sims)))
slr.sims.gg$type = "SLR"
slr2a1.sims.gg = melt(data.frame(t(slr2a1.sims)))
slr2a1.sims.gg$type = "SLR2alpha1"
slr2a0.5.sims.gg = melt(data.frame(t(slr2a0.5.sims)))
slr2a0.5.sims.gg$type = "SLR2alpha0.5"
data.gg = rbind(cl.sims.gg, slr.sims.gg, slr2a1.sims.gg, slr2a0.5.sims.gg)
data.gg = dplyr::filter(data.gg, variable %in% metrics)
data.gg$type = factor(data.gg$type, levels = c("CompLasso", "SLR", "SLR2alpha1", "SLR2alpha0.5"))
ggplot(data.gg, aes(x = type, y = value, color = type)) + 
  facet_wrap(vars(variable), scales = "free_y") + 
  geom_boxplot() + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_blank())

# zoom in to PEtr, PEte
# plot PEtr
PEtr.gg = data.gg[data.gg$variable == "PEtr", ]
# PEtr.gg$value = log(PEtr.gg$value)
plt.PEtr = ggplot(PEtr.gg, aes(x = type, y = value, color = type)) + 
  geom_boxplot() +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")
# plot PEte
PEte.gg = data.gg[data.gg$variable == "PEte", ]
# PEte.gg$value = log(PEte.gg$value)
plt.PEte = ggplot(PEte.gg, aes(x = type, y = value, color = type)) + 
  geom_boxplot() +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
               geom = "errorbar", width = 0.75, 
               linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", shape = 17, size = 2, 
               color = "red") +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")
# plot both
ggarrange(plt.PEtr, plt.PEte)

