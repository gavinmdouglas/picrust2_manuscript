### Code to make figure plotting accuracy metrics across HSP tools on each 16S validation dataset.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories (and save wilcoxon output to RDS).

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# HMP:
hmp_pathabun_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "hmp_pathabun_acc_metrics.rds",
                                                                    dataset_name = "HMP",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

hmp_pathabun_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "hmp_pathabun_acc_metrics.rds",
                                                                 dataset_name = "HMP",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)


hmp_pathabun_precision <- hmp_pathabun_precision_outlist[[1]]
hmp_pathabun_precision_wilcoxon <- hmp_pathabun_precision_outlist[[2]]

hmp_pathabun_recall <- hmp_pathabun_recall_outlist[[1]]
hmp_pathabun_recall_wilcoxon <- hmp_pathabun_recall_outlist[[2]]


# Mammal:

mammal_pathabun_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "mammal_pathabun_acc_metrics.rds",
                                                                    dataset_name = "Mammal",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

mammal_pathabun_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "mammal_pathabun_acc_metrics.rds",
                                                                 dataset_name = "Mammal",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)



mammal_pathabun_precision <- mammal_pathabun_precision_outlist[[1]]
mammal_pathabun_precision_wilcoxon <- mammal_pathabun_precision_outlist[[2]]

mammal_pathabun_recall <- mammal_pathabun_recall_outlist[[1]]
mammal_pathabun_recall_wilcoxon <- mammal_pathabun_recall_outlist[[2]]


# Ocean:

ocean_pathabun_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "ocean_pathabun_acc_metrics.rds",
                                                                    dataset_name = "Ocean",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

ocean_pathabun_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "ocean_pathabun_acc_metrics.rds",
                                                                 dataset_name = "Ocean",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)

ocean_pathabun_precision <- ocean_pathabun_precision_outlist[[1]]
ocean_pathabun_precision_wilcoxon <- ocean_pathabun_precision_outlist[[2]]

ocean_pathabun_recall <- ocean_pathabun_recall_outlist[[1]]
ocean_pathabun_recall_wilcoxon <- ocean_pathabun_recall_outlist[[2]]



# Soil (Blueberry):

blueberry_pathabun_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "blueberry_pathabun_acc_metrics.rds",
                                                                    dataset_name = "Soil (Blueberry)",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

blueberry_pathabun_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "blueberry_pathabun_acc_metrics.rds",
                                                                 dataset_name = "Soil (Blueberry)",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)

blueberry_pathabun_precision <- blueberry_pathabun_precision_outlist[[1]]
blueberry_pathabun_precision_wilcoxon <- blueberry_pathabun_precision_outlist[[2]]

blueberry_pathabun_recall <- blueberry_pathabun_recall_outlist[[1]]
blueberry_pathabun_recall_wilcoxon <- blueberry_pathabun_recall_outlist[[2]]


# Precision
combined_pathabun_precision <- rbind(hmp_pathabun_precision, mammal_pathabun_precision,
                         ocean_pathabun_precision, blueberry_pathabun_precision)

combined_pathabun_precision_wilcoxon <- rbind(hmp_pathabun_precision_wilcoxon, mammal_pathabun_precision_wilcoxon,
                                  ocean_pathabun_precision_wilcoxon, blueberry_pathabun_precision_wilcoxon)

combined_pathabun_precision$cat <- as.character(combined_pathabun_precision$cat)
combined_pathabun_precision$cat <- factor(combined_pathabun_precision$cat,
                                    levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                             "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_pathabun_precision_melt <- melt(combined_pathabun_precision)

pathabun_precision_boxplots <- ggplot(combined_pathabun_precision_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.2, 1)) +
  ylab(c("Precision")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))


# Recall
combined_pathabun_recall <- rbind(hmp_pathabun_recall, mammal_pathabun_recall,
                               ocean_pathabun_recall, blueberry_pathabun_recall)

combined_pathabun_recall_wilcoxon <- rbind(hmp_pathabun_recall_wilcoxon, mammal_pathabun_recall_wilcoxon,
                                        ocean_pathabun_recall_wilcoxon, blueberry_pathabun_recall_wilcoxon)

combined_pathabun_recall$cat <- as.character(combined_pathabun_recall$cat)
combined_pathabun_recall$cat <- factor(combined_pathabun_recall$cat,
                                 levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                          "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_pathabun_recall_melt <- melt(combined_pathabun_recall)

pathabun_recall_boxplots <- ggplot(combined_pathabun_recall_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.2, 1)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))



# Plot boxplots (9x7):
plot_grid(pathabun_precision_boxplots,
          pathabun_recall_boxplots,
          nrow=2,
          ncol=1,
          labels=c("A", "B"))


# Report mean and sd per dataset.

# For precision:
mean(hmp_pathabun_precision[which(hmp_pathabun_precision$cat == "NSTI=2"), "precision"])
sd(hmp_pathabun_precision[which(hmp_pathabun_precision$cat == "NSTI=2"), "precision"])

mean(mammal_pathabun_precision[which(mammal_pathabun_precision$cat == "NSTI=2"), "precision"])
sd(mammal_pathabun_precision[which(mammal_pathabun_precision$cat == "NSTI=2"), "precision"])

mean(ocean_pathabun_precision[which(ocean_pathabun_precision$cat == "NSTI=2"), "precision"])
sd(ocean_pathabun_precision[which(ocean_pathabun_precision$cat == "NSTI=2"), "precision"])

mean(blueberry_pathabun_precision[which(blueberry_pathabun_precision$cat == "NSTI=2"), "precision"])
sd(blueberry_pathabun_precision[which(blueberry_pathabun_precision$cat == "NSTI=2"), "precision"])

