### Code to make figure plotting accuracy metrics across HSP tools on each 16S validation dataset.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories (and save wilcoxon output to RDS).

rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# HMP:
hmp_ko_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "hmp_ko_acc_metrics.rds",
                                                                    dataset_name = "HMP",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

hmp_ko_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "hmp_ko_acc_metrics.rds",
                                                                 dataset_name = "HMP",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)


hmp_ko_precision <- hmp_ko_precision_outlist[[1]]
hmp_ko_precision_wilcoxon <- hmp_ko_precision_outlist[[2]]

hmp_ko_recall <- hmp_ko_recall_outlist[[1]]
hmp_ko_recall_wilcoxon <- hmp_ko_recall_outlist[[2]]


# Mammal:

mammal_ko_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "mammal_ko_acc_metrics.rds",
                                                                    dataset_name = "Mammal",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

mammal_ko_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "mammal_ko_acc_metrics.rds",
                                                                 dataset_name = "Mammal",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)



mammal_ko_precision <- mammal_ko_precision_outlist[[1]]
mammal_ko_precision_wilcoxon <- mammal_ko_precision_outlist[[2]]

mammal_ko_recall <- mammal_ko_recall_outlist[[1]]
mammal_ko_recall_wilcoxon <- mammal_ko_recall_outlist[[2]]


# Ocean:

ocean_ko_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "ocean_ko_acc_metrics.rds",
                                                                    dataset_name = "Ocean",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

ocean_ko_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "ocean_ko_acc_metrics.rds",
                                                                 dataset_name = "Ocean",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)

ocean_ko_precision <- ocean_ko_precision_outlist[[1]]
ocean_ko_precision_wilcoxon <- ocean_ko_precision_outlist[[2]]

ocean_ko_recall <- ocean_ko_recall_outlist[[1]]
ocean_ko_recall_wilcoxon <- ocean_ko_recall_outlist[[2]]



# Soil (Blueberry):

blueberry_ko_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "blueberry_ko_acc_metrics.rds",
                                                                    dataset_name = "Soil (Blueberry)",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

blueberry_ko_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "blueberry_ko_acc_metrics.rds",
                                                                 dataset_name = "Soil (Blueberry)",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)

blueberry_ko_precision <- blueberry_ko_precision_outlist[[1]]
blueberry_ko_precision_wilcoxon <- blueberry_ko_precision_outlist[[2]]

blueberry_ko_recall <- blueberry_ko_recall_outlist[[1]]
blueberry_ko_recall_wilcoxon <- blueberry_ko_recall_outlist[[2]]


# Precision
combined_ko_precision <- rbind(hmp_ko_precision, mammal_ko_precision,
                         ocean_ko_precision, blueberry_ko_precision)

combined_ko_precision_wilcoxon <- rbind(hmp_ko_precision_wilcoxon, mammal_ko_precision_wilcoxon,
                                  ocean_ko_precision_wilcoxon, blueberry_ko_precision_wilcoxon)

combined_ko_precision$cat <- as.character(combined_ko_precision$cat)
combined_ko_precision$cat <- factor(combined_ko_precision$cat,
                                             levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ko_precision_melt <- melt(combined_ko_precision)

ko_precision_boxplots <- ggplot(combined_ko_precision_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.2, 1)) +
  ylab(c("Precision")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


# Recall
combined_ko_recall <- rbind(hmp_ko_recall, mammal_ko_recall,
                               ocean_ko_recall, blueberry_ko_recall)

combined_ko_recall_wilcoxon <- rbind(hmp_ko_recall_wilcoxon, mammal_ko_recall_wilcoxon,
                                        ocean_ko_recall_wilcoxon, blueberry_ko_recall_wilcoxon)

combined_ko_recall$cat <- as.character(combined_ko_recall$cat)
combined_ko_recall$cat <- factor(combined_ko_recall$cat,
                                    levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2",
                                             "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ko_recall_melt <- melt(combined_ko_recall)

ko_recall_boxplots <- ggplot(combined_ko_recall_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.2, 1)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))



# Plot boxplots:
plot_grid(ko_precision_boxplots,
          ko_recall_boxplots,
          nrow=2,
          ncol=1,
          labels=c("A", "B"))


# Report mean and sd per dataset.

# For precision:
mean(hmp_ko_precision[which(hmp_ko_precision$cat == "NSTI=2"), "precision"])
sd(hmp_ko_precision[which(hmp_ko_precision$cat == "NSTI=2"), "precision"])

mean(mammal_ko_precision[which(mammal_ko_precision$cat == "NSTI=2"), "precision"])
sd(mammal_ko_precision[which(mammal_ko_precision$cat == "NSTI=2"), "precision"])

mean(ocean_ko_precision[which(ocean_ko_precision$cat == "NSTI=2"), "precision"])
sd(ocean_ko_precision[which(ocean_ko_precision$cat == "NSTI=2"), "precision"])

mean(blueberry_ko_precision[which(blueberry_ko_precision$cat == "NSTI=2"), "precision"])
sd(blueberry_ko_precision[which(blueberry_ko_precision$cat == "NSTI=2"), "precision"])

