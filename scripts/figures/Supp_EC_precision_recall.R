### Code to make figure plotting accuracy metrics across HSP tools on each 16S validation dataset.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories (and save wilcoxon output to RDS).

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# HMP:
hmp_ec_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "hmp_ec_acc_metrics.rds",
                                                                    dataset_name = "HMP",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

hmp_ec_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "hmp_ec_acc_metrics.rds",
                                                                 dataset_name = "HMP",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)


hmp_ec_precision <- hmp_ec_precision_outlist[[1]]
hmp_ec_precision_wilcoxon <- hmp_ec_precision_outlist[[2]]

hmp_ec_recall <- hmp_ec_recall_outlist[[1]]
hmp_ec_recall_wilcoxon <- hmp_ec_recall_outlist[[2]]


# Mammal:

mammal_ec_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "mammal_ec_acc_metrics.rds",
                                                                    dataset_name = "Mammal",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

mammal_ec_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "mammal_ec_acc_metrics.rds",
                                                                 dataset_name = "Mammal",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)



mammal_ec_precision <- mammal_ec_precision_outlist[[1]]
mammal_ec_precision_wilcoxon <- mammal_ec_precision_outlist[[2]]

mammal_ec_recall <- mammal_ec_recall_outlist[[1]]
mammal_ec_recall_wilcoxon <- mammal_ec_recall_outlist[[2]]


# Ocean:

ocean_ec_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "ocean_ec_acc_metrics.rds",
                                                                    dataset_name = "Ocean",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

ocean_ec_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "ocean_ec_acc_metrics.rds",
                                                                 dataset_name = "Ocean",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)

ocean_ec_precision <- ocean_ec_precision_outlist[[1]]
ocean_ec_precision_wilcoxon <- ocean_ec_precision_outlist[[2]]

ocean_ec_recall <- ocean_ec_recall_outlist[[1]]
ocean_ec_recall_wilcoxon <- ocean_ec_recall_outlist[[2]]



# Soil (Blueberry):

blueberry_ec_precision_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "blueberry_ec_acc_metrics.rds",
                                                                    dataset_name = "Soil (Blueberry)",
                                                                    metric_col="precision",
                                                                    wilcox_cat2ignore = extra_nsti_categories,
                                                                    y_pos_start = 1.03)

blueberry_ec_recall_outlist <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = "blueberry_ec_acc_metrics.rds",
                                                                 dataset_name = "Soil (Blueberry)",
                                                                 metric_col="recall",
                                                                 wilcox_cat2ignore = extra_nsti_categories,
                                                                 y_pos_start = 1.03)

blueberry_ec_precision <- blueberry_ec_precision_outlist[[1]]
blueberry_ec_precision_wilcoxon <- blueberry_ec_precision_outlist[[2]]

blueberry_ec_recall <- blueberry_ec_recall_outlist[[1]]
blueberry_ec_recall_wilcoxon <- blueberry_ec_recall_outlist[[2]]


# Precision
combined_ec_precision <- rbind(hmp_ec_precision, mammal_ec_precision,
                         ocean_ec_precision, blueberry_ec_precision)

combined_ec_precision_wilcoxon <- rbind(hmp_ec_precision_wilcoxon, mammal_ec_precision_wilcoxon,
                                  ocean_ec_precision_wilcoxon, blueberry_ec_precision_wilcoxon)

combined_ec_precision$cat <- as.character(combined_ec_precision$cat)
combined_ec_precision$cat <- factor(combined_ec_precision$cat,
                                    levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                             "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ec_precision_melt <- melt(combined_ec_precision)

ec_precision_boxplots <- ggplot(combined_ec_precision_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
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
combined_ec_recall <- rbind(hmp_ec_recall, mammal_ec_recall,
                               ocean_ec_recall, blueberry_ec_recall)

combined_ec_recall_wilcoxon <- rbind(hmp_ec_recall_wilcoxon, mammal_ec_recall_wilcoxon,
                                        ocean_ec_recall_wilcoxon, blueberry_ec_recall_wilcoxon)

combined_ec_recall$cat <- as.character(combined_ec_recall$cat)
combined_ec_recall$cat <- factor(combined_ec_recall$cat,
                                 levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                          "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ec_recall_melt <- melt(combined_ec_recall)

ec_recall_boxplots <- ggplot(combined_ec_recall_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.2, 1)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))



# 11 x 7
plot_grid(ec_precision_boxplots,
          ec_recall_boxplots,
          nrow=2,
          ncol=1,
          labels=c("A", "B"))


# Report mean and sd per dataset.

# For precision:
mean(hmp_ec_precision[which(hmp_ec_precision$cat == "NSTI=2"), "precision"])
sd(hmp_ec_precision[which(hmp_ec_precision$cat == "NSTI=2"), "precision"])

mean(mammal_ec_precision[which(mammal_ec_precision$cat == "NSTI=2"), "precision"])
sd(mammal_ec_precision[which(mammal_ec_precision$cat == "NSTI=2"), "precision"])

mean(ocean_ec_precision[which(ocean_ec_precision$cat == "NSTI=2"), "precision"])
sd(ocean_ec_precision[which(ocean_ec_precision$cat == "NSTI=2"), "precision"])

mean(blueberry_ec_precision[which(blueberry_ec_precision$cat == "NSTI=2"), "precision"])
sd(blueberry_ec_precision[which(blueberry_ec_precision$cat == "NSTI=2"), "precision"])

