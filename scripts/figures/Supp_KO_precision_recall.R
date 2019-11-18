### Code to make figure plotting accuracy metrics across HSP tools on each 16S validation dataset.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories (and save wilcoxon output to RDS).

rm(list=ls(all=TRUE))

library(cowplot)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

ko_precision_raw <- list()
ko_precision <- list()
ko_precision_wilcoxon <- list()

ko_recall_raw <- list()
ko_recall <- list()
ko_recall_wilcoxon <- list()

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroonian", "indian"="Indian", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")

for(d in datasets) {
  
  ko_precision_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_ko_acc_df.rds", sep=""),
                                                              dataset_name = dataset2name[[d]],
                                                              metric_col="precision",
                                                              wilcox_cat2ignore = extra_nsti_categories,
                                                              y_pos_start = 1.03)
  
  ko_precision[[d]] <- ko_precision_raw[[d]][[1]]
  ko_precision_wilcoxon[[d]] <- ko_precision_raw[[d]][[2]]
  
  
  ko_recall_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_ko_acc_df.rds", sep=""),
                                                           dataset_name = dataset2name[[d]],
                                                           metric_col="recall",
                                                           wilcox_cat2ignore = extra_nsti_categories,
                                                           y_pos_start = 1.03)
  
  ko_recall[[d]] <- ko_recall_raw[[d]][[1]]
  ko_recall_wilcoxon[[d]] <- ko_recall_raw[[d]][[2]]
  
}

# Precision
combined_ko_precision <- do.call("rbind", ko_precision)

combined_ko_precision_wilcoxon <- do.call("rbind", ko_precision_wilcoxon)

combined_ko_precision$cat <- as.character(combined_ko_precision$cat)

# Remove some NSTI categories so it's easier to plot.
combined_ko_precision <- combined_ko_precision[-which(combined_ko_precision$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]

combined_ko_precision$cat <- factor(combined_ko_precision$cat,
                                             levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ko_precision_melt <- melt(combined_ko_precision)

combined_ko_precision_melt$dataset <- factor(combined_ko_precision_melt$dataset, levels=c("Cameroonian", "HMP", "Indian", "Mammal", "Ocean", "Primate", "Soil (Blueberry)"))


ko_precision_boxplots <- ggplot(combined_ko_precision_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  scale_y_continuous(breaks=c(0.4, 0.6, 0.8, 1.0), limits=c(0.25, 1.30)) +
  ylab(c("Precision")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(combined_ko_precision_wilcoxon, label = "p_symbol")


# Recall
combined_ko_recall <- do.call("rbind", ko_recall)

combined_ko_recall_wilcoxon <- do.call("rbind", ko_recall_wilcoxon)

combined_ko_recall$cat <- as.character(combined_ko_recall$cat)
combined_ko_recall <- combined_ko_recall[-which(combined_ko_recall$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]



combined_ko_recall$cat <- factor(combined_ko_recall$cat,
                                    levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2",
                                             "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ko_recall_melt <- melt(combined_ko_recall)

combined_ko_recall_melt$dataset <- factor(combined_ko_recall_melt$dataset, levels=c("Cameroonian", "HMP", "Indian", "Mammal", "Ocean", "Primate", "Soil (Blueberry)"))


ko_recall_boxplots <- ggplot(combined_ko_recall_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  scale_y_continuous(breaks=c(0.4, 0.6, 0.8, 1.0), limits=c(0.25, 1.30)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(combined_ko_recall_wilcoxon, label = "p_symbol")

pdf(file = "../../../figures/Supp_KO_precision_recall.pdf", width=16, height=10)

plot_grid(ko_precision_boxplots,
          ko_recall_boxplots,
          nrow=2,
          ncol=1,
          labels=c("a", "b"))
dev.off()



# Calculate basic summary statistics.

dataset_ko_acc_stats <- data.frame(matrix(NA, nrow=7, ncol=4))
colnames(dataset_ko_acc_stats) <- c("precision_PICRUSt2_mean", "precision_PICRUSt2_sd", "recall_PICRUSt2_mean", "recall_PICRUSt2_sd")
rownames(dataset_ko_acc_stats) <- datasets

for(d in datasets) {
  
  d_name <- dataset2name[[d]]
  
  dataset_ko_acc_stats[d, ] <- c(mean(combined_ko_precision[which(combined_ko_precision$dataset == d_name & combined_ko_precision$cat == "NSTI=2"), "precision"]),
                                 sd(combined_ko_precision[which(combined_ko_precision$dataset == d_name & combined_ko_precision$cat == "NSTI=2"), "precision"]),
                                 mean(combined_ko_recall[which(combined_ko_recall$dataset == d_name & combined_ko_recall$cat == "NSTI=2"), "recall"]),
                                 sd(combined_ko_recall[which(combined_ko_recall$dataset == d_name & combined_ko_recall$cat == "NSTI=2"), "recall"]))
  
}

