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

pathabun_precision_raw <- list()
pathabun_precision <- list()
pathabun_precision_wilcoxon <- list()

pathabun_recall_raw <- list()
pathabun_recall <- list()
pathabun_recall_wilcoxon <- list()

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroon", "indian"="India", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")

for(d in datasets) {
  
  pathabun_precision_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_pathabun_acc_df.rds", sep=""),
                                                                   dataset_name = dataset2name[[d]],
                                                                   metric_col="precision",
                                                                   wilcox_cat2ignore = extra_nsti_categories,
                                                                   y_pos_start = 1.03)
  
  pathabun_precision[[d]] <- pathabun_precision_raw[[d]][[1]]
  pathabun_precision_wilcoxon[[d]] <- pathabun_precision_raw[[d]][[2]]
  
  
  pathabun_recall_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_pathabun_acc_df.rds", sep=""),
                                                                dataset_name = dataset2name[[d]],
                                                                metric_col="recall",
                                                                wilcox_cat2ignore = extra_nsti_categories,
                                                                y_pos_start = 1.03)
  
  pathabun_recall[[d]] <- pathabun_recall_raw[[d]][[1]]
  pathabun_recall_wilcoxon[[d]] <- pathabun_recall_raw[[d]][[2]]
  
}

# Precision
combined_pathabun_precision <- do.call("rbind", pathabun_precision)

combined_pathabun_precision_wilcoxon <- do.call("rbind", pathabun_precision_wilcoxon)

combined_pathabun_precision$cat <- as.character(combined_pathabun_precision$cat)

# Remove some NSTI categories so it's easier to plot.
combined_pathabun_precision <- combined_pathabun_precision[-which(combined_pathabun_precision$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]

combined_pathabun_precision$cat <- factor(combined_pathabun_precision$cat,
                                    levels=c("Null", "NSTI=2", "NSTI=1", "NSTI=0.05"))
combined_pathabun_precision_melt <- melt(combined_pathabun_precision)

combined_pathabun_precision_melt$dataset <- factor(combined_pathabun_precision_melt$dataset, levels=c("Cameroon", "Indian", "HMP", "Primate", "Mammal", "Ocean", "Soil (Blueberry)"))


pathabun_precision_boxplots <- ggplot(combined_pathabun_precision_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.5, 1)) +
  ylab(c("Precision")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))


# Recall
combined_pathabun_recall <- do.call("rbind", pathabun_recall)

combined_pathabun_recall_wilcoxon <- do.call("rbind", pathabun_recall_wilcoxon)

combined_pathabun_recall$cat <- as.character(combined_pathabun_recall$cat)
combined_pathabun_recall <- combined_pathabun_recall[-which(combined_pathabun_recall$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]



combined_pathabun_recall$cat <- factor(combined_pathabun_recall$cat,
                                 levels=c("Null", "NSTI=2", "NSTI=1", "NSTI=0.05"))
combined_pathabun_recall_melt <- melt(combined_pathabun_recall)

combined_pathabun_recall_melt$dataset <- factor(combined_pathabun_recall_melt$dataset, levels=c("Cameroon", "Indian", "HMP", "Primate", "Mammal", "Ocean", "Soil (Blueberry)"))


pathabun_recall_boxplots <- ggplot(combined_pathabun_recall_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.5, 1)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))

pdf(file = "../../../figures/Supp_pathabun_precision_recall.pdf", width=13, height=8)

plot_grid(pathabun_precision_boxplots,
          pathabun_recall_boxplots,
          nrow=2,
          ncol=1,
          labels=c("a", "b"))
dev.off()
