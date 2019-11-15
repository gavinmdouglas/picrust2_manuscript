
rm(list=ls(all=TRUE))

library(cowplot)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

ec_precision_raw <- list()
ec_precision <- list()
ec_precision_wilcoxon <- list()

ec_recall_raw <- list()
ec_recall <- list()
ec_recall_wilcoxon <- list()

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroon", "indian"="India", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")

for(d in datasets) {
  
  ec_precision_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_ec_acc_df.rds", sep=""),
                                                                   dataset_name = dataset2name[[d]],
                                                                   metric_col="precision",
                                                                   wilcox_cat2ignore = extra_nsti_categories,
                                                                   y_pos_start = 1.03)
  
  ec_precision[[d]] <- ec_precision_raw[[d]][[1]]
  ec_precision_wilcoxon[[d]] <- ec_precision_raw[[d]][[2]]
  
  
  ec_recall_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_ec_acc_df.rds", sep=""),
                                                                dataset_name = dataset2name[[d]],
                                                                metric_col="recall",
                                                                wilcox_cat2ignore = extra_nsti_categories,
                                                                y_pos_start = 1.03)
  
  ec_recall[[d]] <- ec_recall_raw[[d]][[1]]
  ec_recall_wilcoxon[[d]] <- ec_recall_raw[[d]][[2]]
  
}

# Precision
combined_ec_precision <- do.call("rbind", ec_precision)

combined_ec_precision_wilcoxon <- do.call("rbind", ec_precision_wilcoxon)

combined_ec_precision$cat <- as.character(combined_ec_precision$cat)

# Remove some NSTI categories so it's easier to plot.
combined_ec_precision <- combined_ec_precision[-which(combined_ec_precision$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]

combined_ec_precision$cat <- factor(combined_ec_precision$cat,
                                    levels=c("Null", "PAPRICA", "NSTI=2","NSTI=1","NSTI=0.05"))
combined_ec_precision_melt <- melt(combined_ec_precision)

combined_ec_precision_melt$dataset <- factor(combined_ec_precision_melt$dataset, levels=c("Cameroon", "Indian", "HMP", "Primate", "Mammal", "Ocean", "Soil (Blueberry)"))


ec_precision_boxplots <- ggplot(combined_ec_precision_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.2, 1)) +
  ylab(c("Precision")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


# Recall
combined_ec_recall <- do.call("rbind", ec_recall)

combined_ec_recall_wilcoxon <- do.call("rbind", ec_recall_wilcoxon)

combined_ec_recall$cat <- as.character(combined_ec_recall$cat)
combined_ec_recall <- combined_ec_recall[-which(combined_ec_recall$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]



combined_ec_recall$cat <- factor(combined_ec_recall$cat,
                                 levels=c("Null", "PAPRICA", "NSTI=2",
                                          "NSTI=1", "NSTI=0.05"))
combined_ec_recall_melt <- melt(combined_ec_recall)

combined_ec_recall_melt$dataset <- factor(combined_ec_recall_melt$dataset, levels=c("Cameroon", "Indian", "HMP", "Primate", "Mammal", "Ocean", "Soil (Blueberry)"))


ec_recall_boxplots <- ggplot(combined_ec_recall_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.2, 1)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

pdf(file = "../../../figures/Supp_EC_precision_recall.pdf", width=14, height=8)

plot_grid(ec_precision_boxplots,
          ec_recall_boxplots,
          nrow=2,
          ncol=1,
          labels=c("a", "b"))
dev.off()
