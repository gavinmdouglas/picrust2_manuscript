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

ko_f1_raw <- list()
ko_f1 <- list()
ko_f1_wilcoxon <- list()

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroonian", "indian"="Indian", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")

for(d in datasets) {
  
  ko_f1_raw[[d]] <- parse_acc_metrics_rds_and_calc_wilcoxon(acc_rds = paste(d, "_ko_acc_df.rds", sep=""),
                                                              dataset_name = dataset2name[[d]],
                                                              metric_col="F1",
                                                              wilcox_cat2ignore = extra_nsti_categories,
                                                              y_pos_start = 1.03)
  
  ko_f1[[d]] <- ko_f1_raw[[d]][[1]]
  ko_f1_wilcoxon[[d]] <- ko_f1_raw[[d]][[2]]
  
  
}

# Precision
combined_ko_f1 <- do.call("rbind", ko_f1)

combined_ko_f1_wilcoxon <- do.call("rbind", ko_f1_wilcoxon)

combined_ko_f1$cat <- as.character(combined_ko_f1$cat)

# Remove some NSTI categories so it's easier to plot.
combined_ko_f1 <- combined_ko_f1[-which(combined_ko_f1$cat %in% c("NSTI=1.5", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]

combined_ko_f1$cat <- factor(combined_ko_f1$cat,
                                             levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
combined_ko_f1_melt <- melt(combined_ko_f1)

combined_ko_f1_melt$dataset <- factor(combined_ko_f1_melt$dataset, levels=c("Cameroonian", "HMP", "Indian", "Mammal", "Ocean", "Primate", "Soil (Blueberry)"))


pdf(file = "../../../figures/Supp_KO_F1.pdf", width=16, height=8)

ggplot(combined_ko_f1_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  scale_y_continuous(breaks=c(0.4, 0.6, 0.8, 1.0), limits=c(0.25, 1.30)) +
  ylab(c("F1 Score")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.background = element_rect(fill = "gray90"),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(combined_ko_f1_wilcoxon, label = "p_symbol")

dev.off()


        
        