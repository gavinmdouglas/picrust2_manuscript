### Code for making PAPRICA pathabun boxplots.

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroon", "indian"="India", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")

paprica_pathabun_metrics <- readRDS("paprica_path_metrics.rds")

paprica_pathabun_rho_list <- list()
paprica_pathabun_acc_list <- list()

for(d in datasets) {
 
  paprica_pathabun_metrics[[d]]$spearman$dataset <- dataset2name[[d]]
  paprica_pathabun_metrics[[d]]$acc$dataset <- dataset2name[[d]]
  
  paprica_pathabun_rho_list[[d]] <- paprica_pathabun_metrics[[d]]$spearman
  paprica_pathabun_acc_list[[d]] <- paprica_pathabun_metrics[[d]]$acc
}

combined_pathabun_paprica_rho <- do.call("rbind", paprica_pathabun_rho_list)
combined_pathabun_paprica_acc <- do.call("rbind", paprica_pathabun_acc_list)

combined_pathabun_paprica_rho_melt <- melt(combined_pathabun_paprica_rho)
combined_pathabun_paprica_acc_melt <- melt(combined_pathabun_paprica_acc)

combined_pathabun_paprica_acc_melt_precision <- combined_pathabun_paprica_acc_melt[which(combined_pathabun_paprica_acc_melt$variable == "precision"), ]
combined_pathabun_paprica_acc_melt_recall <- combined_pathabun_paprica_acc_melt[which(combined_pathabun_paprica_acc_melt$variable == "recall"), ]

combined_pathabun_paprica_rho_plot <- ggplot(combined_pathabun_paprica_rho_melt, aes(x=cat, y=value, fill=cat)) +
                                            geom_boxplot(outlier.shape = NA) +
                                            geom_quasirandom(size=0.1, dodge.width=0.8) +
                                            ylim(c(0.5, 1)) +
                                            ylab(c("Spearman Correlation Coefficient")) +
                                            xlab("") +
                                            guides(fill=FALSE) +
                                            facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                            theme(panel.background = element_rect(fill = "gray90"), axis.line = element_line(colour = "black"),
                                                  axis.text.x=element_text(angle=45, hjust=1)) +
                                            scale_fill_manual(values=c("light grey", "#F8766D"))

combined_pathabun_paprica_precision_plot <- ggplot(combined_pathabun_paprica_acc_melt_precision, aes(x=category, y=value, fill=category)) +
                                                    geom_boxplot(outlier.shape = NA) +
                                                    geom_quasirandom(size=0.1, dodge.width=0.8) +
                                                    ylim(c(0.5, 1)) +
                                                    ylab(c("Precision")) +
                                                    xlab("") +
                                                    guides(fill=FALSE) +
                                                    facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                                    theme(panel.background = element_rect(fill = "gray90"), axis.line = element_line(colour = "black"),
                                                          axis.text.x=element_text(angle=45, hjust=1)) +
                                                    scale_fill_manual(values=c("light grey", "#F8766D"))

combined_pathabun_paprica_recall_plot <- ggplot(combined_pathabun_paprica_acc_melt_recall, aes(x=category, y=value, fill=category)) +
                                                geom_boxplot(outlier.shape = NA) +
                                                geom_quasirandom(size=0.1, dodge.width=0.8) +
                                                ylim(c(0.5, 1)) +
                                                ylab(c("Recall")) +
                                                xlab("") +
                                                guides(fill=FALSE) +
                                                facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                                theme(panel.background = element_rect(fill = "gray90"), axis.line = element_line(colour = "black"),
                                                      axis.text.x=element_text(angle=45, hjust=1)) +
                                                scale_fill_manual(values=c("light grey", "#F8766D"))

pdf(file = "../../../figures/Supp_pathabun_paprica.pdf", width=14, height=13)

plot_grid(combined_pathabun_paprica_rho_plot,
          combined_pathabun_paprica_precision_plot,
          combined_pathabun_paprica_recall_plot,
          labels = c("a", "b", "c"),
          ncol=1, nrow=3)

dev.off()
