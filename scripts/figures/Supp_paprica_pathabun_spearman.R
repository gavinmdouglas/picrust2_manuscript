### Code for making PAPRICA pathabun spearman boxplots.

rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

hmp_pathabun_paprica_rho <- readRDS("hmp_pathabun_PAPRICA_spearman_df.rds")
hmp_pathabun_paprica_acc <- readRDS("hmp_pathabun_PAPRICA_acc_metrics.rds")

mammal_pathabun_paprica_rho <- readRDS("mammal_pathabun_PAPRICA_spearman_df.rds")
mammal_pathabun_paprica_acc <- readRDS("mammal_pathabun_PAPRICA_acc_metrics.rds")

ocean_pathabun_paprica_rho <- readRDS("ocean_pathabun_PAPRICA_spearman_df.rds")
ocean_pathabun_paprica_acc <- readRDS("ocean_pathabun_PAPRICA_acc_metrics.rds")

blueberry_pathabun_paprica_rho <- readRDS("blueberry_pathabun_PAPRICA_spearman_df.rds")
blueberry_pathabun_paprica_acc <- readRDS("blueberry_pathabun_PAPRICA_acc_metrics.rds")


hmp_pathabun_paprica_rho$dataset <- "HMP"
hmp_pathabun_paprica_acc$dataset <- "HMP"

mammal_pathabun_paprica_rho$dataset <- "Mammal"
mammal_pathabun_paprica_acc$dataset <- "Mammal"

ocean_pathabun_paprica_rho$dataset <- "Ocean"
ocean_pathabun_paprica_acc$dataset <- "Ocean"

blueberry_pathabun_paprica_rho$dataset <- "Soil (Blueberry)"
blueberry_pathabun_paprica_acc$dataset <- "Soil (Blueberry)"

combined_pathabun_paprica_rho <- rbind(hmp_pathabun_paprica_rho, mammal_pathabun_paprica_rho,
                                       ocean_pathabun_paprica_rho, blueberry_pathabun_paprica_rho)

combined_pathabun_paprica_acc <- rbind(hmp_pathabun_paprica_acc, mammal_pathabun_paprica_acc,
                                       ocean_pathabun_paprica_acc, blueberry_pathabun_paprica_acc)

combined_pathabun_paprica_rho_melt <- melt(combined_pathabun_paprica_rho)
combined_pathabun_paprica_acc_melt <- melt(combined_pathabun_paprica_acc)

combined_pathabun_paprica_acc_melt_precision <- combined_pathabun_paprica_acc_melt[which(combined_pathabun_paprica_acc_melt$variable == "precision"), ]
combined_pathabun_paprica_acc_melt_recall <- combined_pathabun_paprica_acc_melt[which(combined_pathabun_paprica_acc_melt$variable == "recall"), ]

combined_pathabun_paprica_rho_plot <- ggplot(combined_pathabun_paprica_rho_melt, aes(x=cat, y=value, fill=cat)) + geom_boxplot() +
                                            ylim(c(0.5, 1)) +
                                            ylab(c("Spearman Correlation Coefficient")) +
                                            xlab("") +
                                            guides(fill=FALSE) +
                                            facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                  axis.text.x=element_text(angle=45, hjust=1)) +
                                            scale_fill_manual(values=c("light grey", "#F8766D"))

combined_pathabun_paprica_precision_plot <- ggplot(combined_pathabun_paprica_acc_melt_precision, aes(x=category, y=value, fill=category)) + geom_boxplot() +
                                                    ylim(c(0.5, 1)) +
                                                    ylab(c("Precision")) +
                                                    xlab("") +
                                                    guides(fill=FALSE) +
                                                    facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                          axis.text.x=element_text(angle=45, hjust=1)) +
                                                    scale_fill_manual(values=c("light grey", "#F8766D"))

combined_pathabun_paprica_recall_plot <- ggplot(combined_pathabun_paprica_acc_melt_recall, aes(x=category, y=value, fill=category)) + geom_boxplot() +
  ylim(c(0.5, 1)) +
  ylab(c("Recall")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#F8766D"))


plot_grid(combined_pathabun_paprica_rho_plot, combined_pathabun_paprica_precision_plot, combined_pathabun_paprica_recall_plot,
          labels = c("A", "B", "C"), ncol=2, nrow=2)
