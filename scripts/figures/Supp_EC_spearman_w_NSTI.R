### Code to make figure contrasting EC correlations on each 16S validation dataset.
### Include all NSTI cut-offs in these plots.

rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# HMP:
hmp_ec_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "hmp_ec_spearman_df.rds",
                                                      dataset_name = "HMP",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)


hmp_ec_rho <- hmp_ec_rho_outlist[[1]]
hmp_ec_rho_wilcoxon <- hmp_ec_rho_outlist[[2]]


# Mammal:
mammal_ec_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "mammal_ec_spearman_df.rds",
                                                      dataset_name = "Mammal",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)

mammal_ec_rho <- mammal_ec_rho_outlist[[1]]
mammal_ec_rho_wilcoxon <- mammal_ec_rho_outlist[[2]]


# Ocean:
ocean_ec_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "ocean_ec_spearman_df.rds",
                                                      dataset_name = "Ocean",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)

ocean_ec_rho <- ocean_ec_rho_outlist[[1]]
ocean_ec_rho_wilcoxon <- ocean_ec_rho_outlist[[2]]


# Soil (Blueberry):
blueberry_ec_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "blueberry_ec_spearman_df.rds",
                                                      dataset_name = "Soil (Blueberry)",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)


blueberry_ec_rho <- blueberry_ec_rho_outlist[[1]]
blueberry_ec_rho_wilcoxon <- blueberry_ec_rho_outlist[[2]]



# Make plot for each dataset.

hmp_ec_rho$cat <- as.character(hmp_ec_rho$cat)
hmp_ec_rho$cat <- factor(hmp_ec_rho$cat,
                         levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                  "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

hmp_ec_rho_melt <- melt(hmp_ec_rho)

hmp_ec_spearman_boxplots <- ggplot(hmp_ec_rho_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.4, 1.31)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = c(0.6, 0.85), legend.background = element_rect(color = "black", 
                                                                          fill = "white", size = 0.3, linetype = "solid"),
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(hmp_ec_rho_wilcoxon, label = "p_symbol")


mammal_ec_rho$cat <- as.character(mammal_ec_rho$cat)
mammal_ec_rho$cat <- factor(mammal_ec_rho$cat,
                            levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                     "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

mammal_ec_rho_melt <- melt(mammal_ec_rho)

mammal_ec_spearman_boxplots <- ggplot(mammal_ec_rho_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.4, 1.31)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = c(0.6, 0.85), legend.background = element_rect(color = "black", 
                                                                         fill = "white", size = 0.3, linetype = "solid"),
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(mammal_ec_rho_wilcoxon, label = "p_symbol")


ocean_ec_rho$cat <- as.character(ocean_ec_rho$cat)
ocean_ec_rho$cat <- factor(ocean_ec_rho$cat,
                           levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                    "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

ocean_ec_rho_melt <- melt(ocean_ec_rho)

ocean_ec_spearman_boxplots <- ggplot(ocean_ec_rho_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.4, 1.31)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = c(0.6, 0.85), legend.background = element_rect(color = "black", 
                                                                         fill = "white", size = 0.3, linetype = "solid"),
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(ocean_ec_rho_wilcoxon, label = "p_symbol")

blueberry_ec_rho$cat <- as.character(blueberry_ec_rho$cat)
blueberry_ec_rho$cat <- factor(blueberry_ec_rho$cat,
                               levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5",
                                        "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

blueberry_ec_rho_melt <- melt(blueberry_ec_rho)

blueberry_ec_spearman_boxplots <- ggplot(blueberry_ec_rho_melt, aes(x=cat, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.4, 1.31)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = c(0.6, 0.85), legend.background = element_rect(color = "black", 
                                                                         fill = "white", size = 0.3, linetype = "solid"),
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(blueberry_ec_rho_wilcoxon, label = "p_symbol")

plot_grid(hmp_ec_spearman_boxplots,
          mammal_ec_spearman_boxplots,
          ocean_ec_spearman_boxplots,
          blueberry_ec_spearman_boxplots,
          labels=c("A", "B", "C", "D"),
          nrow=2,
          ncol=2)
