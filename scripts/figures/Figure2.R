### Code to make figure contrasting KO correlations on each 16S validation dataset.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories.

rm(list=ls())

library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# HMP:
hmp_ko_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "hmp_ko_spearman_df.rds",
                                                      dataset_name = "HMP",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)


hmp_ko_rho <- hmp_ko_rho_outlist[[1]]
hmp_ko_rho_wilcoxon <- hmp_ko_rho_outlist[[2]]


# Mammal:
mammal_ko_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "mammal_ko_spearman_df.rds",
                                                      dataset_name = "Mammal",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)

mammal_ko_rho <- mammal_ko_rho_outlist[[1]]
mammal_ko_rho_wilcoxon <- mammal_ko_rho_outlist[[2]]


# Ocean:
ocean_ko_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "ocean_ko_spearman_df.rds",
                                                      dataset_name = "Ocean",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)

ocean_ko_rho <- ocean_ko_rho_outlist[[1]]
ocean_ko_rho_wilcoxon <- ocean_ko_rho_outlist[[2]]


# Soil (Blueberry):
blueberry_ko_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "blueberry_ko_spearman_df.rds",
                                                      dataset_name = "Soil (Blueberry)",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.97)


blueberry_ko_rho <- blueberry_ko_rho_outlist[[1]]
blueberry_ko_rho_wilcoxon <- blueberry_ko_rho_outlist[[2]]



# Rho
combined_ko_rho <- rbind(hmp_ko_rho, mammal_ko_rho,
                         ocean_ko_rho, blueberry_ko_rho)

combined_ko_rho_wilcoxon <- rbind(hmp_ko_rho_wilcoxon, mammal_ko_rho_wilcoxon,
                                  ocean_ko_rho_wilcoxon, blueberry_ko_rho_wilcoxon)

combined_ko_rho_no_nsti <- combined_ko_rho
combined_ko_rho_no_nsti$cat <- as.character(combined_ko_rho_no_nsti$cat)
combined_ko_rho_no_nsti <- combined_ko_rho_no_nsti[-which(combined_ko_rho_no_nsti$cat %in% extra_nsti_categories) ,]
combined_ko_rho_no_nsti[which(combined_ko_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2 (ASVs)"
combined_ko_rho_no_nsti[which(combined_ko_rho_no_nsti$cat == "NSTI=2 (GG)"), "cat"] <- "PICRUSt2 (GG)"
combined_ko_rho_no_nsti$cat <- factor(combined_ko_rho_no_nsti$cat,
                                      levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2 (GG)", "PICRUSt2 (ASVs)"))

combined_ko_rho_wilcoxon_no_nsti <- combined_ko_rho_wilcoxon
combined_ko_rho_wilcoxon_no_nsti[which(combined_ko_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2 (ASVs)"
combined_ko_rho_wilcoxon_no_nsti[which(combined_ko_rho_wilcoxon_no_nsti$group2 == "NSTI=2 (GG)"), "group2"] <- "PICRUSt2 (GG)"

combined_ko_rho_no_nsti_melt <- melt(combined_ko_rho_no_nsti)

combined_ko_rho_no_nsti_melt$dataset

# Saved as 12 x 6.
ggplot(combined_ko_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.5, 1.31)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = c(0.05, 0.85), legend.background = element_rect(color = "black", 
                                                                          fill = "white", size = 0.3, linetype = "solid"),
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(combined_ko_rho_wilcoxon_no_nsti, label = "p_symbol")


# Get mean and sd for each dataset.
mean(hmp_ko_rho[which(hmp_ko_rho$cat=="NSTI=2"), "metric"])
sd(hmp_ko_rho[which(hmp_ko_rho$cat=="NSTI=2"), "metric"])
# 0.8591101 0.0356301

mean(mammal_ko_rho[which(mammal_ko_rho$cat=="NSTI=2"), "metric"])
sd(mammal_ko_rho[which(mammal_ko_rho$cat=="NSTI=2"), "metric"])
# 0.8299978 0.01844903

mean(ocean_ko_rho[which(ocean_ko_rho$cat=="NSTI=2"), "metric"])
sd(ocean_ko_rho[which(ocean_ko_rho$cat=="NSTI=2"), "metric"])
# 0.8735957 0.006982587

mean(blueberry_ko_rho[which(blueberry_ko_rho$cat=="NSTI=2"), "metric"])
sd(blueberry_ko_rho[which(blueberry_ko_rho$cat=="NSTI=2"), "metric"])
# 0.8545425 0.01729161

