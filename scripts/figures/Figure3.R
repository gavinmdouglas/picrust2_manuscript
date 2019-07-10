### Code to make figure plotting pathabun spearman boxplots and also to plot correlations for IMG phenotypes.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories.

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

### PANEL A - pathabun spearman correlations.

extra_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# HMP:
hmp_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "hmp_pathabun_spearman_df.rds",
                                                            dataset_name = "HMP",
                                                            wilcox_cat2ignore = extra_categories,
                                                            y_pos_start = 0.97)

hmp_pathabun_rho <- hmp_pathabun_rho_outlist[[1]]
hmp_pathabun_rho_wilcoxon <- hmp_pathabun_rho_outlist[[2]]


# Mammal:
mammal_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "mammal_pathabun_spearman_df.rds",
                                                               dataset_name = "Mammal",
                                                               wilcox_cat2ignore = extra_categories,
                                                               y_pos_start = 0.97)

mammal_pathabun_rho <- mammal_pathabun_rho_outlist[[1]]
mammal_pathabun_rho_wilcoxon <- mammal_pathabun_rho_outlist[[2]]


# Ocean:
ocean_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "ocean_pathabun_spearman_df.rds",
                                                              dataset_name = "Ocean",
                                                              wilcox_cat2ignore = extra_categories,
                                                              y_pos_start = 0.97)

ocean_pathabun_rho <- ocean_pathabun_rho_outlist[[1]]
ocean_pathabun_rho_wilcoxon <- ocean_pathabun_rho_outlist[[2]]


# Soil (Blueberry):
blueberry_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "blueberry_pathabun_spearman_df.rds",
                                                                  dataset_name = "Soil (Blueberry)",
                                                                  wilcox_cat2ignore = extra_categories,
                                                                  y_pos_start = 0.97)


blueberry_pathabun_rho <- blueberry_pathabun_rho_outlist[[1]]
blueberry_pathabun_rho_wilcoxon <- blueberry_pathabun_rho_outlist[[2]]



# Rho
combined_pathabun_rho <- rbind(hmp_pathabun_rho, mammal_pathabun_rho,
                               ocean_pathabun_rho, blueberry_pathabun_rho)

combined_pathabun_rho_wilcoxon <- rbind(hmp_pathabun_rho_wilcoxon, mammal_pathabun_rho_wilcoxon,
                                        ocean_pathabun_rho_wilcoxon, blueberry_pathabun_rho_wilcoxon)

combined_pathabun_rho_no_nsti <- combined_pathabun_rho
combined_pathabun_rho_no_nsti$cat <- as.character(combined_pathabun_rho_no_nsti$cat)
combined_pathabun_rho_no_nsti <- combined_pathabun_rho_no_nsti[-which(combined_pathabun_rho_no_nsti$cat %in% extra_categories) ,]
combined_pathabun_rho_no_nsti[which(combined_pathabun_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2 (ASVs)"
combined_pathabun_rho_no_nsti[which(combined_pathabun_rho_no_nsti$cat == "NSTI=2 (GG)"), "cat"] <- "PICRUSt2 (GG)"
combined_pathabun_rho_no_nsti$cat <- factor(combined_pathabun_rho_no_nsti$cat,
                                            levels=c("Null", "PICRUSt2 (GG)", "PICRUSt2 (ASVs)"))

combined_pathabun_rho_wilcoxon_no_nsti <- combined_pathabun_rho_wilcoxon
combined_pathabun_rho_wilcoxon_no_nsti[which(combined_pathabun_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2 (ASVs)"
combined_pathabun_rho_wilcoxon_no_nsti[which(combined_pathabun_rho_wilcoxon_no_nsti$group2 == "NSTI=2 (GG)"), "group2"] <- "PICRUSt2 (GG)"

combined_pathabun_rho_no_nsti_melt <- melt(combined_pathabun_rho_no_nsti)

combined_pathabun_rho_wilcoxon_no_nsti$p_symbol[which(combined_pathabun_rho_wilcoxon_no_nsti$p_symbol == "ns")] <- ""

combined_pathabun_rho_wilcoxon_no_nsti$to_plot <- combined_pathabun_rho_wilcoxon_no_nsti$raw_p

low_num <- which(combined_pathabun_rho_wilcoxon_no_nsti$to_plot < 0.001)

combined_pathabun_rho_wilcoxon_no_nsti$to_plot[-low_num] <- round(combined_pathabun_rho_wilcoxon_no_nsti$to_plot[-low_num], digits = 3)

combined_pathabun_rho_wilcoxon_no_nsti$to_plot[low_num] <- formatC(combined_pathabun_rho_wilcoxon_no_nsti$to_plot[low_num], format = "e", digits = 2)


combined_pathabun_rho_wilcoxon_no_nsti$clean_p <- paste("P=",
                                                        combined_pathabun_rho_wilcoxon_no_nsti$to_plot,
                                                        combined_pathabun_rho_wilcoxon_no_nsti$p_symbol,
                                                        sep="")

# Clean up a p-value by hand that does not need to be in scientific notation:
combined_ec_rho_wilcoxon_no_nsti$clean_p[which(combined_ec_rho_wilcoxon_no_nsti$clean_p == "P=7.81e-03*")] <- "P=0.00781*"

pathabun_rho_boxplots <- ggplot(combined_pathabun_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.6, 1.05)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  stat_pvalue_manual(combined_pathabun_rho_wilcoxon_no_nsti, label = "clean_p")


### PANEL B - IMG PHENOTYPES HOLDOUT VALIDATIONS.

combined_acc_by_phenotype <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/IMG_pheno_LOOCV_metrics.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, comment.char="")

combined_acc_by_phenotype_subset <- combined_acc_by_phenotype[,c("sample", "acc", "precision", "recall", "category")]
colnames(combined_acc_by_phenotype_subset) <- c("Phenotype", "Accuracy", "Precision", "Recall", "Category")
combined_acc_by_phenotype_subset <- combined_acc_by_phenotype_subset[with(combined_acc_by_phenotype_subset, order(Category, Phenotype)),]

combined_acc_by_phenotype_subset$Category <- factor(combined_acc_by_phenotype_subset$Category, levels=c("Null", "PICRUSt2"))
combined_acc_by_phenotype_subset_melt <- melt(combined_acc_by_phenotype_subset)

phenotype_accuracy_wilcox <- wilcox.test(combined_acc_by_phenotype_subset$Accuracy[which(combined_acc_by_phenotype_subset$Category=="Null")],
                                         combined_acc_by_phenotype_subset$Accuracy[which(combined_acc_by_phenotype_subset$Category=="PICRUSt2")], paired=TRUE)

phenotype_precision_wilcox <- wilcox.test(combined_acc_by_phenotype_subset$Precision[which(combined_acc_by_phenotype_subset$Category=="Null")],
                                          combined_acc_by_phenotype_subset$Precision[which(combined_acc_by_phenotype_subset$Category=="PICRUSt2")], paired=TRUE)

phenotype_recall_wilcox <- wilcox.test(combined_acc_by_phenotype_subset$Recall[which(combined_acc_by_phenotype_subset$Category=="Null")],
                                       combined_acc_by_phenotype_subset$Recall[which(combined_acc_by_phenotype_subset$Category=="PICRUSt2")], paired=TRUE)

phenotype_wilcox_p_df <- data.frame(variable=c("Accuracy", "Precision", "Recall"),
                                    group1="Null",
                                    group2="PICRUSt2",
                                    pval=c(phenotype_accuracy_wilcox$p.value,
                                           phenotype_precision_wilcox$p.value,
                                           phenotype_recall_wilcox$p.value),
                                    y.position=1.03,
                                    Category=NA)

phenotype_wilcox_p_df$p_symbol <- "ns"
phenotype_wilcox_p_df[which(phenotype_wilcox_p_df$pval < 0.05), "p_symbol"] <- "*"
phenotype_wilcox_p_df[which(phenotype_wilcox_p_df$pval < 0.001), "p_symbol"] <- "**"


phenotype_wilcox_p_df$clean_p <- paste("P=",
                                       formatC(phenotype_wilcox_p_df$pval, format = "e", digits = 0),
                                       phenotype_wilcox_p_df$p_symbol,
                                       sep="")

IMG_pheno_boxplots <- ggplot(combined_acc_by_phenotype_subset_melt, aes(x=Category, y=value, fill=Category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") +
  xlab("") +
  ylim(c(0, 1.05)) +
  ylab(c("Phenotype Prediction Performance")) +
  guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  stat_pvalue_manual(phenotype_wilcox_p_df, label = "clean_p", tip.length = 0.01) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1))

# Plot figure (11x5)
plot_grid(pathabun_rho_boxplots, IMG_pheno_boxplots, labels = c("A", "B"), ncol=2, nrow=1, rel_widths = c(1, 0.5))
