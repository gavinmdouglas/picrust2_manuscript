### Code to make figure plotting pathabun spearman boxplots and also to plot correlations for IMG phenotypes.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories.

rm(list=ls(all.names=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

### PANEL A - pathabun spearman correlations.

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

pathabun_rho_outlist_raw <- list()
pathabun_rho <- list()
pathabun_rho_wilcoxon <- list()

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroonian", "indian"="Indian", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")


for(d in datasets) {
  
  pathabun_rho_outlist_raw[[d]] <- parse_rho_rds_and_calc_wilcoxon(rho_rds = paste(d, "_pathabun_spearman_df.rds", sep=""),
                                                             dataset_name = dataset2name[[d]],
                                                             wilcox_cat2ignore = extra_nsti_categories,
                                                             y_pos_start = 0.97,
                                                             dist_to_add=0.03)
  pathabun_rho[[d]] <- pathabun_rho_outlist_raw[[d]][[1]]
  pathabun_rho_wilcoxon[[d]] <- pathabun_rho_outlist_raw[[d]][[2]]
  
}

combined_pathabun_rho <- do.call("rbind", pathabun_rho)
combined_pathabun_rho_wilcoxon <- do.call("rbind", pathabun_rho_wilcoxon)

combined_pathabun_rho <- combined_pathabun_rho[-which(combined_pathabun_rho$cat == "Scrambled"), ]
combined_pathabun_rho_wilcoxon <- combined_pathabun_rho_wilcoxon[-which(combined_pathabun_rho_wilcoxon$group2 == "Scrambled"), ]

combined_pathabun_rho_no_nsti <- combined_pathabun_rho
combined_pathabun_rho_no_nsti$cat <- as.character(combined_pathabun_rho_no_nsti$cat)
combined_pathabun_rho_no_nsti <- combined_pathabun_rho_no_nsti[-which(combined_pathabun_rho_no_nsti$cat %in% extra_nsti_categories) ,]
combined_pathabun_rho_no_nsti[which(combined_pathabun_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2"
combined_pathabun_rho_no_nsti$cat <- factor(combined_pathabun_rho_no_nsti$cat,
                                            levels=c("Null", "PICRUSt2"))

combined_pathabun_rho_wilcoxon_no_nsti <- combined_pathabun_rho_wilcoxon
combined_pathabun_rho_wilcoxon_no_nsti[which(combined_pathabun_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2"

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

pathabun_rho_boxplots <- ggplot(combined_pathabun_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0.485, 1)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  guides(fill=FALSE) +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "gray95"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  stat_pvalue_manual(combined_pathabun_rho_wilcoxon_no_nsti, label = "clean_p") 


### PANEL B - IMG PHENOTYPES HOLDOUT VALIDATIONS.

combined_acc_by_phenotype <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/IMG_pheno_LOOCV_metrics.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, comment.char="")

combined_acc_by_phenotype_subset <- combined_acc_by_phenotype[,c("sample", "F1", "precision", "recall", "category")]
colnames(combined_acc_by_phenotype_subset) <- c("Phenotype", "F1 Score", "Precision", "Recall", "Category")
combined_acc_by_phenotype_subset <- combined_acc_by_phenotype_subset[with(combined_acc_by_phenotype_subset, order(Category, Phenotype)),]

combined_acc_by_phenotype_subset$Category <- factor(combined_acc_by_phenotype_subset$Category, levels=c("Null", "PICRUSt2"))
combined_acc_by_phenotype_subset_melt <- melt(combined_acc_by_phenotype_subset)

phenotype_F1_wilcox <- wilcox.test(combined_acc_by_phenotype_subset$`F1 Score`[which(combined_acc_by_phenotype_subset$Category=="Null")],
                                         combined_acc_by_phenotype_subset$`F1 Score`[which(combined_acc_by_phenotype_subset$Category=="PICRUSt2")], paired=TRUE)

phenotype_precision_wilcox <- wilcox.test(combined_acc_by_phenotype_subset$Precision[which(combined_acc_by_phenotype_subset$Category=="Null")],
                                          combined_acc_by_phenotype_subset$Precision[which(combined_acc_by_phenotype_subset$Category=="PICRUSt2")], paired=TRUE)

phenotype_recall_wilcox <- wilcox.test(combined_acc_by_phenotype_subset$Recall[which(combined_acc_by_phenotype_subset$Category=="Null")],
                                       combined_acc_by_phenotype_subset$Recall[which(combined_acc_by_phenotype_subset$Category=="PICRUSt2")], paired=TRUE)

phenotype_wilcox_p_df <- data.frame(variable=c("F1 Score", "Precision", "Recall"),
                                    group1="Null",
                                    group2="PICRUSt2",
                                    pval=c(phenotype_F1_wilcox$p.value,
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

combined_acc_by_phenotype_subset_melt$value[which(is.na(combined_acc_by_phenotype_subset_melt$value))] <- 0

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
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "gray95", colour = NA),
        panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1))

# Plot figure (11x5)

pdf(file = "../../../figures/Figure3.pdf", width=14.5, height=5)

plot_grid(pathabun_rho_boxplots, IMG_pheno_boxplots,
          labels = c("a", "b"), ncol=2, nrow=1, rel_widths = c(1, 0.5))

dev.off()
