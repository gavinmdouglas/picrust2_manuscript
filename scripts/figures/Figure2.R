### Code to make figure contrasting KO correlations on each 16S validation dataset.
### These metrics were calculated by assuming the metagenomics sequencing was the "gold standard".
### Also tests for statistical significance between these categories.

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("../../../scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metics per dataset

ko_rho_outlist_raw <- list()
ko_rho <- list()
ko_rho_wilcoxon <- list()

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroonian", "indian"="Indian", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")


for(d in datasets) {
 
  ko_rho_outlist_raw[[d]] <- parse_rho_rds_and_calc_wilcoxon(rho_rds = paste(d, "_ko_spearman_df.rds", sep=""),
                                                             dataset_name = dataset2name[[d]],
                                                             wilcox_cat2ignore = extra_nsti_categories,
                                                             y_pos_start = 0.94,
                                                             dist_to_add=0.03)
  ko_rho[[d]] <- ko_rho_outlist_raw[[d]][[1]]
  ko_rho_wilcoxon[[d]] <- ko_rho_outlist_raw[[d]][[2]]
   
}

combined_ko_rho <- do.call("rbind", ko_rho)
combined_ko_rho_wilcoxon <- do.call("rbind", ko_rho_wilcoxon)

combined_ko_rho_no_nsti <- combined_ko_rho
combined_ko_rho_no_nsti$cat <- as.character(combined_ko_rho_no_nsti$cat)
combined_ko_rho_no_nsti <- combined_ko_rho_no_nsti[-which(combined_ko_rho_no_nsti$cat %in% extra_nsti_categories), ]
combined_ko_rho_no_nsti[which(combined_ko_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2"

combined_ko_rho_no_nsti <- combined_ko_rho_no_nsti[-which(combined_ko_rho_no_nsti$cat == "Scrambled"), ]

combined_ko_rho_no_nsti$cat <- factor(combined_ko_rho_no_nsti$cat,
                                      levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))

combined_ko_rho_wilcoxon_no_nsti <- combined_ko_rho_wilcoxon
combined_ko_rho_wilcoxon_no_nsti[which(combined_ko_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2"
combined_ko_rho_wilcoxon_no_nsti <- combined_ko_rho_wilcoxon_no_nsti[-which(combined_ko_rho_wilcoxon_no_nsti$group2 == "Scrambled"), ]
combined_ko_rho_wilcoxon_no_nsti$group2

combined_ko_rho_no_nsti_melt <- melt(combined_ko_rho_no_nsti)

combined_ko_rho_no_nsti_melt$dataset <- factor(combined_ko_rho_no_nsti_melt$dataset, levels=c("Cameroonian", "Indian", "HMP", "Primate", "Mammal", "Ocean", "Soil (Blueberry)"))

combined_ko_rho_no_nsti_melt[which(combined_ko_rho_no_nsti_melt$cat == "Scrambled\nASVs"), "Database"] <- "PICRUSt2"

ko_rho_boxplots <- ggplot(combined_ko_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
                        geom_boxplot(outlier.shape = NA) +
                        geom_quasirandom(size=0.1) +
                        scale_y_continuous(breaks=c(0.6, 0.8, 1.0), limits=c(0.485, 1.10)) +
                        ylab(c("Spearman Correlation Coefficient")) +
                        xlab("") +
                        facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black"),
                              axis.text.x=element_text(angle=45, hjust=1),
                              legend.position = c(0.07, 0.2), legend.background = element_rect(color = "black", 
                                                                                                fill = "white", size = 0.3, linetype = "solid"),
                              legend.title = element_text(colour="black", size=8, face="bold"),
                              legend.text = element_text(colour="black", size=8)) +
                        scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
                        stat_pvalue_manual(data = combined_ko_rho_wilcoxon_no_nsti, label = "p_symbol", bracket.size = 0.2, tip.length = 0.01, label.size = 2.5)



wilcoxon_out_musicc_perf_0.05 <- readRDS(file="../DA_concordance/ko_wilcoxon_out_musicc_perf_0.05.rds")

wilcoxon_out_musicc_perf_0.05$hmp$metric <- rownames(wilcoxon_out_musicc_perf_0.05$hmp)
wilcoxon_out_musicc_perf_0.05_hmp <- melt(wilcoxon_out_musicc_perf_0.05$hmp)
wilcoxon_out_musicc_perf_0.05_hmp$dataset <- "HMP\n22 supragingival plaque vs 36 tongue dorsum"

wilcoxon_out_musicc_perf_0.05$cameroon$metric <- rownames(wilcoxon_out_musicc_perf_0.05$cameroon)
wilcoxon_out_musicc_perf_0.05_cameroon <- melt(wilcoxon_out_musicc_perf_0.05$cameroon)
wilcoxon_out_musicc_perf_0.05_cameroon$dataset <- "Cameroonian\n19 parasite+ vs 36 parasite- individuals"

wilcoxon_out_musicc_perf_0.05$indian$metric <- rownames(wilcoxon_out_musicc_perf_0.05$indian)
wilcoxon_out_musicc_perf_0.05_indian <- melt(wilcoxon_out_musicc_perf_0.05$indian)
wilcoxon_out_musicc_perf_0.05_indian$dataset <- "Indian\n51 individuals from Bhopal vs 38 from Kerala"

wilcoxon_out_musicc_perf_0.05$primate$metric <- rownames(wilcoxon_out_musicc_perf_0.05$primate)
wilcoxon_out_musicc_perf_0.05_primate <- melt(wilcoxon_out_musicc_perf_0.05$primate)
wilcoxon_out_musicc_perf_0.05_primate$dataset <- "Primate\n29 old world vs 29 new world monkeys"

wilcoxon_out_musicc_perf_0.05_combined <- rbind(wilcoxon_out_musicc_perf_0.05_hmp, wilcoxon_out_musicc_perf_0.05_cameroon,
                                                wilcoxon_out_musicc_perf_0.05_indian, wilcoxon_out_musicc_perf_0.05_primate)

wilcoxon_out_musicc_perf_0.05_combined$variable <- as.character(wilcoxon_out_musicc_perf_0.05_combined$variable)
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "mgs_ko_alt")] <- "Alt. MGS"
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "panfp_ko")] <- "PanFP"
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "picrust1_ko")] <- "PICRUSt1"
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "picrust2_ko_nsti2")] <- "PICRUSt2"
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "picrust2_scrambled")] <- "Scrambled\nASVs"
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "piphillin_ko")] <- "Piphillin"
wilcoxon_out_musicc_perf_0.05_combined$variable[which(wilcoxon_out_musicc_perf_0.05_combined$variable == "tax4fun2_ko")] <- "Tax4Fun2"

wilcoxon_out_musicc_perf_0.05_combined$variable <- factor(wilcoxon_out_musicc_perf_0.05_combined$variable, levels=c("Scrambled\nASVs", "Alt. MGS", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))

# wilcoxon_out_musicc_perf_0.05_combined <- wilcoxon_out_musicc_perf_0.05_combined[-which(wilcoxon_out_musicc_perf_0.05_combined$variable %in% c("PanFP", "Tax4Fun2", "PICRUSt1")), ]
# 
# wilcoxon_out_musicc_perf_0.05_combined$variable <- factor(wilcoxon_out_musicc_perf_0.05_combined$variable,
#                                                           levels=c("Scrambled\nASVs", "Alt. MGS", "Piphillin", "PICRUSt2"))

wilcoxon_out_musicc_perf_0.05_combined_subset <- wilcoxon_out_musicc_perf_0.05_combined[which(wilcoxon_out_musicc_perf_0.05_combined$metric %in% c("precision", "recall", "f1")), ]
wilcoxon_out_musicc_perf_0.05_combined_subset[which(wilcoxon_out_musicc_perf_0.05_combined_subset$metric == "precision"), "metric"] <- "Precision"
wilcoxon_out_musicc_perf_0.05_combined_subset[which(wilcoxon_out_musicc_perf_0.05_combined_subset$metric == "recall"), "metric"] <- "Recall"
wilcoxon_out_musicc_perf_0.05_combined_subset[which(wilcoxon_out_musicc_perf_0.05_combined_subset$metric == "f1"), "metric"] <- "F1 Score"

wilcoxon_out_musicc_perf_0.05_combined_subset$metric <- factor(wilcoxon_out_musicc_perf_0.05_combined_subset$metric, levels=c("Precision", "Recall", "F1 Score"))

DA_acc_barplot <- ggplot(wilcoxon_out_musicc_perf_0.05_combined_subset,
                             aes(x=variable, y=value, fill=metric)) +
                        geom_bar(stat="identity") +
                        scale_fill_manual(values=c("#0072B2", "#E69F00", "#009E73")) +
                        facet_grid(metric ~ dataset, scales = "free", space = "free", switch="y") +
                        theme(panel.grid.minor = element_blank(),
                              panel.grid.major.x = element_blank(),
                              axis.line = element_line(colour = "black"),
                              axis.text.x=element_text(angle=45, hjust=1),
                              legend.position = "none",
                              panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(0.75,"lines")) +
                        ylab("") +
                        xlab("") +
                        scale_y_continuous(expand = c(0, 0), limits = c(0, 1))


pdf(file = "../../../figures/Figure2.pdf", width=12, height=8)

plot_grid(ko_rho_boxplots, DA_acc_barplot, nrow=2, labels=c('a', 'b'), rel_heights = c(1, 0.85))

dev.off()


# Calculate basic summary statistics.

dataset_ko_rho_stats <- data.frame(matrix(NA, nrow=7, ncol=2))
colnames(dataset_ko_rho_stats) <- c("PICRUSt2_mean", "PICRUSt2_sd")
rownames(dataset_ko_rho_stats) <- datasets

for(d in datasets) {
  
  d_name <- dataset2name[[d]]
  
  dataset_ko_rho_stats[d, ] <- c(mean(combined_ko_rho[which(combined_ko_rho$dataset == d_name & combined_ko_rho$cat == "NSTI=2"), "metric"]),
                                 sd(combined_ko_rho[which(combined_ko_rho$dataset == d_name & combined_ko_rho$cat == "NSTI=2"), "metric"]))
  
}

# HMP dataset PICRUSt2 vs Piphillin is difficult to interpret based on means / medians since the distributions are so similar.
# Looking at the distribution of per-sample differences is the best way to interpret how they differ - which shows
# that PICRUSt2 does *slightly* better.
summary(combined_ko_rho[which(combined_ko_rho$dataset == "HMP" & combined_ko_rho$cat == "NSTI=2"), "metric"] -
        combined_ko_rho[which(combined_ko_rho$dataset == "HMP" & combined_ko_rho$cat == "Piphillin"), "metric"])
