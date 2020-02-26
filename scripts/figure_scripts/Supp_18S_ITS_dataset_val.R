### Plot spearman correlations for Blueberry soil and wine fermentation EC/pathabun predictions vs MGS.
### Also plot how spearman correlations change as % eukaryotic varies.
### Also possibly plot correlation with example fungi-specific EC number.

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

# Read in metrics and prep per dataset.
# Blueberry:
blue_ec_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_ec_scc_metrics.rds",
                                                      dataset_name = "Soil (Blueberry)",
                                                      wilcox_cat2ignore = extra_nsti_categories,
                                                      y_pos_start = 0.75)

blue_ec_rho <- blue_ec_rho_outlist[[1]]
blue_ec_rho_wilcoxon <- blue_ec_rho_outlist[[2]]


# wine:
wine_ec_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "saved_RDS/18S_ITS_vs_MGS_metrics/wine_ec_scc_metrics.rds",
                                                       dataset_name = "Wine Fermentation",
                                                       wilcox_cat2ignore = extra_nsti_categories,
                                                       y_pos_start = 0.75)

wine_ec_rho <- wine_ec_rho_outlist[[1]]
wine_ec_rho_wilcoxon <- wine_ec_rho_outlist[[2]]

combined_ec_rho <- rbind(blue_ec_rho, wine_ec_rho)
combined_ec_rho_wilcoxon <- rbind(blue_ec_rho_wilcoxon,wine_ec_rho_wilcoxon)

combined_ec_rho_no_nsti <- combined_ec_rho
combined_ec_rho_no_nsti$cat <- as.character(combined_ec_rho_no_nsti$cat)
combined_ec_rho_no_nsti <- combined_ec_rho_no_nsti[-which(combined_ec_rho_no_nsti$cat %in% extra_nsti_categories) ,]
combined_ec_rho_no_nsti[which(combined_ec_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2"
combined_ec_rho_no_nsti$cat <- factor(combined_ec_rho_no_nsti$cat,
                                      levels=c("Null", "PICRUSt2"))

combined_ec_rho_wilcoxon_no_nsti <- combined_ec_rho_wilcoxon
combined_ec_rho_wilcoxon_no_nsti[which(combined_ec_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2"

# Add column of p-values to add to plot.
combined_ec_rho_wilcoxon_no_nsti$clean_p <- paste("P=",
                                                  formatC(combined_ec_rho_wilcoxon_no_nsti$raw_p, format = "e", digits = 2),
                                                  combined_ec_rho_wilcoxon_no_nsti$p_symbol,
                                                  sep="")
# Clean up a p-value by hand that does not need to be in scientific notation:
combined_ec_rho_wilcoxon_no_nsti$clean_p[which(combined_ec_rho_wilcoxon_no_nsti$clean_p == "P=7.81e-03*")] <- "P=0.00781*"

combined_ec_rho_no_nsti_melt <- melt(combined_ec_rho_no_nsti)

ec_scc_boxplots <- ggplot(combined_ec_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1) +
  ylim(c(0, 1)) +
  ylab(c("Spearman Correlation")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free") +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #      panel.background = element_blank(), axis.line = element_line(colour = "black"),
  #      axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  stat_pvalue_manual(combined_ec_rho_wilcoxon_no_nsti, label = "clean_p") +
  ggtitle("EC Numbers")

# Get mean and sd values.
mean(blue_ec_rho[which(blue_ec_rho$cat == "Null"), "metric"])
sd(blue_ec_rho[which(blue_ec_rho$cat == "Null"), "metric"])

mean(blue_ec_rho[which(blue_ec_rho$cat == "NSTI=2"), "metric"])
sd(blue_ec_rho[which(blue_ec_rho$cat == "NSTI=2"), "metric"])


mean(wine_ec_rho[which(wine_ec_rho$cat == "Null"), "metric"])
sd(wine_ec_rho[which(wine_ec_rho$cat == "Null"), "metric"])

mean(wine_ec_rho[which(wine_ec_rho$cat == "NSTI=2"), "metric"])
sd(wine_ec_rho[which(wine_ec_rho$cat == "NSTI=2"), "metric"])


# Do the same thing, but for pathway abundances.
blue_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_pathabun_scc_metrics.rds",
                                                       dataset_name = "Soil (Blueberry)",
                                                       wilcox_cat2ignore = extra_nsti_categories,
                                                       y_pos_start = 0.75)

blue_pathabun_rho <- blue_pathabun_rho_outlist[[1]]
blue_pathabun_rho_wilcoxon <- blue_pathabun_rho_outlist[[2]]

blue_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_pathabun_scc_metrics.rds",
                                                       dataset_name = "Soil (Blueberry)",
                                                       wilcox_cat2ignore = extra_nsti_categories,
                                                       y_pos_start = 0.75)

blue_pathabun_rho <- blue_pathabun_rho_outlist[[1]]
blue_pathabun_rho_wilcoxon <- blue_pathabun_rho_outlist[[2]]


# wine:
wine_pathabun_rho_outlist <- parse_rho_rds_and_calc_wilcoxon(rho_rds = "saved_RDS/18S_ITS_vs_MGS_metrics/wine_pathabun_scc_metrics.rds",
                                                       dataset_name = "Wine Fermentation",
                                                       wilcox_cat2ignore = extra_nsti_categories,
                                                       y_pos_start = 0.75)

wine_pathabun_rho <- wine_pathabun_rho_outlist[[1]]
wine_pathabun_rho_wilcoxon <- wine_pathabun_rho_outlist[[2]]





combined_pathabun_rho <- rbind(blue_pathabun_rho, wine_pathabun_rho)
combined_pathabun_rho_wilcoxon <- rbind(blue_pathabun_rho_wilcoxon,wine_pathabun_rho_wilcoxon)

combined_pathabun_rho_no_nsti <- combined_pathabun_rho
combined_pathabun_rho_no_nsti$cat <- as.character(combined_pathabun_rho_no_nsti$cat)
combined_pathabun_rho_no_nsti <- combined_pathabun_rho_no_nsti[-which(combined_pathabun_rho_no_nsti$cat %in% extra_nsti_categories) ,]
combined_pathabun_rho_no_nsti[which(combined_pathabun_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2"
combined_pathabun_rho_no_nsti$cat <- factor(combined_pathabun_rho_no_nsti$cat,
                                      levels=c("Null", "PICRUSt2"))

combined_pathabun_rho_wilcoxon_no_nsti <- combined_pathabun_rho_wilcoxon
combined_pathabun_rho_wilcoxon_no_nsti[which(combined_pathabun_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2"

combined_pathabun_rho_no_nsti_melt <- melt(combined_pathabun_rho_no_nsti)

combined_pathabun_rho_wilcoxon_no_nsti$clean_p <- paste("P=",
                                                  formatC(combined_pathabun_rho_wilcoxon_no_nsti$raw_p, format = "e", digits = 2),
                                                  combined_pathabun_rho_wilcoxon_no_nsti$p_symbol,
                                                  sep="")

# Clean up a p-value by hand that does not need to be in scientific notation:
combined_pathabun_rho_wilcoxon_no_nsti$clean_p[which(combined_pathabun_rho_wilcoxon_no_nsti$clean_p == "P=7.81e-03*")] <- "P=0.00781*"
combined_pathabun_rho_wilcoxon_no_nsti$clean_p[which(combined_pathabun_rho_wilcoxon_no_nsti$clean_p == "P=5.47e-02ns")] <- "P=0.0547"

pathabun_scc_boxplots <- ggplot(combined_pathabun_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1) +
  ylim(c(0, 1)) +
  ylab(c("Spearman Correlation")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free") +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #      panel.background = element_blank(), axis.line = element_line(colour = "black"),
  #      axis.text.x=element_text(angle=45, hjust=1)) +
  guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  stat_pvalue_manual(combined_pathabun_rho_wilcoxon_no_nsti, label = "clean_p") +
  ggtitle("MetaCyc Pathways")


# Get mean and sd values.
mean(blue_pathabun_rho[which(blue_pathabun_rho$cat == "Null"), "metric"])
sd(blue_pathabun_rho[which(blue_pathabun_rho$cat == "Null"), "metric"])

mean(blue_pathabun_rho[which(blue_pathabun_rho$cat == "NSTI=2"), "metric"])
sd(blue_pathabun_rho[which(blue_pathabun_rho$cat == "NSTI=2"), "metric"])


mean(wine_pathabun_rho[which(wine_pathabun_rho$cat == "Null"), "metric"])
sd(wine_pathabun_rho[which(wine_pathabun_rho$cat == "Null"), "metric"])

mean(wine_pathabun_rho[which(wine_pathabun_rho$cat == "NSTI=2"), "metric"])
sd(wine_pathabun_rho[which(wine_pathabun_rho$cat == "NSTI=2"), "metric"])

### PANEL C
# Read in % of MGS reads mapped to each phyla (based on metaxa2).
# Get percent eukaryotic after removing Metazoa and Virdiplantae.
blue_metaxa2 <- read.table("working_tables/metaxa2_counts/blueberry_metaxa2_phyla_counts.tsv",
                                header=T, row.names=1, sep="\t", check.names = FALSE)
rownames(blue_metaxa2) <- gsub("BB", "", rownames(blue_metaxa2))
rownames(blue_metaxa2) <- gsub("_", "-", rownames(blue_metaxa2))
blue_metaxa2_relab <- data.frame(sweep(blue_metaxa2, 1, rowSums(blue_metaxa2), '/'), check.names = FALSE) * 100
blue_metaxa2_relab_euk <- blue_metaxa2_relab[, grep("^Eukaryota\\|", colnames(blue_metaxa2_relab))]
blue_metaxa2_relab_euk_subset <- blue_metaxa2_relab_euk[, -which(colnames(blue_metaxa2_relab_euk) %in% c("Eukaryota|Metazoa", "Eukaryota|Viridiplantae"))]
blue_metaxa2_relab_euk_subset$total_filt <- rowSums(blue_metaxa2_relab_euk_subset)

blue_ec_rho_nsti2 <- blue_ec_rho[which(blue_ec_rho$cat == "NSTI=2"),]
rownames(blue_ec_rho_nsti2) <- as.character(blue_ec_rho_nsti2$sample_names)

wine_metaxa2 <- read.table("working_tables/metaxa2_counts/wine_metaxa2_phyla_counts.tsv",
                                header=T, row.names=1, sep="\t", check.names = FALSE)
wine_id_map <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/wine_id_mapping.txt",
                          header=T, sep="\t", stringsAsFactors = FALSE)
rownames(wine_id_map) <- wine_id_map$MGS_run
rownames(wine_metaxa2) <- wine_id_map[rownames(wine_metaxa2), "X16S_name"]

wine_metaxa2_relab <- data.frame(sweep(wine_metaxa2, 1, rowSums(wine_metaxa2), '/'), check.names = FALSE) * 100
wine_metaxa2_relab_euk <- wine_metaxa2_relab[, grep("^Eukaryota\\|", colnames(wine_metaxa2_relab))]
wine_metaxa2_relab_euk_subset <- wine_metaxa2_relab_euk[, -which(colnames(wine_metaxa2_relab_euk) %in% c("Eukaryota|Metazoa", "Eukaryota|Viridiplantae"))]
wine_metaxa2_relab_euk_subset$total_filt <- rowSums(wine_metaxa2_relab_euk_subset)

wine_ec_rho_nsti2 <- wine_ec_rho[which(wine_ec_rho$cat == "NSTI=2"),]
rownames(wine_ec_rho_nsti2) <- as.character(wine_ec_rho_nsti2$sample_names)

blue_metaxa2_relab_euk_subset$rho <- blue_ec_rho_nsti2[rownames(blue_metaxa2_relab_euk), "metric"]
wine_metaxa2_relab_euk_subset$rho <- wine_ec_rho_nsti2[rownames(wine_metaxa2_relab_euk), "metric"]

blue_metaxa2_relab_euk_subset$Dataset <- "Soil (Blueberry)"
wine_metaxa2_relab_euk_subset$Dataset <- "Wine Fermentation"

blue_metaxa2_relab_euk_subset <- blue_metaxa2_relab_euk_subset[, c("Dataset", "total_filt", "rho")]
wine_metaxa2_relab_euk_subset <- wine_metaxa2_relab_euk_subset[, c("Dataset", "total_filt", "rho")]

combined_metaxa2_rho <- rbind(blue_metaxa2_relab_euk_subset, wine_metaxa2_relab_euk_subset)

eukaryote_rho_scatterplot <- ggplot(combined_metaxa2_rho, aes(x=total_filt, y=rho, colour=Dataset)) +
  geom_point() +
  ylim(c(0, 1)) +
  ylab(c("Spearman Correlation")) +
  xlab("% Eukaryotic DNA") +
  scale_colour_manual(values=c("deepskyblue3", "firebrick3")) +
  theme(legend.position = c(0.05, 0.8), legend.background = element_rect(color = "black", 
                                                                        fill = "white", size = 0.3, linetype = "solid"))


# Get mean and sd % eukaryotic for each dataset
mean(blue_metaxa2_relab_euk_subset$total_filt)
sd(blue_metaxa2_relab_euk_subset$total_filt)

mean(wine_metaxa2_relab_euk_subset$total_filt)
sd(wine_metaxa2_relab_euk_subset$total_filt)


# PANEL D: Sig ECs in blueberry dataset.
blueberry_sig_ec <- read.table("working_tables/all_blueberry_sig_ecs.tsv",
                                 header=TRUE, sep="\t", check.names=FALSE)

# Restrict to EC numbers with a mean of at least 0.15% relative abundance across all samples.

EC_col <- grep("EC:", colnames(blueberry_sig_ec), value = TRUE)

common_EC <- names(which(colMeans(blueberry_sig_ec[, EC_col]) > 0.15))

blueberry_sig_ec <- blueberry_sig_ec[, c(common_EC, "group", "sample")]

blueberry_sig_ec_melt <- melt(blueberry_sig_ec)

blueberry_sig_ec_melt$group <- as.character(blueberry_sig_ec_melt$group)
blueberry_sig_ec_melt$Type <- factor(gsub("Mng", "", blueberry_sig_ec_melt$group), levels=c("Bulk", "Rhizo", "Root"))

blueberry_sig_ecs_boxplots <- ggplot(blueberry_sig_ec_melt, aes(x=Type, y=value, fill=Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  ylim(c(0, 1.5)) +
  ylab("% Relative Abundance") +
  xlab("") +
  scale_fill_manual("Sample Type", values=c("light grey", "cornflowerblue", "orange")) +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1))


### Plot final figure.
top_row <- plot_grid(ec_scc_boxplots,
                     pathabun_scc_boxplots,
                     eukaryote_rho_scatterplot,
                     labels = c('a', 'b', 'c'),
                     ncol=3,
                     nrow=1)
                    
pdf(file = "../figures/Supp_eukaryotic_validation.pdf", width=13, height=7)

plot_grid(top_row,
          blueberry_sig_ecs_boxplots,
          labels = c('', 'd'),
          ncol = 1,
          nrow=2)

dev.off()
