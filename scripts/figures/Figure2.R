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

dataset2name <- list("cameroon"="Cameroon", "indian"="India", "hmp"="HMP", "mammal"="Mammal",
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
combined_ko_rho_no_nsti <- combined_ko_rho_no_nsti[-which(combined_ko_rho_no_nsti$cat %in% extra_nsti_categories) ,]
combined_ko_rho_no_nsti[which(combined_ko_rho_no_nsti$cat == "NSTI=2"), "cat"] <- "PICRUSt2"
combined_ko_rho_no_nsti$cat <- factor(combined_ko_rho_no_nsti$cat,
                                      levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))

combined_ko_rho_wilcoxon_no_nsti <- combined_ko_rho_wilcoxon
combined_ko_rho_wilcoxon_no_nsti[which(combined_ko_rho_wilcoxon_no_nsti$group1 == "NSTI=2"), "group1"] <- "PICRUSt2"

combined_ko_rho_no_nsti_melt <- melt(combined_ko_rho_no_nsti)

combined_ko_rho_no_nsti_melt$dataset <- factor(combined_ko_rho_no_nsti_melt$dataset, levels=c("Cameroon", "Indian", "HMP", "Primate", "Mammal", "Ocean", "Soil (Blueberry)"))

pdf(file = "../../../figures/Figure2.pdf", width=12, height=6, units = "in")

ggplot(combined_ko_rho_no_nsti_melt, aes(x=cat, y=value, fill=Database)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=0.1) +
  scale_y_continuous(breaks=c(0.6, 0.8, 1.0), limits=c(0.485, 1.06)) +
  ylab(c("Spearman Correlation Coefficient")) +
  xlab("") +
  facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position = c(0.05, 0.12), legend.background = element_rect(color = "black", 
                                                                          fill = "white", size = 0.3, linetype = "solid"),
        legend.title = element_text(colour="black", size=8, face="bold"),
        legend.text = element_text(colour="black", size=8)) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  stat_pvalue_manual(data = combined_ko_rho_wilcoxon_no_nsti, label = "p_symbol", bracket.size = 0.2, tip.length = 0.01, label.size = 3)

dev.off()
