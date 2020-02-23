### Code to make figure contrasting KO correlations on each 16S validation dataset.
### Include all NSTI cut-offs in these plots.

rm(list=ls(all.names=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

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
                                                             dist_to_add=0.05)
  ko_rho[[d]] <- ko_rho_outlist_raw[[d]][[1]]
  ko_rho_wilcoxon[[d]] <- ko_rho_outlist_raw[[d]][[2]]
  ko_rho_wilcoxon[[d]][which(ko_rho_wilcoxon[[d]]$group2 == "Scrambled"), "group2"] <- "Shuffled\nASVs"
  
}


# Make plot for each dataset.
KO_spearman_boxplots <- list()

for(j in 1:length(datasets)) {
  
  d <- datasets[j]
  
  dataset_ko_rho <- ko_rho[[d]]
  
  dataset_ko_rho$cat <- as.character(dataset_ko_rho$cat)
  dataset_ko_rho[which(dataset_ko_rho$cat == "Scrambled"), "cat"] <- "Shuffled\nASVs"
  dataset_ko_rho <- dataset_ko_rho[-which(dataset_ko_rho$cat %in% c("NSTI=1.5",  "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]
  
  dataset_ko_rho$cat <- factor(dataset_ko_rho$cat,
                               levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "Shuffled\nASVs", "NSTI=2", "NSTI=1", "NSTI=0.05"))

  dataset_ko_rho_melt <- melt(dataset_ko_rho)

  dataset_ko_rho_melt[which(dataset_ko_rho_melt$cat == "Shuffled\nASVs"), "Database"] <- "PICRUSt2"
  
  KO_spearman_boxplots[[d]] <- ggplot(dataset_ko_rho_melt, aes(x=cat, y=value, fill=Database)) +
                                                geom_boxplot(outlier.shape = NA) +
                                                geom_quasirandom(size=0.1) +
                                                ylim(c(0.485, 1.31)) +
                                                ylab(c("Spearman Correlation Coefficient")) +
                                                xlab("") +
                                                facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                      axis.text.x=element_text(angle=45, hjust=1),
                                                      legend.position = c(0.13, 0.85), legend.background = element_rect(color = "black", 
                                                                                                                        fill = "white", size = 0.3, linetype = "solid"),
                                                      legend.title = element_text(colour="black", size=8, face="bold"),
                                                      legend.text = element_text(colour="black", size=8)) +
                                                scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
                                                stat_pvalue_manual(ko_rho_wilcoxon[[d]], label = "p_symbol")

}


pdf(file = "../../../figures/Supp_KO_spearman_w_NSTI_and_scrambled.pdf", width=16, height=16)

plot_grid(KO_spearman_boxplots[["cameroon"]],
          KO_spearman_boxplots[["hmp"]],
          KO_spearman_boxplots[["indian"]],
          KO_spearman_boxplots[["mammal"]],
          KO_spearman_boxplots[["ocean"]],
          KO_spearman_boxplots[["primate"]],
          KO_spearman_boxplots[["blueberry"]],
          labels=c("a", "b", "c", "d", "e", "f", "g"),
          nrow=3,
          ncol=3)
dev.off()
