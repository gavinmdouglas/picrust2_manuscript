### Code to make figure contrasting pathway abundance correlations on each 16S validation dataset.
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
                                                             y_pos_start = 0.94,
                                                             dist_to_add=0.05)
  pathabun_rho[[d]] <- pathabun_rho_outlist_raw[[d]][[1]]
  pathabun_rho_wilcoxon[[d]] <- pathabun_rho_outlist_raw[[d]][[2]]
  pathabun_rho_wilcoxon[[d]][which(pathabun_rho_wilcoxon[[d]]$group2 == "Scrambled"), "group2"] <- "Shuffled\nASVs"
}


# Make plot for each dataset.
pathabun_spearman_boxplots <- list()

for(j in 1:length(datasets)) {
  
  d <- datasets[j]
  
  dataset_pathabun_rho <- pathabun_rho[[d]]
  
  dataset_pathabun_rho$cat <- as.character(dataset_pathabun_rho$cat)
  dataset_pathabun_rho[which(dataset_pathabun_rho$cat == "Scrambled"), "cat"] <- "Shuffled\nASVs"
  dataset_pathabun_rho <- dataset_pathabun_rho[-which(dataset_pathabun_rho$cat %in% c("NSTI=1.5",  "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]
  
  dataset_pathabun_rho$cat <- factor(dataset_pathabun_rho$cat,
                               levels=c("Null", "Shuffled\nASVs", "NSTI=2", "NSTI=1", "NSTI=0.05"))
  
  dataset_pathabun_rho_melt <- melt(dataset_pathabun_rho)
  
  dataset_pathabun_rho_melt[which(dataset_pathabun_rho_melt$cat == "Shuffled\nASVs"), "Database"] <- "PICRUSt2"
  
  pathabun_spearman_boxplots[[d]] <- ggplot(dataset_pathabun_rho_melt, aes(x=cat, y=value, fill=Database)) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(size=0.1) +
    ylim(c(0.6, 1.1)) +
    ylab(c("Spearman Correlation Coefficient")) +
    xlab("") +
    facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=45, hjust=1),
          legend.position = "none", legend.background = element_rect(color = "black", 
                                                                            fill = "white", size = 0.3, linetype = "solid"),
          legend.title = element_text(colour="black", size=8, face="bold"),
          legend.text = element_text(colour="black", size=8)) +
    scale_fill_manual(values=c("light grey", "#00BFC4")) +
    stat_pvalue_manual(pathabun_rho_wilcoxon[[d]], label = "p_symbol")
  
}


pdf(file = "../../../figures/Supp_pathabun_spearman_w_NSTI_and_scrambled.pdf", width=12, height=14)

plot_grid(pathabun_spearman_boxplots[["cameroon"]],
          pathabun_spearman_boxplots[["hmp"]],
          pathabun_spearman_boxplots[["indian"]],
          pathabun_spearman_boxplots[["mammal"]],
          pathabun_spearman_boxplots[["ocean"]],
          pathabun_spearman_boxplots[["primate"]],
          pathabun_spearman_boxplots[["blueberry"]],
          labels=c("a", "b", "c", "d", "e", "f", "g"),
          nrow=3,
          ncol=3)
dev.off()


