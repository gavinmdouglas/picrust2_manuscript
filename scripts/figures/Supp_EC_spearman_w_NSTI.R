### Code to make figure contrasting EC correlations on each 16S validation dataset.
### Include all NSTI cut-offs in these plots.

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggbeeswarm)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

extra_nsti_categories <- c("NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")

datasets <- c("hmp", "mammal", "ocean", "blueberry", "indian", "cameroon", "primate")
dataset_names <- c("HMP", "Mammal", "Ocean", "Soil (Blueberry)", "India", "Cameroon", "Primate")

ec_rho_outlist <- list()
ec_rho <- list()
ec_rho_wilcoxon <- list()

for(i in 1:length(datasets)) {
  
  ec_rho_outlist[[datasets[i]]] <- parse_rho_rds_and_calc_wilcoxon(rho_rds = paste(datasets[i], "_ec_spearman_df.rds", sep=""),
                                                                  dataset_name = dataset_names[i],
                                                                  wilcox_cat2ignore = extra_nsti_categories,
                                                                  y_pos_start = 0.9)
  
   ec_rho[[datasets[i]]] <- ec_rho_outlist[[datasets[i]]][[1]]
   ec_rho_wilcoxon[[datasets[i]]] <- ec_rho_outlist[[datasets[i]]][[2]]

}

# Make plot for each dataset.
EC_spearman_boxplots <- list()

for(j in 1:length(datasets)) {
  
  dataset_ec_rho <- ec_rho[[datasets[j]]]
  
  dataset_ec_rho$cat <- as.character(dataset_ec_rho$cat)
  dataset_ec_rho <- dataset_ec_rho[-which(dataset_ec_rho$cat %in% c("NSTI=1.5",  "NSTI=0.5", "NSTI=0.25", "NSTI=0.1")), ]
  
  dataset_ec_rho$cat <- factor(dataset_ec_rho$cat,
                         levels=c("Null", "PAPRICA", "NSTI=2", "NSTI=1", "NSTI=0.05"))

  dataset_ec_rho_melt <- melt(dataset_ec_rho)

  dataset_ec_rho_melt[which(dataset_ec_rho_melt$Database == "Other"), "Database"] <- "PAPRICA"

  EC_spearman_boxplots[[datasets[j]]] <- ggplot(dataset_ec_rho_melt, aes(x=cat, y=value, fill=Database)) +
                                        geom_boxplot(outlier.shape = NA) +
                                        geom_quasirandom(size=0.1) +
                                        scale_y_continuous(breaks=c(0.4, 0.6, 0.8, 1.0), limits=c(0.4, 1)) +
                                        ylab(c("Spearman Correlation Coefficient")) +
                                        xlab("") +
                                        facet_grid(. ~ dataset, scales = "free", space = "free", switch="x") +
                                        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                              axis.text.x=element_text(angle=45, hjust=1),
                                              legend.position = c(0.6, 0.2), legend.background = element_rect(color = "black", 
                                                                                                                fill = "white", size = 0.3, linetype = "solid"),
                                              legend.title = element_text(colour="black", size=8, face="bold"),
                                              legend.text = element_text(colour="black", size=8)) +
                                        scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
                                        stat_pvalue_manual(ec_rho_wilcoxon[[datasets[j]]], label = "p_symbol")

}

pdf(file = "../../../figures/Supp_EC_spearman.pdf", width=15, height=9)

plot_grid(EC_spearman_boxplots[["cameroon"]],
          EC_spearman_boxplots[["hmp"]],
          EC_spearman_boxplots[["indian"]],
          EC_spearman_boxplots[["mammal"]],
          EC_spearman_boxplots[["ocean"]],
          EC_spearman_boxplots[["primate"]],
          EC_spearman_boxplots[["blueberry"]],
          labels=c("a", "b", "c", "d", "e", "f", "g"),
          nrow=2,
          ncol=4)

dev.off()
