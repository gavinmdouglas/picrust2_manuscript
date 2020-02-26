### Code to make figures contasting F1 scores based on differential abundance tests across different statistical methods

rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("../../../scripts/picrust2_ms_functions.R")

parse_perf_to_df_pathabun <- function(in_perf, method_name) {
  
  in_perf$hmp$metric <- rownames(in_perf$hmp)
  in_perf_hmp <- melt(in_perf$hmp)
  in_perf_hmp$dataset <- "HMP subset"
  
  in_perf$cameroon$metric <- rownames(in_perf$cameroon)
  in_perf_cameroon <- melt(in_perf$cameroon)
  in_perf_cameroon$dataset <- "Cameroonian subset"
  
  in_perf$indian$metric <- rownames(in_perf$indian)
  in_perf_indian <- melt(in_perf$indian)
  in_perf_indian$dataset <- "Indian subset"
  
  in_perf$primate$metric <- rownames(in_perf$primate)
  in_perf_primate <- melt(in_perf$primate)
  in_perf_primate$dataset <- "Primate subset"
  
  in_perf_combined <- rbind(in_perf_hmp, in_perf_cameroon,
                            in_perf_indian, in_perf_primate)
  
  in_perf_combined$variable <- as.character(in_perf_combined$variable)
  
  in_perf_combined$method <- method_name
  
  
  categories2remove <- c("picrust2_pathabun_scrambled_median", "picrust2_pathabun_nsti0.05", "picrust2_pathabun_nsti0.1",
                         "picrust2_pathabun_nsti0.25", "picrust2_pathabun_nsti0.5", "picrust2_pathabun_nsti1", "picrust2_pathabun_nsti1.5")
  if(length(which(in_perf_combined$variable %in% categories2remove)) > 0) {
    in_perf_combined <- in_perf_combined[-which(in_perf_combined$variable %in% categories2remove), ]
  }

  in_perf_combined$variable[which(in_perf_combined$variable == "picrust2_pathabun_nsti2")] <- "PICRUSt2"
  in_perf_combined$variable[which(in_perf_combined$variable == "picrust2_pathabun_scrambled_mean")] <- "Shuffled\nASVs"
  

  in_perf_combined$variable <- factor(in_perf_combined$variable,
                                      levels=c("Shuffled\nASVs", "PICRUSt2"))

  
  return(in_perf_combined)
}
  
pathabun_da_out <- readRDS(file="../DA_concordance/pathabun_da_out.rds")

wilcoxon_out_perf_0.05_combined <- parse_perf_to_df_pathabun(pathabun_da_out$wilcoxon_out_perf_0.05, "Wilcoxon")
aldex2_perf_0.05_combined <- parse_perf_to_df_pathabun(pathabun_da_out$aldex2_out_perf_0.05, "ALDEx2")
deseq2_perf_0.05_combined <- parse_perf_to_df_pathabun(pathabun_da_out$deseq2_out_perf_0.05, "DESeq2")
deseq2_GMeans_perf_0.05_combined <- parse_perf_to_df_pathabun(pathabun_da_out$deseq2_GMmeans_out_perf_0.05, "DESeq2 (GMeans)")
fishers_perf_0.05_combined <- parse_perf_to_df_pathabun(pathabun_da_out$fishers_out_perf_0.05, "Fisher's Exact (Prevalence)")

wilcoxon_out_perf_0.05_combined_f1 <- wilcoxon_out_perf_0.05_combined[which(wilcoxon_out_perf_0.05_combined$metric == "f1"), ]
aldex2_perf_0.05_combined_f1 <- aldex2_perf_0.05_combined[which(aldex2_perf_0.05_combined$metric == "f1"), ]
deseq2_perf_0.05_combined_f1 <- deseq2_perf_0.05_combined[which(deseq2_perf_0.05_combined$metric == "f1"), ]
deseq2_GMeans_perf_0.05_combined_f1 <- deseq2_GMeans_perf_0.05_combined[which(deseq2_GMeans_perf_0.05_combined$metric == "f1"), ]
fishers_perf_0.05_combined_f1 <- fishers_perf_0.05_combined[which(fishers_perf_0.05_combined$metric == "f1"), ]

pathabun_f1_dfs <- list(Wilcoxon=wilcoxon_out_perf_0.05_combined_f1,
                  ALDEx2=aldex2_perf_0.05_combined_f1,
                  DESeq2=deseq2_perf_0.05_combined_f1,
                  DESeq2_gmeans=deseq2_GMeans_perf_0.05_combined_f1,
                  Fisher=fishers_perf_0.05_combined_f1)


pathabun_f1_titles <- list(Wilcoxon="Wilcoxon test (relative abundances)",
                      ALDEx2="ALDEx2",
                      DESeq2="DESeq2",
                      DESeq2_gmeans="DESeq2 (GMeans)",
                      Fisher="Fisher's exact test (prevalence)")

pathabun_f1_plots <- list()


for(test in names(pathabun_f1_dfs)) {
  pathabun_f1_plots[[test]] <- ggplot(pathabun_f1_dfs[[test]],
                                aes(x=variable, y=value)) +
    geom_bar(stat="identity", fill="#009E73", na.rm=TRUE) +
    facet_grid(~ dataset, scales = "free", space = "free", switch="y") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=45, hjust=1),
          legend.position = "none",
          panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(0.75,"lines")) +
    ylab("F1 Score") +
    xlab("") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ggtitle(pathabun_f1_titles[[test]])
  
  if(length(which(is.na(pathabun_f1_dfs[[test]]$value))) > 0) {
    
    na_df <- pathabun_f1_dfs[[test]][which(is.na(pathabun_f1_dfs[[test]]$value)), ]
    na_df$group1 <- na_df$variable
    na_df$group2 <- na_df$variable
    na_df$y.position <- 0.05
    na_df$p_symbol <- "NA"
  
    pathabun_f1_plots[[test]] <- pathabun_f1_plots[[test]] + 
    stat_pvalue_manual(data = na_df, label = "p_symbol",
                       bracket.size = 0, tip.length = 0, label.size = 2.5)
    
      
  }
  
  
  
}

pdf(file = "../../../figures/Supp_DA_F1_pathabun.pdf", width=8, height=9)

plot_grid(pathabun_f1_plots$Wilcoxon, pathabun_f1_plots$ALDEx2, pathabun_f1_plots$DESeq2, pathabun_f1_plots$Fisher, nrow=4)

dev.off()

