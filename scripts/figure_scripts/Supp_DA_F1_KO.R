### Code to make figures contasting F1 scores based on differential abundance tests across different statistical methods

rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/")

source("../../../scripts/picrust2_ms_functions.R")

# First for KOs
parse_perf_to_df <- function(in_perf, method_name) {
  
  in_perf$hmp$metric <- rownames(in_perf$hmp)
  in_perf_hmp <- melt(in_perf$hmp)
  in_perf_hmp$dataset <- "HMP\n22 supragingival plaque vs 36 tongue dorsum"
  
  in_perf$cameroon$metric <- rownames(in_perf$cameroon)
  in_perf_cameroon <- melt(in_perf$cameroon)
  in_perf_cameroon$dataset <- "Cameroonian\n19 parasite+ vs 36 parasite- individuals"
  
  in_perf$indian$metric <- rownames(in_perf$indian)
  in_perf_indian <- melt(in_perf$indian)
  in_perf_indian$dataset <- "Indian\n51 individuals from Bhopal vs 38 from Kerala"
  
  in_perf$primate$metric <- rownames(in_perf$primate)
  in_perf_primate <- melt(in_perf$primate)
  in_perf_primate$dataset <- "Primate\n29 old world vs 29 new world monkeys"
  
  in_perf_combined <- rbind(in_perf_hmp, in_perf_cameroon,
                            in_perf_indian, in_perf_primate)
  
  in_perf_combined$variable <- as.character(in_perf_combined$variable)
  
  in_perf_combined$method <- method_name
  
  
  categories2remove <- c("picrust2_ko_scrambled_median", "picrust2_ko_nsti0.05", "picrust2_ko_nsti0.1",
                         "picrust2_ko_nsti0.25", "picrust2_ko_nsti0.5", "picrust2_ko_nsti1", "picrust2_ko_nsti1.5")
  if(length(which(in_perf_combined$variable %in% categories2remove)) > 0) {
    in_perf_combined <- in_perf_combined[-which(in_perf_combined$variable %in% categories2remove), ]
  }

  in_perf_combined$variable[which(in_perf_combined$variable == "mgs_ko_alt")] <- "Alt. MGS"
  in_perf_combined$variable[which(in_perf_combined$variable == "picrust1_ko")] <- "PICRUSt1"
  in_perf_combined$variable[which(in_perf_combined$variable == "picrust2_ko_nsti2")] <- "PICRUSt2"
  in_perf_combined$variable[which(in_perf_combined$variable == "picrust2_ko_scrambled_mean")] <- "Shuffled\nASVs"
  in_perf_combined$variable[which(in_perf_combined$variable == "piphillin_ko")] <- "Piphillin"
  
  other_counter = 0
  if("panfp_ko" %in% in_perf_combined$variable) {
    in_perf_combined$variable[which(in_perf_combined$variable == "panfp_ko")] <- "PanFP"
    other_counter = other_counter + 1
  }
  
  if("tax4fun2_ko" %in% in_perf_combined$variable) {
    in_perf_combined$variable[which(in_perf_combined$variable == "tax4fun2_ko")] <- "Tax4Fun2"
    other_counter = other_counter + 1
  }
  
  if(other_counter > 0) {
  in_perf_combined$variable <- factor(in_perf_combined$variable,
                                       levels=c("Shuffled\nASVs", "Alt. MGS", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))
  } else {
    in_perf_combined$variable <- factor(in_perf_combined$variable,
                                        levels=c("Shuffled\nASVs", "Alt. MGS", "Piphillin", "PICRUSt1", "PICRUSt2"))
  }
  
  return(in_perf_combined)
}
  
ko_da_out <- readRDS(file="../DA_concordance/ko_da_out.rds")

wilcoxon_out_musicc_perf_0.05_combined <- parse_perf_to_df(ko_da_out$wilcoxon_musicc_out_perf_0.05, "Wilcoxon")
aldex2_perf_0.05_combined <- parse_perf_to_df(ko_da_out$aldex2_out_perf_0.05, "ALDEx2")
deseq2_perf_0.05_combined <- parse_perf_to_df(ko_da_out$deseq2_out_perf_0.05, "DESeq2")
deseq2_GMeans_perf_0.05_combined <- parse_perf_to_df(ko_da_out$deseq2_GMmeans_out_perf_0.05, "DESeq2 (GMeans)")
fishers_perf_0.05_combined <- parse_perf_to_df(ko_da_out$fishers_out_perf_0.05, "Fisher's Exact (Prevalence)")

wilcoxon_out_musicc_perf_0.05_combined_f1 <- wilcoxon_out_musicc_perf_0.05_combined[which(wilcoxon_out_musicc_perf_0.05_combined$metric == "f1"), ]
aldex2_perf_0.05_combined_f1 <- aldex2_perf_0.05_combined[which(aldex2_perf_0.05_combined$metric == "f1"), ]
deseq2_perf_0.05_combined_f1 <- deseq2_perf_0.05_combined[which(deseq2_perf_0.05_combined$metric == "f1"), ]
deseq2_GMeans_perf_0.05_combined_f1 <- deseq2_GMeans_perf_0.05_combined[which(deseq2_GMeans_perf_0.05_combined$metric == "f1"), ]
fishers_perf_0.05_combined_f1 <- fishers_perf_0.05_combined[which(fishers_perf_0.05_combined$metric == "f1"), ]

ko_f1_dfs <- list(Wilcoxon=wilcoxon_out_musicc_perf_0.05_combined_f1,
                  ALDEx2=aldex2_perf_0.05_combined_f1,
                  DESeq2=deseq2_perf_0.05_combined_f1,
                  DESeq2_gmeans=deseq2_GMeans_perf_0.05_combined_f1,
                  Fisher=fishers_perf_0.05_combined_f1)


ko_f1_titles <- list(Wilcoxon="Wilcoxon test",
                      ALDEx2="ALDEx2",
                      DESeq2="DESeq2",
                      DESeq2_gmeans="DESeq2 (GMeans)",
                      Fisher="Fisher's exact test (prevalence)")

ko_f1_plots <- list()


for(test in names(ko_f1_dfs)) {
  ko_f1_plots[[test]] <- ggplot(ko_f1_dfs[[test]],
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
    ggtitle(ko_f1_titles[[test]])
  
  if(length(which(is.na(ko_f1_dfs[[test]]$value))) > 0) {
    
    na_df <- ko_f1_dfs[[test]][which(is.na(ko_f1_dfs[[test]]$value)), ]
    na_df$group1 <- na_df$variable
    na_df$group2 <- na_df$variable
    na_df$y.position <- 0.05
    na_df$p_symbol <- "NA"
  
    ko_f1_plots[[test]] <- ko_f1_plots[[test]] + 
    stat_pvalue_manual(data = na_df, label = "p_symbol",
                       bracket.size = 0, tip.length = 0, label.size = 2.5)
    
      
  }
  
  
  
}

pdf(file = "../../../figures/Supp_DA_F1_KO.pdf", width=12, height=9)

plot_grid(ko_f1_plots$Wilcoxon, ko_f1_plots$ALDEx2, ko_f1_plots$DESeq2, ko_f1_plots$Fisher, nrow=4)

dev.off()

