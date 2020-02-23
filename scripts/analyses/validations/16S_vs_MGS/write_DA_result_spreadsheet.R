### Commands to write the DA results to an excel spreadsheet, which is easier to explore.

rm(list=ls(all.names=TRUE))

library(openxlsx)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/DA_concordance/")

ko_da_out <- readRDS(file = "ko_da_out.rds")
pathabun_da_out <- readRDS(file = "pathabun_da_out.rds")

perf_list_to_dataframe <- function(in_list) {
  
  for(dataset in names(in_list)) {
    orig_colnames <- colnames(in_list[[dataset]])
    in_list[[dataset]][, "dataset"] <- dataset
    in_list[[dataset]][, "metric"] <- rownames(in_list[[dataset]])
    in_list[[dataset]] <- in_list[[dataset]][, c("dataset", "metric", orig_colnames)]
  }
   
  combined_df <- do.call("rbind", in_list)
  rownames(combined_df) <- NULL
  
  return(combined_df)
}

tmp <- ko_da_out$wilcoxon_musicc_out_perf_0.05)

full_DA <- list(KO_wilcoxon_MUSICC_normalized=perf_list_to_dataframe(ko_da_out$wilcoxon_musicc_out_perf_0.05),
                KO_ALDEx2=perf_list_to_dataframe(ko_da_out$aldex2_out_perf_0.05),
                KO_DESeq2=perf_list_to_dataframe(ko_da_out$deseq2_out_perf_0.05),
                KO_DESeq2_GMeans=perf_list_to_dataframe(ko_da_out$deseq2_GMmeans_out_perf_0.05),
                KO_prevalence=perf_list_to_dataframe(ko_da_out$fishers_out_perf_0.05),
                
                pathabun_wilcoxon=perf_list_to_dataframe(pathabun_da_out$wilcoxon_out_perf_0.05),
                pathabun_ALDEx2=perf_list_to_dataframe(pathabun_da_out$aldex2_out_perf_0.05),
                pathabun_DESeq2=perf_list_to_dataframe(pathabun_da_out$deseq2_out_perf_0.05),
                pathabun_DESeq2_GMeans=perf_list_to_dataframe(pathabun_da_out$deseq2_GMmeans_out_perf_0.05),
                pathabun_prevalence=perf_list_to_dataframe(pathabun_da_out$fishers_out_perf_0.05))
                
write.xlsx(full_DA, file = "../../DA_validation_summary.xlsx")

