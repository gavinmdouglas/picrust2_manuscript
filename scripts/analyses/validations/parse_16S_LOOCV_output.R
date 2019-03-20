library(Hmisc)

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")
setwd("/home/gavin/projects/picrust_pipeline/16S_LOOCV/")

make_tabular_LOOCV_metrics <- function(outfolder, num_subsets, reffile_path) {
  
  LOOCV_out <- data.frame(matrix(NA, nrow=num_subsets, ncol=7))
  colnames(LOOCV_out) <- c("level", "taxon", "mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")
  
  reffile <- data.frame(t(read.table(reffile_path, header=T, sep="\t", row.names=1)))
  
  LOOCV_subdirs <- list.dirs(path = outfolder, recursive = FALSE)
  
  rep_i = 1
  
  for (subdir in LOOCV_subdirs) {
    LOOCV_outfiles <- list.files(subdir, full.names = TRUE)
    
    subdir_clean <- sub(".*/", "", subdir)
    
    for(outfile in LOOCV_outfiles) {
      
      outfile_clean <- sub("_loocv.txt", "", outfile)
      outfile_clean <- sub(".*/", "", outfile_clean)
      
      loocv_pred <- data.frame(t(read.table(outfile, header=T, sep="\t", row.names=1)))
      
      # Remove NSTI row.
      nsti_row_i <- which(rownames(loocv_pred) == "metadata_NSTI")
      nsti_row <- as.numeric(loocv_pred[nsti_row_i, ])
      loocv_pred <- loocv_pred[-nsti_row_i, , drop=FALSE]
      
      reffile_subset <- reffile[, colnames(loocv_pred), drop=FALSE]
      
      rho_out <- cor_all_cols(tab1 = loocv_pred, tab2 = reffile_subset, cat_string="LOOCV", metric="spearman")
      
      acc_out <- calc_accuracy_metrics(reffile_subset, loocv_pred, category="LOOCV")
      
      LOOCV_out[rep_i, c("level", "taxon")] <- c(subdir_clean, outfile_clean)
      LOOCV_out[rep_i, c("mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")] <- c(mean(rho_out$metric),
                                                                                                     mean(acc_out$acc),
                                                                                                     mean(acc_out$precision),
                                                                                                     mean(acc_out$recall),
                                                                                                     mean(nsti_row))
      rep_i = rep_i + 1
    }
    
  }
  return(LOOCV_out)
}

# Get table of all metrics per subset.
# Note that the num_subsets is the number of lines in the "groupings" files like 16S_taxa_groupings_rand100.tsv
LOOCV_out_16S <- make_tabular_LOOCV_metrics(outfolder="LOOCV_out/",
                                            num_subsets = 664,
                                            reffile_path="/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/untracked/ko.txt")

#LOOCV_out_16S_saved <- LOOCV_out_16S

LOOCV_out_16S$level <- capitalize(LOOCV_out_16S$level)

LOOCV_out_16S$level <- factor(LOOCV_out_16S$level, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Assembly"))

# Get null metrics based on comparing each individual assembly with database.
ko_16S <- read.table("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/untracked/ko.txt",
                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)
assembly_16S <- data.frame(t(ko_16S), check.names = FALSE)
ko_16S_null_df <- generate_null_mean_db_funcs(db = ko_16S, tab = assembly_16S)
ko_16S_null_df_round <- round(ko_16S_null_df  - 0.00000001)

ko_16S_null_assembly_scc <- cor_all_cols(tab1 = ko_16S_null_df, tab2 = assembly_16S, cat_string="Assembly Null", metric="spearman")
ko_16S_null_assembly_metrics <- calc_accuracy_metrics(assembly_16S, ko_16S_null_df_round, category="Assembly Null")

rownames(ko_16S_null_assembly_scc) <- as.character(ko_16S_null_assembly_scc$sample_names)
rownames(ko_16S_null_assembly_metrics) <- as.character(ko_16S_null_assembly_metrics$sample)

ko_16S_null_out <- cbind(ko_16S_null_assembly_scc, ko_16S_null_assembly_metrics)
ko_16S_null_out <- ko_16S_null_out[, c("cat", "sample", "metric", "acc", "precision", "recall")]
ko_16S_null_out$mean_nsti <- NA
colnames(ko_16S_null_out) <- c("level", "taxon", "mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")
rownames(ko_16S_null_out) <- NULL

LOOCV_out_16S_w_null <- rbind(LOOCV_out_16S, ko_16S_null_out)

saveRDS(LOOCV_out_16S_w_null, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_LOOCV_metrics.rds")

