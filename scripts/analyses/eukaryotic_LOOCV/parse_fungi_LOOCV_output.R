library(Hmisc)

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")
setwd("/home/gavin/projects/picrust_pipeline/fungal_genomes/LOOCV")

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
# Note that the num_subsets is the number of lines in the "groupings" files like fungi_18S_taxa_groupings.tsv
LOOCV_out_18S <- make_tabular_LOOCV_metrics(outfolder="fungi_18S",
                                            num_subsets = 761,
                                            reffile_path="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_18S_counts.txt")

LOOCV_out_18S$level <- capitalize(LOOCV_out_18S$level)

LOOCV_out_18S$level <- factor(LOOCV_out_18S$level, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Assembly"))

# Get null metrics based on comparing each individual assembly with database.
ec_18S <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_18S_counts.txt",
                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)
assembly_18S <- data.frame(t(ec_18S), check.names = FALSE)
ec_18S_null_df <- generate_null_mean_db_funcs(db = ec_18S, tab = assembly_18S)
ec_18S_null_df_round <- round(ec_18S_null_df  - 0.00000001)

ec_18S_null_assembly_scc <- cor_all_cols(tab1 = ec_18S_null_df, tab2 = assembly_18S, cat_string="Assembly Null", metric="spearman")
ec_18S_null_assembly_metrics <- calc_accuracy_metrics(assembly_18S, ec_18S_null_df_round, category="Assembly Null")

rownames(ec_18S_null_assembly_scc) <- as.character(ec_18S_null_assembly_scc$sample_names)
rownames(ec_18S_null_assembly_metrics) <- as.character(ec_18S_null_assembly_metrics$sample)

ec_18S_null_out <- cbind(ec_18S_null_assembly_scc, ec_18S_null_assembly_metrics)
ec_18S_null_out <- ec_18S_null_out[, c("cat", "sample", "metric", "acc", "precision", "recall")]
ec_18S_null_out$mean_nsti <- NA
colnames(ec_18S_null_out) <- c("level", "taxon", "mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")
rownames(ec_18S_null_out) <- NULL

LOOCV_out_18S_w_null <- rbind(LOOCV_out_18S, ec_18S_null_out)

saveRDS(LOOCV_out_18S_w_null, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/18S_LOOCV_metrics.rds")



# Also run for ITS:
LOOCV_out_ITS <- make_tabular_LOOCV_metrics(outfolder="fungi_ITS",
                                            num_subsets = 672,
                                            reffile_path="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_ITS_counts.txt")

LOOCV_out_ITS$level <- capitalize(LOOCV_out_ITS$level)

LOOCV_out_ITS$level <- factor(LOOCV_out_ITS$level, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Assembly"))

# Get null metrics based on comparing each individual assembly with database.
ec_ITS <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_ITS_counts.txt",
                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)
assembly_ITS <- data.frame(t(ec_ITS), check.names = FALSE)
ec_ITS_null_df <- generate_null_mean_db_funcs(db = ec_ITS, tab = assembly_ITS)
ec_ITS_null_df_round <- round(ec_ITS_null_df  - 0.00000001)

ec_ITS_null_assembly_scc <- cor_all_cols(tab1 = ec_ITS_null_df, tab2 = assembly_ITS, cat_string="Assembly Null", metric="spearman")
ec_ITS_null_assembly_metrics <- calc_accuracy_metrics(assembly_ITS, ec_ITS_null_df_round, category="Assembly Null")

rownames(ec_ITS_null_assembly_scc) <- as.character(ec_ITS_null_assembly_scc$sample_names)
rownames(ec_ITS_null_assembly_metrics) <- as.character(ec_ITS_null_assembly_metrics$sample)

ec_ITS_null_out <- cbind(ec_ITS_null_assembly_scc, ec_ITS_null_assembly_metrics)
ec_ITS_null_out <- ec_ITS_null_out[, c("cat", "sample", "metric", "acc", "precision", "recall")]
ec_ITS_null_out$mean_nsti <- NA
colnames(ec_ITS_null_out) <- c("level", "taxon", "mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")
rownames(ec_ITS_null_out) <- NULL

LOOCV_out_ITS_w_null <- rbind(LOOCV_out_ITS, ec_ITS_null_out)

saveRDS(LOOCV_out_ITS_w_null, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/ITS_LOOCV_metrics.rds")
