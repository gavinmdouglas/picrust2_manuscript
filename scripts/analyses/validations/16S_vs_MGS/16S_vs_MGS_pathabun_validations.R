### Comparing predicted MetaCyc pathway abundances based on 16S to MGS "gold standard".
### Also compared to Paprica, which output EC number predictions.

rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

compute_pathabun_validation_metrics <- function(dataset_infiles, in_db, out_prefix=NULL, save_RDS=FALSE) {

  if (save_RDS && is.null(out_prefix)) {
    stop("Out prefix options needs to be specificed if saving results to RDS")
  }
  
  # Generate random tables based on subsampling database to calculate null distributions.
  dataset_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = in_db, tab = dataset_infiles$all_pathabun$mgs_pathabun)
  
  dataset_pathabun_mgs_null_df_round <- round(dataset_pathabun_mgs_null_df - 0.00000001)
  
  # Get spearman correlations for each table with MGS:
  dataset_pathabun_mgs_null <- cor_all_cols(tab1 = dataset_pathabun_mgs_null_df, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="Null", metric="spearman")
  dataset_pathabun_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti2, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2", metric="spearman")
  dataset_pathabun_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti1.5, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1.5", metric="spearman")
  dataset_pathabun_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti1, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1", metric="spearman")
  dataset_pathabun_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.5, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.5", metric="spearman")
  dataset_pathabun_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.25, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.25", metric="spearman")
  dataset_pathabun_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.1, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.1", metric="spearman")
  dataset_pathabun_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.05, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.05", metric="spearman")
  dataset_pathabun_picrust2_scrambled_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_pathabun$picrust2_pathabun_scrambled_mean, tab2 = dataset_infiles$all_pathabun$mgs_pathabun, cat_string="Scrambled", metric="spearman")
  
  # Make combined dfs of spearman correlation coefficient subsets of interest:
  dataset_pathabun_spearman_df <- rbind(dataset_pathabun_mgs_null,
                                        dataset_pathabun_picrust2_nsti2_vs_mgs,
                                        dataset_pathabun_picrust2_nsti1.5_vs_mgs,
                                        dataset_pathabun_picrust2_nsti1_vs_mgs,
                                        dataset_pathabun_picrust2_nsti0.5_vs_mgs,
                                        dataset_pathabun_picrust2_nsti0.25_vs_mgs,
                                        dataset_pathabun_picrust2_nsti0.1_vs_mgs,
                                        dataset_pathabun_picrust2_nsti0.05_vs_mgs,
                                        dataset_pathabun_picrust2_scrambled_vs_mgs)
  
  dataset_picrust2_pathabun_nsti2_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti2, category="NSTI=2")
  dataset_picrust2_pathabun_nsti1.5_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti1.5, category="NSTI=1.5")
  dataset_picrust2_pathabun_nsti1_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti1, category="NSTI=1")
  dataset_picrust2_pathabun_nsti0.5_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.5, category="NSTI=0.5")
  dataset_picrust2_pathabun_nsti0.25_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.25, category="NSTI=0.25")
  dataset_picrust2_pathabun_nsti0.1_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.1, category="NSTI=0.1")
  dataset_picrust2_pathabun_nsti0.05_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_nsti0.05, category="NSTI=0.05")
  dataset_null_pathabun_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_pathabun_mgs_null_df_round, category="Null")
  dataset_picrust2_pathabun_scrambled_metrics <- calc_accuracy_metrics(dataset_infiles$all_pathabun$mgs_pathabun, dataset_infiles$all_pathabun$picrust2_pathabun_scrambled_median, category="Scrambled")
  
  dataset_pathabun_acc_df <- rbind(dataset_null_pathabun_metrics,
                                   dataset_picrust2_pathabun_nsti2_metrics,
                                   dataset_picrust2_pathabun_nsti1.5_metrics, 
                                   dataset_picrust2_pathabun_nsti1_metrics,
                                   dataset_picrust2_pathabun_nsti0.5_metrics,
                                   dataset_picrust2_pathabun_nsti0.25_metrics,
                                   dataset_picrust2_pathabun_nsti0.1_metrics,
                                   dataset_picrust2_pathabun_nsti0.05_metrics,
                                   dataset_picrust2_pathabun_scrambled_metrics)
  
  
  dataset_pathabun_acc_df$category <- factor(dataset_pathabun_acc_df$category, levels=c("Null", "NSTI=2", "NSTI=1.5", "NSTI=1",
                                                                                        "NSTI=0.5","NSTI=0.25", "NSTI=0.1", "NSTI=0.05", "Scrambled"))
  
  if(save_RDS) {
    saveRDS(object = dataset_pathabun_spearman_df, file = paste(out_prefix, "_pathabun_spearman_df.rds", sep=""))
    saveRDS(object = dataset_pathabun_acc_df, file = paste(out_prefix, "_pathabun_acc_df.rds", sep=""))
  }
  
  return(list("null_df"=dataset_pathabun_mgs_null_df,
              "spearman_df"=dataset_pathabun_spearman_df,
              "acc_df"=dataset_pathabun_acc_df))
}


pathabun_db <- data.frame(t(read.table(gzfile("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_16S_pathway/path_abun_unstrat.tsv"),
                                    row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

### Loop over all dataset names: read in predictions (restrict to overlapping samples only,
### and get subsets with all possible ECs that overlap across tools filled in) and compute performance metrics.
datasets <- c("hmp", "mammal", "ocean", "blueberry", "indian", "cameroon", "primate")

pathabun_metrics_out <- list()

for(dataset in datasets) {
  pathabun_metrics_out[[dataset]] <- list()
  pathabun_metrics_out[[dataset]][["infiles"]] <- read_in_pathway_predictions(dataset)
  
  pathabun_metrics_out[[dataset]][["metrics"]] <- compute_pathabun_validation_metrics(dataset_infiles = pathabun_metrics_out[[dataset]][["infiles"]],
                                                                                      in_db = pathabun_db,
                                                                                      save_RDS = TRUE,
                                                                                      out_prefix = paste("../saved_RDS/16S_vs_MGS_metrics/", dataset, sep=""))
}

