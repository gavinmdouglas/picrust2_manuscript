### Comparing predicted functions based on 16S to MGS "gold standard".
### Comparisons made between KEGG orthologs predicted by each tool to MGS.
### RDS files saved for spearman correlations and accuracy metrics.

rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

compute_ko_validation_metrics <- function(dataset_infiles, in_db, ko_subset=NULL, out_prefix=NULL, save_RDS=FALSE) {
  
  if (save_RDS && is.null(out_prefix)) {
    stop("Out prefix options needs to be specificed if saving results to RDS")
  }
  
  if (! is.null(ko_subset)) {
    ko_to_keep <- ko_subset[which(ko_subset %in% rownames(dataset_infiles$all_kos_overlap$mgs_ko))]
    
    for(category in names(dataset_infiles$all_kos_overlap)) {
      dataset_infiles$all_kos_overlap[[category]] <- dataset_infiles$all_kos_overlap[[category]][ko_to_keep, ]
    }
  }
  
  # Generate random tables based on subsampling database to calculate null distributions.
  dataset_ko_mgs_null_df <- generate_null_mean_db_funcs(db = in_db, tab = dataset_infiles$all_kos_overlap$mgs_ko)
  
  dataset_ko_mgs_null_df_round <- round(dataset_ko_mgs_null_df - 0.00000001)
  
  # Get spearman correlations for each table with MGS:
  dataset_ko_mgs_null <- cor_all_cols(tab1 = dataset_ko_mgs_null_df, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="Null", metric="spearman")
  dataset_ko_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti2, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2", metric="spearman")
  dataset_ko_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti1.5, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1.5", metric="spearman")
  dataset_ko_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti1, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1", metric="spearman")
  dataset_ko_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.5, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.5", metric="spearman")
  dataset_ko_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.25, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.25", metric="spearman")
  dataset_ko_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.1, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.1", metric="spearman")
  dataset_ko_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.05, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.05", metric="spearman")
  dataset_ko_picrust1_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$picrust1_ko, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="PICRUSt1", metric="spearman")
  dataset_ko_panfp_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$panfp_ko, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="PanFP", metric="spearman")
  dataset_ko_piphillin_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$piphillin_ko, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="Piphillin", metric="spearman")
  dataset_ko_tax4fun2_vs_mgs <- cor_all_cols(tab1 = dataset_infiles$all_kos_overlap$tax4fun2_ko, tab2 = dataset_infiles$all_kos_overlap$mgs_ko, cat_string="Tax4Fun2", metric="spearman")
  
  # Make combined dfs of spearman correlation coefficient subsets of interest:
  dataset_ko_spearman_df <- rbind(dataset_ko_mgs_null,
                                  dataset_ko_tax4fun2_vs_mgs,
                                  dataset_ko_panfp_vs_mgs,
                                  dataset_ko_piphillin_vs_mgs,
                                  dataset_ko_picrust1_vs_mgs,
                                  dataset_ko_picrust2_nsti2_vs_mgs,
                                  dataset_ko_picrust2_nsti1.5_vs_mgs,
                                  dataset_ko_picrust2_nsti1_vs_mgs,
                                  dataset_ko_picrust2_nsti0.5_vs_mgs,
                                  dataset_ko_picrust2_nsti0.25_vs_mgs,
                                  dataset_ko_picrust2_nsti0.1_vs_mgs,
                                  dataset_ko_picrust2_nsti0.05_vs_mgs)
  
  dataset_picrust2_ko_nsti2_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti2, category="NSTI=2")
  dataset_picrust2_ko_nsti1.5_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti1.5, category="NSTI=1.5")
  dataset_picrust2_ko_nsti1_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti1, category="NSTI=1")
  dataset_picrust2_ko_nsti0.5_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.5, category="NSTI=0.5")
  dataset_picrust2_ko_nsti0.25_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.25, category="NSTI=0.25")
  dataset_picrust2_ko_nsti0.1_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.1, category="NSTI=0.1")
  dataset_picrust2_ko_nsti0.05_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust2_ko_nsti0.05, category="NSTI=0.05")
  dataset_picrust1_ko_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$picrust1_ko, category="PICRUSt1")
  dataset_piphillin_ko_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$piphillin, category="Piphillin")
  dataset_panfp_ko_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$panfp, category="PanFP")
  dataset_tax4fun2_ko_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_infiles$all_kos_overlap$tax4fun2, category="Tax4Fun2")
  dataset_null_ko_metrics <- calc_accuracy_metrics(dataset_infiles$all_kos_overlap$mgs_ko, dataset_ko_mgs_null_df_round, category="Null")
  
  dataset_ko_acc_df <- rbind(dataset_null_ko_metrics,
                             dataset_tax4fun2_ko_metrics,
                             dataset_panfp_ko_metrics,
                             dataset_piphillin_ko_metrics,
                             dataset_picrust1_ko_metrics,
                             dataset_picrust2_ko_nsti2_metrics,
                             dataset_picrust2_ko_nsti1.5_metrics, 
                             dataset_picrust2_ko_nsti1_metrics,
                             dataset_picrust2_ko_nsti0.5_metrics,
                             dataset_picrust2_ko_nsti0.25_metrics,
                             dataset_picrust2_ko_nsti0.1_metrics,
                             dataset_picrust2_ko_nsti0.05_metrics)
  
  
  dataset_ko_acc_df$category <- factor(dataset_ko_acc_df$category, levels=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2",
                                                                            "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
  
  
  if(save_RDS) {
    saveRDS(object = dataset_ko_spearman_df, file = paste(out_prefix, "_ko_spearman_df.rds", sep=""))
    saveRDS(object = dataset_ko_acc_df, file = paste(out_prefix, "_ko_acc_df.rds", sep=""))
  }
  
  return(list("null_df"=dataset_ko_mgs_null_df,
              "spearman_df"=dataset_ko_spearman_df,
              "acc_df"=dataset_ko_acc_df))
}

### Read in database files (used for calculating null distribution).
ko_db <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)

### Loop over all dataset names: read in predictions (restrict to overlapping samples only,
### and get subsets with all possible KOs that overlap across tools filled in) and compute performance metrics.

datasets <- c("hmp", "mammal", "ocean", "blueberry", "crossbiome", "mat", "indian", "cameroon", "primate")

for(dataset in datasets) {
  tmp <- compute_ko_validation_metrics(dataset_infiles = read_in_ko_predictions(dataset),
                                       in_db = ko_db,
                                       save_RDS = TRUE,
                                       out_prefix = paste("../saved_RDS/16S_vs_MGS_metrics/", dataset, sep=""))
}
