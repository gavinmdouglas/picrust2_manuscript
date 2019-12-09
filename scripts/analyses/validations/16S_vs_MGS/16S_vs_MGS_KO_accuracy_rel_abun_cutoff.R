# Generate Spearman and accuracy metrics as before, but instead of saving the results run three analyses:
# (1) Determine the mean and sd in relative abundance of false positive KOs for each approach.
# (2) Determine the fold-difference in FPs and FNs per dataset between PICRUSt2 and Piphillin.
# (3) Calculate Spearman's R value between # of MGS annotated reads (in humann2 output) and per-sample precision.

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
datasets <- c("hmp", "mammal", "ocean", "blueberry", "indian", "cameroon", "primate")

ko_metrics_out <- list()

# Read in and convert all prediction tables to relative abundance
for(dataset in datasets) {
  ko_metrics_out[[dataset]] <- list()
  ko_metrics_out[[dataset]][["infiles"]] <- read_in_ko_predictions(dataset)
  
  for(n in names(ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap)) {
    ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap_rel <- list()
    ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap_rel[[n]] <- data.frame(sweep(ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap[[n]],  2, colSums(ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap[[n]]) , '/')) * 100
    
    ko_metrics_out[[dataset]][["metrics"]] <- compute_ko_validation_metrics(dataset_infiles = ko_metrics_out[[dataset]][["infiles"]],
                                                                            in_db = ko_db,
                                                                            save_RDS = FALSE)
  }
  
}

# Get breakdown of relative abundance of FP KOs.
categories <- c("picrust2_ko_nsti2", "piphillin_ko")
FP_KO_rel_abun <- data.frame(matrix(NA, nrow=7, ncol=5))
colnames(FP_KO_rel_abun) <- c("dataset_mean_picrust2", "dataset_mean_piphillin", "mean_piphillin_vs_mean_picrust2", "sd_mean_piphillin_vs_mean_picrust2", "fold_wilcox_p")
rownames(FP_KO_rel_abun) <- datasets

for(dataset in datasets) {
  
  mgs_table <- ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap_rel[["mgs_ko"]]
  picrust2_table <- ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap_rel[["picrust2_ko_nsti2"]]
  piphillin_table <- ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap_rel[["piphillin_ko"]]

  picrust2_mean_sample_FP_KO_relabun <- c()
  piphillin_mean_sample_FP_KO_relabun <- c()
    
    for(sample in colnames(mgs_table)) {
      true_pos_i <- which(mgs_table[, sample] > 0)

      picrust2_pos_i <- which(picrust2_table[, sample] > 0)
      picrust2_false_pos_i <- picrust2_pos_i[which(! picrust2_pos_i %in% true_pos_i)]

      piphillin_pos_i <- which(piphillin_table[, sample] > 0)
      piphillin_false_pos_i <- piphillin_pos_i[which(! piphillin_pos_i %in% true_pos_i)]
        
      picrust2_mean_sample_FP_KO_relabun <- c(picrust2_mean_sample_FP_KO_relabun, mean(picrust2_table[picrust2_false_pos_i, sample]))
      piphillin_mean_sample_FP_KO_relabun <- c(piphillin_mean_sample_FP_KO_relabun, mean(piphillin_table[piphillin_false_pos_i, sample]))
    }
  
  FP_KO_rel_abun[dataset, "dataset_mean_picrust2"] <- mean(picrust2_mean_sample_FP_KO_relabun)
  FP_KO_rel_abun[dataset, "dataset_mean_piphillin"] <- mean(piphillin_mean_sample_FP_KO_relabun)
  FP_KO_rel_abun[dataset, "mean_piphillin_vs_mean_picrust2"] <- mean(piphillin_mean_sample_FP_KO_relabun / picrust2_mean_sample_FP_KO_relabun)
  FP_KO_rel_abun[dataset, "sd_mean_piphillin_vs_mean_picrust2"] <- sd(piphillin_mean_sample_FP_KO_relabun / picrust2_mean_sample_FP_KO_relabun)
  FP_KO_rel_abun[dataset, "fold_wilcox_p"] <- wilcox.test(piphillin_mean_sample_FP_KO_relabun, picrust2_mean_sample_FP_KO_relabun)$p.value
  
}


# Get breakdown of how many KOs with non-zero relative abundance are below a 0.00001% cutoff
categories <- c("picrust2_ko_nsti2", "picrust1_ko", "panfp_ko", "piphillin_ko", "tax4fun2_ko", "mgs_ko")
categories_mean <- paste("mean", categories, sep="_")
categories_sd <- paste("sd", categories, sep="_")
num_present_KOs_below_cutoff <- data.frame(matrix(NA, nrow=7, ncol=12))
colnames(num_present_KOs_below_cutoff) <- c(categories_mean, categories_sd)
rownames(num_present_KOs_below_cutoff) <- datasets

for(dataset in datasets) {
  for(cat in categories) {
  
    pred_table <- ko_metrics_out[[dataset]][["infiles"]]$all_kos_overlap_rel[[cat]]
    pred_table_percent_below <- c()
    
    for(sample in colnames(pred_table)) {
      positive_i <- which(pred_table[, sample] > 0)
      percentage_below <- (length(which(pred_table[positive_i, sample] < 0.0001)) / length(positive_i)) * 100
      pred_table_percent_below <- c(pred_table_percent_below, percentage_below)
    }

    num_present_KOs_below_cutoff[dataset, paste("mean", cat, sep="_")] <- mean(pred_table_percent_below)
    num_present_KOs_below_cutoff[dataset, paste("sd", cat, sep="_")] <- sd(pred_table_percent_below)
    
  }
}


# Get table of FPs and FNs by PICRUSt2 and Piphillin
false_calls <- data.frame(matrix(NA, nrow=7, ncol=4))
colnames(false_calls) <- c("picrust2_FP", "picrust2_FN", "piphillin_FP", "piphillin_FN")
rownames(false_calls) <- datasets

for(dataset in datasets) {
  false_calls[dataset, "picrust2_FP"] <- mean(ko_metrics_out[[dataset]][["metrics"]]$acc_df[which(ko_metrics_out[[dataset]][["metrics"]]$acc_df$category == "NSTI=2"), "FP"])
  false_calls[dataset, "picrust2_FN"] <- mean(ko_metrics_out[[dataset]][["metrics"]]$acc_df[which(ko_metrics_out[[dataset]][["metrics"]]$acc_df$category == "NSTI=2"), "FN"])
  false_calls[dataset, "piphillin_FP"] <- mean(ko_metrics_out[[dataset]][["metrics"]]$acc_df[which(ko_metrics_out[[dataset]][["metrics"]]$acc_df$category == "Piphillin"), "FP"])
  false_calls[dataset, "piphillin_FN"] <- mean(ko_metrics_out[[dataset]][["metrics"]]$acc_df[which(ko_metrics_out[[dataset]][["metrics"]]$acc_df$category == "Piphillin"), "FN"])
}

false_calls$picrust2_FP_percent <- (false_calls$picrust2_FP / 6220) * 100
false_calls$picrust2_FN_percent <- (false_calls$picrust2_FN / 6220) * 100
false_calls$piphillin_FP_percent <- (false_calls$piphillin_FP / 6220) * 100
false_calls$piphillin_FN_percent <- (false_calls$piphillin_FN / 6220) * 100


# Get mean Spearman rho between MGS annotated depth and PICRUSt2 precision per dataset (along with P value).
# Do this for all combined and just for all combined.

depth_vs_precision <- data.frame(matrix(NA, nrow=7, ncol=2))
colnames(depth_vs_precision) <- c("rho", "p")
rownames(depth_vs_precision) <- datasets

all_precision <- c()
all_depth <- c()

for(dataset in datasets) {
  ko_sums <- colSums(ko_metrics_out[[dataset]]$infiles$all_kos_overlap$mgs_ko)

  picrust2_subset <- ko_metrics_out[[dataset]]$metrics$acc_df[which(ko_metrics_out[[dataset]]$metrics$acc_df$category == "PICRUSt2"), ]
  
  picrust2_subset$mgs_ko_sums <- colSums(ko_metrics_out[[dataset]]$infiles$all_kos_overlap$mgs_ko[, picrust2_subset$sample])
  
  picrust2_spearman <- cor.test(colSums(ko_metrics_out[[dataset]]$infiles$all_kos_overlap$mgs_ko[, picrust2_subset$sample]), picrust2_subset$precision, method="spearman")
  
  all_precision <- c(all_precision, picrust2_subset$precision)
  all_depth <- c(all_depth, colSums(ko_metrics_out[[dataset]]$infiles$all_kos_overlap$mgs_ko[, picrust2_subset$sample]))
  
  depth_vs_precision[dataset, "rho"] <- picrust2_spearman$estimate
  depth_vs_precision[dataset, "p"] <- picrust2_spearman$p.value

}

depth_vs_precision_spearman <- cor.test(all_depth, all_precision, method="spearman")

