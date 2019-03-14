### Comparing predicted functions based on 16S to MGS "gold standard".
### Comparisons made between KEGG orthologs predicted by each tool to MGS.
### RDS files saved for spearman correlations and accuracy metrics.

library(ggplot2)
library(reshape2)

setwd("/home/gavin/projects/picrust2_manuscript/data/16S_validation/hmmalign/")
source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")
read_in_ko_predictions_TEST <- function(dataset) {
  
  # Read in all possible ECs output by each tool.
  possible_picrust2_kos <- read.table("possible_KOs/PICRUSt2_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_mgs_kos <- read.table("possible_KOs/humann2_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_picrust1_kos <- read.table("possible_KOs/PICRUSt1_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_piphillin_kos <- read.table("possible_KOs/piphillin_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_panfp_kos <- read.table("possible_KOs/PanFP_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_tax4fun_kos <- read.table("possible_KOs/Tax4Fun_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
  
  # Identify subset of ECs that could have been output by all approaches.
  overlapping_possible_kos <- possible_picrust2_kos[which(possible_picrust2_kos %in% possible_picrust1_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_mgs_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_piphillin_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_panfp_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_tax4fun_kos)]
  
  # Read in all KO prediction tables (and MGS).
  picrust2_ko_nsti2_gg_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti2_GGonly.tsv", sep="")
  picrust2_ko_nsti2_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti2.tsv", sep="")
  picrust2_ko_nsti1.5_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti1.5.tsv", sep="")
  picrust2_ko_nsti1_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti1.tsv", sep="")
  picrust2_ko_nsti0.5_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.5.tsv", sep="")
  picrust2_ko_nsti0.25_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.25.tsv", sep="")
  picrust2_ko_nsti0.1_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.1.tsv", sep="")
  picrust2_ko_nsti0.05_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.05.tsv", sep="")
  
  picrust2_ko_nsti2_gg <- read.table(picrust2_ko_nsti2_gg_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti2 <- read.table(picrust2_ko_nsti2_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti1.5 <- read.table(picrust2_ko_nsti1.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti1 <- read.table(picrust2_ko_nsti1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.5 <- read.table(picrust2_ko_nsti0.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.25 <- read.table(picrust2_ko_nsti0.25_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.1 <- read.table(picrust2_ko_nsti0.1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.05 <- read.table(picrust2_ko_nsti0.05_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  picrust1_ko_file <- paste("picrust1_out/", dataset, "_picrust1_ko.tsv", sep="")
  picrust1_ko <- read.table(picrust1_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")
  
  # The other tools don't output KOs that are missing in all samples, so remove these rows.
  picrust1_ko <- picrust1_ko[-which(rowSums(picrust1_ko) == 0),]
  
  panfp_ko_file <- paste("panfp_out/", dataset, "_panfp_ko.tsv", sep="")
  panfp_ko <- read.table(panfp_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
  
  piphillin_ko_file <- paste("piphillin_out/", dataset, "_piphillin_ko.tsv", sep="")
  piphillin_ko <- read.table(piphillin_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
  
  tax4fun_ko_file <- paste("tax4fun_out/", dataset, "_tax4fun_ko.tsv", sep="")
  tax4fun_ko <- read.table(tax4fun_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
  
  mgs_ko_file <- paste("../../mgs_validation/", dataset, "/humann2_ko_unstrat.tsv", sep="")
  mgs_ko <- read.table(mgs_ko_file, header=T, sep="\t", row.names=1)
  
  rows2remove <- which(rownames(mgs_ko) %in% c("UNMAPPED", "UNGROUPED"))
  if(length(rows2remove) > 0) {
    mgs_ko <- mgs_ko[-rows2remove,]
  }
  
  # Determine overlapping samples.
  start_num_samples <- length(colnames(picrust2_ko_nsti2))
  overlapping_samples <- colnames(picrust1_ko)[which(colnames(picrust1_ko) %in% colnames(picrust2_ko_nsti2))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(tax4fun_ko))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(mgs_ko))]
  
  print("original num samples:")
  print(start_num_samples)
  print("num final overlapping samples:")
  print(length(overlapping_samples))
  
  # Subset to overlapping samples only.
  picrust2_ko_nsti2_gg <- picrust2_ko_nsti2_gg[, overlapping_samples]
  picrust2_ko_nsti2 <- picrust2_ko_nsti2[, overlapping_samples]
  picrust2_ko_nsti1.5 <- picrust2_ko_nsti1.5[, overlapping_samples]
  picrust2_ko_nsti1 <- picrust2_ko_nsti1[, overlapping_samples]
  picrust2_ko_nsti0.5 <- picrust2_ko_nsti0.5[, overlapping_samples]
  picrust2_ko_nsti0.25 <- picrust2_ko_nsti0.25[, overlapping_samples]
  picrust2_ko_nsti0.1 <- picrust2_ko_nsti0.1[, overlapping_samples]
  picrust2_ko_nsti0.05 <- picrust2_ko_nsti0.05[, overlapping_samples]
  picrust1_ko <- picrust1_ko[, overlapping_samples]
  panfp_ko <- panfp_ko[, overlapping_samples]
  piphillin_ko <- piphillin_ko[, overlapping_samples]
  tax4fun_ko <- tax4fun_ko[, overlapping_samples]
  mgs_ko <- mgs_ko[, overlapping_samples]
  
  # Add in missing KOs to each
  picrust2_ko_nsti2_gg_all <- add_missing_funcs(picrust2_ko_nsti2_gg, possible_picrust2_kos)
  picrust2_ko_nsti2_all <- add_missing_funcs(picrust2_ko_nsti2, possible_picrust2_kos)
  picrust2_ko_nsti1.5_all <- add_missing_funcs(picrust2_ko_nsti1.5, possible_picrust2_kos)
  picrust2_ko_nsti1_all <- add_missing_funcs(picrust2_ko_nsti1, possible_picrust2_kos)
  picrust2_ko_nsti0.5_all <- add_missing_funcs(picrust2_ko_nsti0.5, possible_picrust2_kos)
  picrust2_ko_nsti0.25_all <- add_missing_funcs(picrust2_ko_nsti0.25, possible_picrust2_kos)
  picrust2_ko_nsti0.1_all <- add_missing_funcs(picrust2_ko_nsti0.1, possible_picrust2_kos)
  picrust2_ko_nsti0.05_all <- add_missing_funcs(picrust2_ko_nsti0.05, possible_picrust2_kos)
  picrust1_ko_all <- add_missing_funcs(picrust1_ko, possible_picrust1_kos)
  panfp_ko_all <- add_missing_funcs(panfp_ko, possible_panfp_kos)
  piphillin_ko_all <- add_missing_funcs(piphillin_ko, possible_piphillin_kos)
  tax4fun_ko_all <- add_missing_funcs(tax4fun_ko, possible_tax4fun_kos)
  mgs_ko_all <- add_missing_funcs(mgs_ko, possible_mgs_kos)
  
  # Create list of each set of predictions (with missing KOs added and not).
  nonzero_kos <- list(picrust2_ko_nsti2_gg=picrust2_ko_nsti2_gg,
                      picrust2_ko_nsti2=picrust2_ko_nsti2,
                      picrust2_ko_nsti1.5=picrust2_ko_nsti1.5,
                      picrust2_ko_nsti1=picrust2_ko_nsti1,
                      picrust2_ko_nsti0.5=picrust2_ko_nsti0.5,
                      picrust2_ko_nsti0.25=picrust2_ko_nsti0.25,
                      picrust2_ko_nsti0.1=picrust2_ko_nsti0.1,
                      picrust2_ko_nsti0.05=picrust2_ko_nsti0.05,
                      picrust1_ko=picrust1_ko,
                      panfp_ko=panfp_ko,
                      piphillin_ko=piphillin_ko,
                      tax4fun_ko=tax4fun_ko,
                      mgs_ko=mgs_ko)
  
  all_kos <-     list(picrust2_ko_nsti2_gg=picrust2_ko_nsti2_gg_all,
                      picrust2_ko_nsti2=picrust2_ko_nsti2_all,
                      picrust2_ko_nsti1.5=picrust2_ko_nsti1.5_all,
                      picrust2_ko_nsti1=picrust2_ko_nsti1_all,
                      picrust2_ko_nsti0.5=picrust2_ko_nsti0.5_all,
                      picrust2_ko_nsti0.25=picrust2_ko_nsti0.25_all,
                      picrust2_ko_nsti0.1=picrust2_ko_nsti0.1_all,
                      picrust2_ko_nsti0.05=picrust2_ko_nsti0.05_all,
                      picrust1_ko=picrust1_ko_all,
                      panfp_ko=panfp_ko_all,
                      piphillin_ko=piphillin_ko_all,
                      tax4fun_ko=tax4fun_ko_all,
                      mgs_ko=mgs_ko_all)
  
  
  # Also define list with only KOs overlapping across all tools.
  all_kos_overlap <- list(picrust2_ko_nsti2_gg=picrust2_ko_nsti2_gg_all[overlapping_possible_kos,],
                          picrust2_ko_nsti2=picrust2_ko_nsti2_all[overlapping_possible_kos,],
                          picrust2_ko_nsti1.5=picrust2_ko_nsti1.5_all[overlapping_possible_kos,],
                          picrust2_ko_nsti1=picrust2_ko_nsti1_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.5=picrust2_ko_nsti0.5_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.25=picrust2_ko_nsti0.25_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.1=picrust2_ko_nsti0.1_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.05=picrust2_ko_nsti0.05_all[overlapping_possible_kos,],
                          picrust1_ko=picrust1_ko_all[overlapping_possible_kos,],
                          panfp_ko=panfp_ko_all[overlapping_possible_kos,],
                          piphillin_ko=piphillin_ko_all[overlapping_possible_kos,],
                          tax4fun_ko=tax4fun_ko_all[overlapping_possible_kos,],
                          mgs_ko=mgs_ko_all[overlapping_possible_kos,])
  
  # Return a list with these 3 lists as different indices.
  return(list(nonzero=nonzero_kos, all_kos=all_kos, all_kos_overlap=all_kos_overlap))
}

### Read in database files (used for calculating null distribution).
ko <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ko.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)

# Read in all tables, restrict to overlapping samples only, and get subsets with all possible KOs filled in and with
# only KOs overlapping between all tools. Focus analyses on KOs overlapping between all tools, but outputted others as well
# for sanity checks.
mammal_infiles <- read_in_ko_predictions_TEST("mammal")


# Generate random tables based on subsampling database to calculate null distributions.
mammal_ko_mgs_null_df <- generate_null_mean_db_funcs(db = ko, tab = mammal_infiles$all_kos_overlap$mgs_ko)


mammal_ko_mgs_null_df_round <- round(mammal_ko_mgs_null_df - 0.00000001)



mammal_ko_mgs_null <- cor_all_cols(tab1 = mammal_ko_mgs_null_df, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="Null", metric="spearman")
mammal_ko_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2 (GG)", metric="spearman")
mammal_ko_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti2, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2", metric="spearman")
mammal_ko_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti1.5, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1.5", metric="spearman")
mammal_ko_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti1, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1", metric="spearman")
mammal_ko_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.5, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.5", metric="spearman")
mammal_ko_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.25, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.25", metric="spearman")
mammal_ko_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.1, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.1", metric="spearman")
mammal_ko_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.05, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.05", metric="spearman")
mammal_ko_picrust1_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$picrust1_ko, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="PICRUSt1", metric="spearman")
mammal_ko_panfp_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$panfp_ko, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="PanFP", metric="spearman")
mammal_ko_piphillin_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$piphillin_ko, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="Piphillin", metric="spearman")
mammal_ko_tax4fun_vs_mgs <- cor_all_cols(tab1 = mammal_infiles$all_kos_overlap$tax4fun_ko, tab2 = mammal_infiles$all_kos_overlap$mgs_ko, cat_string="Tax4Fun", metric="spearman")


mammal_ko_spearman_df <- rbind(mammal_ko_mgs_null, mammal_ko_tax4fun_vs_mgs, mammal_ko_panfp_vs_mgs, mammal_ko_piphillin_vs_mgs,
                            mammal_ko_picrust1_vs_mgs, mammal_ko_picrust2_nsti2_gg_vs_mgs, mammal_ko_picrust2_nsti2_vs_mgs,
                            mammal_ko_picrust2_nsti1.5_vs_mgs, mammal_ko_picrust2_nsti1_vs_mgs, mammal_ko_picrust2_nsti0.5_vs_mgs,
                            mammal_ko_picrust2_nsti0.25_vs_mgs, mammal_ko_picrust2_nsti0.1_vs_mgs, mammal_ko_picrust2_nsti0.05_vs_mgs)


mammal_picrust2_ko_nsti2_gg_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, category="NSTI=2 (GG)")
mammal_picrust2_ko_nsti2_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti2, category="NSTI=2")
mammal_picrust2_ko_nsti1.5_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti1.5, category="NSTI=1.5")
mammal_picrust2_ko_nsti1_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti1, category="NSTI=1")
mammal_picrust2_ko_nsti0.5_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.5, category="NSTI=0.5")
mammal_picrust2_ko_nsti0.25_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.25, category="NSTI=0.25")
mammal_picrust2_ko_nsti0.1_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.1, category="NSTI=0.1")
mammal_picrust2_ko_nsti0.05_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust2_ko_nsti0.05, category="NSTI=0.05")
mammal_picrust1_ko_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$picrust1_ko, category="PICRUSt1")
mammal_piphillin_ko_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$piphillin, category="Piphillin")
mammal_panfp_ko_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$panfp, category="PanFP")
mammal_tax4fun_ko_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_infiles$all_kos_overlap$tax4fun, category="Tax4Fun")
mammal_null_ko_metrics <- calc_accuracy_metrics(mammal_infiles$all_kos_overlap$mgs_ko, mammal_ko_mgs_null_df_round, category="Null")



mammal_acc <- rbind(mammal_null_ko_metrics,
                 mammal_tax4fun_ko_metrics,
                 mammal_panfp_ko_metrics,
                 mammal_piphillin_ko_metrics,
                 mammal_picrust1_ko_metrics,
                 mammal_picrust2_ko_nsti2_gg_metrics,
                 mammal_picrust2_ko_nsti2_metrics,
                 mammal_picrust2_ko_nsti1.5_metrics, 
                 mammal_picrust2_ko_nsti1_metrics,
                 mammal_picrust2_ko_nsti0.5_metrics,
                 mammal_picrust2_ko_nsti0.25_metrics,
                 mammal_picrust2_ko_nsti0.1_metrics,
                 mammal_picrust2_ko_nsti0.05_metrics)


mammal_ko_spearman <- mammal_ko_spearman_df
colnames(mammal_ko_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
mammal_ko_spearman <- mammal_ko_spearman[with(mammal_ko_spearman, order(Category, Sample)),]
mammal_ko_spearman$Database <- "Other"
mammal_ko_spearman$Database[grep("NSTI" , mammal_ko_spearman$Category)] <- "PICRUSt2"
mammal_ko_spearman$Database[which(mammal_ko_spearman$Category=="Null")] <- "Null"
mammal_ko_spearman$Category <- factor(mammal_ko_spearman$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_ko_spearman_melt <- melt(mammal_ko_spearman)

mammal_spearman_boxplots <- ggplot(mammal_ko_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Mammal")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))




mammal_acc_metrics <- mammal_acc

mammal_acc_metrics_subset <- mammal_acc_metrics[,c("sample", "precision", "recall", "category")]

colnames(mammal_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")

mammal_acc_metrics_subset$Database <- "Other"

mammal_acc_metrics_subset$Database[grep("NSTI" , mammal_acc_metrics_subset$Category)] <- "PICRUSt2"

mammal_acc_metrics_subset$Database[which(mammal_acc_metrics_subset$Category=="Null")] <- "Null"

mammal_acc_metrics_subset$Category <- factor(mammal_acc_metrics_subset$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

mammal_acc_metrics_subset_melt <- melt(mammal_acc_metrics_subset)

mammal_acc_metrics_boxplots <- ggplot(mammal_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.25, 1.0)) + ylab(c("")) + ggtitle("Mammal") + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))
