### Functions used for analyses in PICRUSt2 manuscript.

library(castor)
library(parallel)
library(stringi)

parse_rho_rds_and_calc_wilcoxon <- function(rho_rds, dataset_name, wilcox_cat2ignore, y_pos_start=0.97, dist_to_add=0.06) {
  
  rho_metrics <- readRDS(rho_rds)
  rho_metrics <- rho_metrics[with(rho_metrics, order(cat, sample_names)),]
  
  rho_metrics$Database <- "Other"
  rho_metrics$Database[grep("NSTI" , rho_metrics$cat)] <- "PICRUSt2"
  rho_metrics$Database[which(rho_metrics$cat=="Null")] <- "Null"
  
  rho_metrics$dataset <- dataset_name
  
  rho_metrics_wilcoxon <- paired_wilcoxon_vs_category(in_df = rho_metrics,
                                                  category_col = "cat",
                                                  metric_col = "metric",
                                                  y_pos_start=y_pos_start,
                                                  dist_to_add=dist_to_add,
                                                  categories2exclude = wilcox_cat2ignore)
  
  rho_metrics_wilcoxon$dataset <- dataset_name
  rho_metrics_wilcoxon$Database <- NA

  return(list(rho_metrics, rho_metrics_wilcoxon))
}

parse_acc_metrics_rds_and_calc_wilcoxon <- function(acc_rds, metric_col, dataset_name, wilcox_cat2ignore, y_pos_start=1, dist_to_add=0.06, levels2use=c("Null", "Tax4Fun2", "PanFP", "Piphillin", "PICRUSt1",
                                                                                                                                      "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5",
                                                                                                                                      "NSTI=0.25", "NSTI=0.1", "NSTI=0.05")) {
  
  acc_metrics <- readRDS(acc_rds)
  acc_metrics <- acc_metrics[, c("category", "sample", metric_col)]
  colnames(acc_metrics) <- c("cat", "sample_names", metric_col)
  
  acc_metrics$cat <- factor(acc_metrics$cat, levels=)
  
  acc_metrics <- acc_metrics[with(acc_metrics, order(cat, sample_names)),]
  
  acc_metrics$Database <- "Other"
  acc_metrics$Database[grep("NSTI" , acc_metrics$cat)] <- "PICRUSt2"
  acc_metrics$Database[which(acc_metrics$cat=="Null")] <- "Null"
  
  acc_metrics$dataset <- dataset_name
  
  acc_metrics_wilcoxon <- paired_wilcoxon_vs_category(in_df = acc_metrics,
                                                  category_col = "cat",
                                                  metric_col = metric_col,
                                                  y_pos_start=y_pos_start,
                                                  dist_to_add=dist_to_add,
                                                  categories2exclude = wilcox_cat2ignore)
  
  acc_metrics_wilcoxon$dataset <- dataset_name
  acc_metrics_wilcoxon$Database <- NA

  return(list(acc_metrics, acc_metrics_wilcoxon))
}


paired_wilcoxon_vs_category <- function(in_df, category_col, metric_col, focal_cat="NSTI=2", y_pos_start=0.85, categories2exclude=c(), dist_to_add=0.06) {
  
  categories <- levels(in_df[, category_col])
  
  categories <- categories[-which(categories == focal_cat)]
  
  for(c in categories2exclude) {
    categories <- categories[-which(categories == c)]
  }
  
  out_df <- data.frame(matrix(NA, nrow=length(categories), ncol=6))
  
  colnames(out_df) <- c("group1", "group2", "raw_p", "W", "p_symbol", "y.position")
  
  current_y_pos = y_pos_start
  
  i = 1
  for(category in categories) {
    out_df[i, c("group1", "group2")] <- c(focal_cat, category)
    
    wilcoxon_out <- wilcox.test(in_df[which(in_df[, category_col] == focal_cat), metric_col],
                                in_df[which(in_df[, category_col] == category), metric_col],
                                paired=TRUE)
    
    out_df[i, c("raw_p", "W")] <- c(wilcoxon_out$p.value, wilcoxon_out$statistic)
    
    if(wilcoxon_out$p.value >= 0.05) {
      out_df[i, "p_symbol"] <- "ns"
    } else if(wilcoxon_out$p.value >= 0.001) {
      out_df[i, "p_symbol"] <- "*"
    } else if(wilcoxon_out$p.value < 0.001) {
      out_df[i, "p_symbol"] <- "**"
    } else {
      stop(c("ERROR with this pvalue:", wilcoxon_out$p.value))
    }
    
    out_df[i, "y.position"] <- current_y_pos
    
    i = i + 1
    
    current_y_pos = current_y_pos + dist_to_add
  }
  return(out_df)
}

read_table_check_exists <- function(filename, ...) {
 
  if(!file.exists(filename)){
    print("The following file was not found:")
    print(filename)
  }
  
  return(read.table(filename, ...))
}

# Add in rows for missing functions with all set to 0.
add_missing_funcs <- function(in_df, all_funcs) {
  
  # Check if any functions in dataframe aren't in all_funcs.
  if(length(which(! rownames(in_df) %in% all_funcs)) > 0) {
    stop("error - rows in dataframe not found in all_funcs")
  }
  
  missing_funcs <- all_funcs[which(! all_funcs %in% rownames(in_df))]
  
  if(length(missing_funcs) == 0) {
    return(in_df)   
  }
  
  missing_df <- data.frame(matrix(0, nrow=length(missing_funcs), ncol=ncol(in_df)))
  
  rownames(missing_df) <- missing_funcs
  colnames(missing_df) <- colnames(in_df)
  
  return(rbind(in_df, missing_df))
}


read_in_pathway_predictions <- function(dataset) {
  
  # Note that the possible pathways output by PICRUSt2 is a subset of those output by HUMAnN2 already.
  possible_picrust2_pathways <- read.table("possible_path/picrust2_path.txt", header=F, stringsAsFactors = FALSE)$V1
  
  # Read in all pathway prediction tables (and MGS).
  picrust2_pathabun_nsti2_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti2.0.tsv", sep="")
  picrust2_pathabun_nsti1.5_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti1.5.tsv", sep="")
  picrust2_pathabun_nsti1_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti1.0.tsv", sep="")
  picrust2_pathabun_nsti0.5_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti0.5.tsv", sep="")
  picrust2_pathabun_nsti0.25_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti0.25.tsv", sep="")
  picrust2_pathabun_nsti0.1_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti0.1.tsv", sep="")
  picrust2_pathabun_nsti0.05_file <- paste("picrust2_out/", dataset, "_picrust2_path_nsti0.05.tsv", sep="")
  
  picrust2_pathabun_scrambled_mean_file <- paste("picrust2_scrambled_out/", dataset, "_pathway_mean.tsv", sep="")
  picrust2_pathabun_scrambled_median_file <- paste("picrust2_scrambled_out/", dataset, "_pathway_median.tsv", sep="")
  
  picrust2_pathabun_nsti2 <- read_table_check_exists(picrust2_pathabun_nsti2_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_nsti1.5 <- read_table_check_exists(picrust2_pathabun_nsti1.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_nsti1 <- read_table_check_exists(picrust2_pathabun_nsti1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_nsti0.5 <- read_table_check_exists(picrust2_pathabun_nsti0.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_nsti0.25 <- read_table_check_exists(picrust2_pathabun_nsti0.25_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_nsti0.1 <- read_table_check_exists(picrust2_pathabun_nsti0.1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_nsti0.05 <- read_table_check_exists(picrust2_pathabun_nsti0.05_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  picrust2_pathabun_scrambled_mean <- read_table_check_exists(picrust2_pathabun_scrambled_mean_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_pathabun_scrambled_median <- read_table_check_exists(picrust2_pathabun_scrambled_median_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  mgs_pathabun_file <- paste("../mgs_validation/", dataset, "/humann2_pathabun_unstrat.tsv", sep="")
  mgs_pathabun <- read_table_check_exists(mgs_pathabun_file, header=T, sep="\t", row.names=1)
  
  # Only keep MGS pathways that are in PICRUSt2 prokaryotic set.
  pathabun_rows2keep <- which(rownames(mgs_pathabun) %in% possible_picrust2_pathways)
  mgs_pathabun <- mgs_pathabun[pathabun_rows2keep,] 
  
  # Determine overlapping samples.
  start_num_samples <- length(colnames(picrust2_pathabun_nsti2))
  overlapping_samples <- colnames(picrust2_pathabun_nsti2)[which(colnames(picrust2_pathabun_nsti2) %in% colnames(mgs_pathabun))]
  
  print("original num samples:")
  print(start_num_samples)
  print("num final overlapping samples:")
  print(length(overlapping_samples))
  
  # Subset to overlapping samples only.
  picrust2_pathabun_nsti2 <- picrust2_pathabun_nsti2[, overlapping_samples]
  picrust2_pathabun_nsti1.5 <- picrust2_pathabun_nsti1.5[, overlapping_samples]
  picrust2_pathabun_nsti1 <- picrust2_pathabun_nsti1[, overlapping_samples]
  picrust2_pathabun_nsti0.5 <- picrust2_pathabun_nsti0.5[, overlapping_samples]
  picrust2_pathabun_nsti0.25 <- picrust2_pathabun_nsti0.25[, overlapping_samples]
  picrust2_pathabun_nsti0.1 <- picrust2_pathabun_nsti0.1[, overlapping_samples]
  picrust2_pathabun_nsti0.05 <- picrust2_pathabun_nsti0.05[, overlapping_samples]
  
  picrust2_pathabun_scrambled_mean <- picrust2_pathabun_scrambled_mean[, overlapping_samples]
  picrust2_pathabun_scrambled_median <- picrust2_pathabun_scrambled_median[, overlapping_samples]
  
  mgs_pathabun <- mgs_pathabun[, overlapping_samples]
  
  # Add in missing pathways to tables.
  picrust2_pathabun_nsti2_all <- add_missing_funcs(picrust2_pathabun_nsti2, possible_picrust2_pathways)
  picrust2_pathabun_nsti1.5_all <- add_missing_funcs(picrust2_pathabun_nsti1.5, possible_picrust2_pathways)
  picrust2_pathabun_nsti1_all <- add_missing_funcs(picrust2_pathabun_nsti1, possible_picrust2_pathways)
  picrust2_pathabun_nsti0.5_all <- add_missing_funcs(picrust2_pathabun_nsti0.5, possible_picrust2_pathways)
  picrust2_pathabun_nsti0.25_all <- add_missing_funcs(picrust2_pathabun_nsti0.25, possible_picrust2_pathways)
  picrust2_pathabun_nsti0.1_all <- add_missing_funcs(picrust2_pathabun_nsti0.1, possible_picrust2_pathways)
  picrust2_pathabun_nsti0.05_all <- add_missing_funcs(picrust2_pathabun_nsti0.05, possible_picrust2_pathways)
  
  picrust2_pathabun_scrambled_mean_all <- add_missing_funcs(picrust2_pathabun_scrambled_mean, possible_picrust2_pathways)
  picrust2_pathabun_scrambled_median_all <- add_missing_funcs(picrust2_pathabun_scrambled_median, possible_picrust2_pathways)

  mgs_pathabun_all <- add_missing_funcs(mgs_pathabun, possible_picrust2_pathways)

  # Create list of each set of predictions (with missing pathways added and not).
  nonzero_pathabuns <- list(picrust2_pathabun_nsti2=picrust2_pathabun_nsti2,
                            picrust2_pathabun_nsti1.5=picrust2_pathabun_nsti1.5,
                            picrust2_pathabun_nsti1=picrust2_pathabun_nsti1,
                            picrust2_pathabun_nsti0.5=picrust2_pathabun_nsti0.5,
                            picrust2_pathabun_nsti0.25=picrust2_pathabun_nsti0.25,
                            picrust2_pathabun_nsti0.1=picrust2_pathabun_nsti0.1,
                            picrust2_pathabun_nsti0.05=picrust2_pathabun_nsti0.05,
                            picrust2_pathabun_scrambled_mean=picrust2_pathabun_scrambled_mean,
                            picrust2_pathabun_scrambled_median=picrust2_pathabun_scrambled_median,
                            mgs_pathabun=mgs_pathabun)
  
  all_pathabuns <- list(picrust2_pathabun_nsti2=picrust2_pathabun_nsti2_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_nsti1.5=picrust2_pathabun_nsti1.5_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_nsti1=picrust2_pathabun_nsti1_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_nsti0.5=picrust2_pathabun_nsti0.5_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_nsti0.25=picrust2_pathabun_nsti0.25_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_nsti0.1=picrust2_pathabun_nsti0.1_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_nsti0.05=picrust2_pathabun_nsti0.05_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_scrambled_mean=picrust2_pathabun_scrambled_mean_all[possible_picrust2_pathways, ],
                        picrust2_pathabun_scrambled_median=picrust2_pathabun_scrambled_median_all[possible_picrust2_pathways, ],
                        mgs_pathabun=mgs_pathabun_all[possible_picrust2_pathways, ])
  
  return(list(nonzero_pathabun=nonzero_pathabuns, all_pathabun=all_pathabuns))
}


read_in_ec_predictions <- function(dataset) {
  
  # Read in possible ECs from each pipeline.
  possible_picrust2_ecs <- read_table_check_exists("possible_ec/picrust2_ec.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_mgs_ecs <- read_table_check_exists("possible_ec/humann2_ec.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_paprica_ecs <- read_table_check_exists("possible_ec/paprica_ec.txt", header=F, stringsAsFactors = FALSE)$V1

  # Identify ecs overlapping in all.
  overlapping_possible_ecs <- possible_picrust2_ecs[which(possible_picrust2_ecs %in% possible_mgs_ecs)]
  overlapping_possible_ecs <- overlapping_possible_ecs[which(overlapping_possible_ecs %in% possible_paprica_ecs)]
  
  # Read in all EC prediction tables (and MGS).
  picrust2_ec_nsti2_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti2.0.tsv", sep="")
  picrust2_ec_nsti1.5_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti1.5.tsv", sep="")
  picrust2_ec_nsti1_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti1.0.tsv", sep="")
  picrust2_ec_nsti0.5_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti0.5.tsv", sep="")
  picrust2_ec_nsti0.25_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti0.25.tsv", sep="")
  picrust2_ec_nsti0.1_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti0.1.tsv", sep="")
  picrust2_ec_nsti0.05_file <- paste("picrust2_out/", dataset, "_picrust2_ec_nsti0.05.tsv", sep="")

  picrust2_ec_scrambled_mean_file <- paste("picrust2_scrambled_out/", dataset, "_ec_mean.tsv", sep="")
  picrust2_ec_scrambled_median_file <- paste("picrust2_scrambled_out/", dataset, "_ec_median.tsv", sep="")
  
  picrust2_ec_nsti2 <- read_table_check_exists(picrust2_ec_nsti2_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_nsti1.5 <- read_table_check_exists(picrust2_ec_nsti1.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_nsti1 <- read_table_check_exists(picrust2_ec_nsti1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_nsti0.5 <- read_table_check_exists(picrust2_ec_nsti0.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_nsti0.25 <- read_table_check_exists(picrust2_ec_nsti0.25_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_nsti0.1 <- read_table_check_exists(picrust2_ec_nsti0.1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_nsti0.05 <- read_table_check_exists(picrust2_ec_nsti0.05_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  picrust2_ec_scrambled_mean <- read_table_check_exists(picrust2_ec_scrambled_mean_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ec_scrambled_median <- read_table_check_exists(picrust2_ec_scrambled_median_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  paprica_ec_file <- paste("paprica_out/", dataset, "_paprica_ec.csv", sep="")
  paprica_ec <- data.frame(t(read_table_check_exists(paprica_ec_file, header=T, sep=",", stringsAsFactors = FALSE,
                                        row.names=1, quote="", comment.char="", check.names=FALSE)))
  
  # The other tools don't output ECs that are missing in all samples, so remove these rows.
  paprica_ec <- paprica_ec[-which(rowSums(paprica_ec) == 0),]
  rownames(paprica_ec) <- gsub("^", "EC:", rownames(paprica_ec))
  
  mgs_ec_file <- paste("../mgs_validation/", dataset, "/humann2_ec_unstrat.tsv", sep="")

  mgs_ec <- read_table_check_exists(mgs_ec_file, header=T, sep="\t", row.names=1)
  
  rows2remove <- which(rownames(mgs_ec) %in% c("EC:UNMAPPED", "EC:UNGROUPED"))
  if(length(rows2remove) > 0) {
    mgs_ec <- mgs_ec[-rows2remove,]
  }
  
  # Determine overlapping samples.
  start_num_samples <- length(colnames(picrust2_ec_nsti2))
  overlapping_samples <- colnames(picrust2_ec_nsti2)[which(colnames(picrust2_ec_nsti2) %in% colnames(paprica_ec))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(mgs_ec))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(picrust2_ec_nsti2))]
  
  print("original num samples:")
  print(start_num_samples)
  print("num final overlapping samples:")
  print(length(overlapping_samples))
  
  # Subset to overlapping samples only.
  picrust2_ec_nsti2 <- picrust2_ec_nsti2[, overlapping_samples]
  picrust2_ec_nsti1.5 <- picrust2_ec_nsti1.5[, overlapping_samples]
  picrust2_ec_nsti1 <- picrust2_ec_nsti1[, overlapping_samples]
  picrust2_ec_nsti0.5 <- picrust2_ec_nsti0.5[, overlapping_samples]
  picrust2_ec_nsti0.25 <- picrust2_ec_nsti0.25[, overlapping_samples]
  picrust2_ec_nsti0.1 <- picrust2_ec_nsti0.1[, overlapping_samples]
  picrust2_ec_nsti0.05 <- picrust2_ec_nsti0.05[, overlapping_samples]
  
  picrust2_ec_scrambled_mean <- picrust2_ec_scrambled_mean[, overlapping_samples]
  picrust2_ec_scrambled_median <- picrust2_ec_scrambled_median[, overlapping_samples]
  
  paprica_ec <- paprica_ec[, overlapping_samples]
  mgs_ec <- mgs_ec[, overlapping_samples]
  
  # Add in missing ECs to PICRUSt2 tables.
  picrust2_ec_nsti2_all <- add_missing_funcs(picrust2_ec_nsti2, possible_picrust2_ecs)
  picrust2_ec_nsti1.5_all <- add_missing_funcs(picrust2_ec_nsti1.5, possible_picrust2_ecs)
  picrust2_ec_nsti1_all <- add_missing_funcs(picrust2_ec_nsti1, possible_picrust2_ecs)
  picrust2_ec_nsti0.5_all <- add_missing_funcs(picrust2_ec_nsti0.5, possible_picrust2_ecs)
  picrust2_ec_nsti0.25_all <- add_missing_funcs(picrust2_ec_nsti0.25, possible_picrust2_ecs)
  picrust2_ec_nsti0.1_all <- add_missing_funcs(picrust2_ec_nsti0.1, possible_picrust2_ecs)
  picrust2_ec_nsti0.05_all <- add_missing_funcs(picrust2_ec_nsti0.05, possible_picrust2_ecs)
  
  picrust2_ec_scrambled_mean_all <- add_missing_funcs(picrust2_ec_scrambled_mean, possible_picrust2_ecs)
  picrust2_ec_scrambled_median_all <- add_missing_funcs(picrust2_ec_scrambled_median, possible_picrust2_ecs)

  paprica_ec_all <- add_missing_funcs(paprica_ec, possible_paprica_ecs)
  mgs_ec_all <- add_missing_funcs(mgs_ec, possible_mgs_ecs)
  
  # Create list of each set of predictions (with missing ECs added and not).
  nonzero_ecs <- list(picrust2_ec_nsti2=picrust2_ec_nsti2,
                      picrust2_ec_nsti1.5=picrust2_ec_nsti1.5,
                      picrust2_ec_nsti1=picrust2_ec_nsti1,
                      picrust2_ec_nsti0.5=picrust2_ec_nsti0.5,
                      picrust2_ec_nsti0.25=picrust2_ec_nsti0.25,
                      picrust2_ec_nsti0.1=picrust2_ec_nsti0.1,
                      picrust2_ec_nsti0.05=picrust2_ec_nsti0.05,
                      picrust2_ec_scrambled_mean=picrust2_ec_scrambled_mean,
                      picrust2_ec_scrambled_median=picrust2_ec_scrambled_median,
                      paprica_ec=paprica_ec,
                      mgs_ec=mgs_ec)
  
  # Also define list with only KOs overlapping across all tools.
  all_ecs_overlap <- list(picrust2_ec_nsti2=picrust2_ec_nsti2_all[overlapping_possible_ecs,],
                          picrust2_ec_nsti1.5=picrust2_ec_nsti1.5_all[overlapping_possible_ecs,],
                          picrust2_ec_nsti1=picrust2_ec_nsti1_all[overlapping_possible_ecs,],
                          picrust2_ec_nsti0.5=picrust2_ec_nsti0.5_all[overlapping_possible_ecs,],
                          picrust2_ec_nsti0.25=picrust2_ec_nsti0.25_all[overlapping_possible_ecs,],
                          picrust2_ec_nsti0.1=picrust2_ec_nsti0.1_all[overlapping_possible_ecs,],
                          picrust2_ec_nsti0.05=picrust2_ec_nsti0.05_all[overlapping_possible_ecs,],
                          picrust2_ec_scrambled_mean=picrust2_ec_scrambled_mean_all[overlapping_possible_ecs,],
                          picrust2_ec_scrambled_median=picrust2_ec_scrambled_median_all[overlapping_possible_ecs,],
                          paprica_ec=paprica_ec_all[overlapping_possible_ecs,],
                          mgs_ec=mgs_ec_all[overlapping_possible_ecs,])
  
  # Return a list with these 3 lists as different indices.
  return(list(nonzero=nonzero_ecs, all_ecs_overlap=all_ecs_overlap))
}



read_in_ko_predictions <- function(dataset, verbose=TRUE) {
  
  if(verbose) {
    print("Reading in possible KOs")    
  }
  
  # Read in all possible ko output by each tool.
  possible_picrust2_kos <- read_table_check_exists("possible_ko/PICRUSt2_ko.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_mgs_kos <- read_table_check_exists("possible_ko/humann2_ko.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_picrust1_kos <- read_table_check_exists("possible_ko/PICRUSt1_ko.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_piphillin_kos <- read_table_check_exists("possible_ko/piphillin_ko.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_panfp_kos <- read_table_check_exists("possible_ko/PanFP_ko.txt", header=F, stringsAsFactors = FALSE)$V1
  possible_tax4fun2_kos <- read_table_check_exists("possible_ko/Tax4Fun2_ko.txt", header=F, stringsAsFactors = FALSE)$V1
  
  # Identify subset of KOs that could have been output by all approaches.

  overlapping_possible_kos <- possible_picrust2_kos[which(possible_picrust2_kos %in% possible_picrust1_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_mgs_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_piphillin_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_panfp_kos)]
  overlapping_possible_kos <- overlapping_possible_kos[which(overlapping_possible_kos %in% possible_tax4fun2_kos)]

  picrust2_ko_nsti2_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti2.0.tsv", sep="")
  picrust2_ko_nsti1.5_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti1.5.tsv", sep="")
  picrust2_ko_nsti1_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti1.0.tsv", sep="")
  picrust2_ko_nsti0.5_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.5.tsv", sep="")
  picrust2_ko_nsti0.25_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.25.tsv", sep="")
  picrust2_ko_nsti0.1_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.1.tsv", sep="")
  picrust2_ko_nsti0.05_file <- paste("picrust2_out/", dataset, "_picrust2_ko_nsti0.05.tsv", sep="")

  picrust2_ko_scrambled_mean_file <- paste("picrust2_scrambled_out/", dataset, "_ko_mean.tsv", sep="")
  picrust2_ko_scrambled_median_file <- paste("picrust2_scrambled_out/", dataset, "_ko_median.tsv", sep="")
  
  if(verbose) {
    print("Reading in unstratified KO tables")    
  }

  picrust2_ko_nsti2 <- read_table_check_exists(picrust2_ko_nsti2_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti1.5 <- read_table_check_exists(picrust2_ko_nsti1.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti1 <- read_table_check_exists(picrust2_ko_nsti1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.5 <- read_table_check_exists(picrust2_ko_nsti0.5_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.25 <- read_table_check_exists(picrust2_ko_nsti0.25_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.1 <- read_table_check_exists(picrust2_ko_nsti0.1_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_nsti0.05 <- read_table_check_exists(picrust2_ko_nsti0.05_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  picrust2_ko_scrambled_mean <- read_table_check_exists(picrust2_ko_scrambled_mean_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  picrust2_ko_scrambled_median <- read_table_check_exists(picrust2_ko_scrambled_median_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
  
  if(verbose) {
    print("Reading in alternative HSP tool tables")    
  }
  
  picrust1_ko_file <- paste("picrust1_out/", dataset, "_picrust1_ko.tsv", sep="")
  picrust1_ko <- read_table_check_exists(picrust1_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
  
  # The other tools don't output KOs that are missing in all samples, so remove these rows.
  picrust1_ko <- picrust1_ko[-which(rowSums(picrust1_ko) == 0),]
  
  panfp_ko_file <- paste("panfp_out/", dataset, "_panfp_ko.tsv", sep="")

  panfp_ko <- read_table_check_exists(panfp_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
  
  piphillin_ko_file <- paste("piphillin_out/", dataset, "_piphillin_ko.tsv", sep="")
  piphillin_ko <- read_table_check_exists(piphillin_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
  
  tax4fun2_ko_file <- paste("tax4fun2_out/", dataset, "_tax4fun2_ko.tsv", sep="")
  tax4fun2_ko <- read_table_check_exists(tax4fun2_ko_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

  if(verbose) {
    print("Reading in MGS table")    
  }
  
  mgs_ko_file <- paste("../mgs_validation/", dataset, "/humann2_ko_unstrat.tsv", sep="")
  mgs_ko <- read_table_check_exists(mgs_ko_file, header=T, sep="\t", row.names=1)
  
  rows2remove <- which(rownames(mgs_ko) %in% c("UNMAPPED", "UNGROUPED"))
  if(length(rows2remove) > 0) {
    mgs_ko <- mgs_ko[-rows2remove,]
  }
  
  # Determine overlapping samples.
  start_num_samples <- length(colnames(picrust2_ko_nsti2))
  overlapping_samples <- colnames(picrust1_ko)[which(colnames(picrust1_ko) %in% colnames(picrust2_ko_nsti2))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(panfp_ko))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(piphillin_ko))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(tax4fun2_ko))]
  overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(mgs_ko))]
  
  print("original num samples:")
  print(start_num_samples)
  print("num final overlapping samples:")
  print(length(overlapping_samples))
  
  # Subset to overlapping samples only.
  picrust2_ko_nsti2 <- picrust2_ko_nsti2[, overlapping_samples]
  picrust2_ko_nsti1.5 <- picrust2_ko_nsti1.5[, overlapping_samples]
  picrust2_ko_nsti1 <- picrust2_ko_nsti1[, overlapping_samples]
  picrust2_ko_nsti0.5 <- picrust2_ko_nsti0.5[, overlapping_samples]
  picrust2_ko_nsti0.25 <- picrust2_ko_nsti0.25[, overlapping_samples]
  picrust2_ko_nsti0.1 <- picrust2_ko_nsti0.1[, overlapping_samples]
  picrust2_ko_nsti0.05 <- picrust2_ko_nsti0.05[, overlapping_samples]
  picrust2_ko_scrambled_mean <- picrust2_ko_scrambled_mean[, overlapping_samples]
  picrust2_ko_scrambled_median <- picrust2_ko_scrambled_median[, overlapping_samples]
  picrust1_ko <- picrust1_ko[, overlapping_samples]
  panfp_ko <- panfp_ko[, overlapping_samples]
  piphillin_ko <- piphillin_ko[, overlapping_samples]
  tax4fun2_ko <- tax4fun2_ko[, overlapping_samples]
  mgs_ko <- mgs_ko[, overlapping_samples]
  
  # Add in missing KOs to each.
  picrust2_ko_nsti2_all <- add_missing_funcs(picrust2_ko_nsti2, possible_picrust2_kos)
  picrust2_ko_nsti1.5_all <- add_missing_funcs(picrust2_ko_nsti1.5, possible_picrust2_kos)
  picrust2_ko_nsti1_all <- add_missing_funcs(picrust2_ko_nsti1, possible_picrust2_kos)
  picrust2_ko_nsti0.5_all <- add_missing_funcs(picrust2_ko_nsti0.5, possible_picrust2_kos)
  picrust2_ko_nsti0.25_all <- add_missing_funcs(picrust2_ko_nsti0.25, possible_picrust2_kos)
  picrust2_ko_nsti0.1_all <- add_missing_funcs(picrust2_ko_nsti0.1, possible_picrust2_kos)
  picrust2_ko_nsti0.05_all <- add_missing_funcs(picrust2_ko_nsti0.05, possible_picrust2_kos)
  picrust2_ko_scrambled_mean_all <- add_missing_funcs(picrust2_ko_scrambled_mean, possible_picrust2_kos)
  picrust2_ko_scrambled_median_all <- add_missing_funcs(picrust2_ko_scrambled_median, possible_picrust2_kos)
  picrust1_ko_all <- add_missing_funcs(picrust1_ko, possible_picrust1_kos)
  panfp_ko_all <- add_missing_funcs(panfp_ko, possible_panfp_kos)
  piphillin_ko_all <- add_missing_funcs(piphillin_ko, possible_piphillin_kos)
  tax4fun2_ko_all <- add_missing_funcs(tax4fun2_ko, possible_tax4fun2_kos)
  mgs_ko_all <- add_missing_funcs(mgs_ko, possible_mgs_kos)
  
  # Create list of each set of predictions (with missing KOs added and not).
  nonzero_kos <- list(picrust2_ko_nsti2=picrust2_ko_nsti2,
                      picrust2_ko_nsti1.5=picrust2_ko_nsti1.5,
                      picrust2_ko_nsti1=picrust2_ko_nsti1,
                      picrust2_ko_nsti0.5=picrust2_ko_nsti0.5,
                      picrust2_ko_nsti0.25=picrust2_ko_nsti0.25,
                      picrust2_ko_nsti0.1=picrust2_ko_nsti0.1,
                      picrust2_ko_nsti0.05=picrust2_ko_nsti0.05,
                      picrust2_ko_scrambled_mean=picrust2_ko_scrambled_mean,
                      picrust2_ko_scrambled_median=picrust2_ko_scrambled_median,
                      picrust1_ko=picrust1_ko,
                      panfp_ko=panfp_ko,
                      piphillin_ko=piphillin_ko,
                      tax4fun2_ko=tax4fun2_ko,
                      mgs_ko=mgs_ko)
  
  # Also define list with only KOs overlapping across all tools.
  all_kos_overlap <- list(picrust2_ko_nsti2=picrust2_ko_nsti2_all[overlapping_possible_kos,],
                          picrust2_ko_nsti1.5=picrust2_ko_nsti1.5_all[overlapping_possible_kos,],
                          picrust2_ko_nsti1=picrust2_ko_nsti1_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.5=picrust2_ko_nsti0.5_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.25=picrust2_ko_nsti0.25_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.1=picrust2_ko_nsti0.1_all[overlapping_possible_kos,],
                          picrust2_ko_nsti0.05=picrust2_ko_nsti0.05_all[overlapping_possible_kos,],
                          picrust2_ko_scrambled_mean=picrust2_ko_scrambled_mean_all[overlapping_possible_kos,],
                          picrust2_ko_scrambled_median=picrust2_ko_scrambled_median_all[overlapping_possible_kos,],
                          picrust1_ko=picrust1_ko_all[overlapping_possible_kos,],
                          panfp_ko=panfp_ko_all[overlapping_possible_kos,],
                          piphillin_ko=piphillin_ko_all[overlapping_possible_kos,],
                          tax4fun2_ko=tax4fun2_ko_all[overlapping_possible_kos,],
                          mgs_ko=mgs_ko_all[overlapping_possible_kos,])
  
  # Return a list with these 3 lists as different indices.
  return(list(nonzero=nonzero_kos, all_kos_overlap=all_kos_overlap))
}

# Also calculate accuracy metrics for each category as well.
# Accuracy metrics based on presence/absence calls of functions.
# Df1 is assumed to be the gold standard.
calc_accuracy_metrics <- function(df1, df2, category) {
  
  # Subset only to columns and rows that overlap between both and convert to present (TRUE) and absent (FALSE)
  cols2keep <- colnames(df1)[which(colnames(df1) %in% colnames(df2))]
  rows2keep <- rownames(df1)[which(rownames(df1) %in% rownames(df2))]
  

  df1 <- df1[rows2keep, cols2keep, drop=FALSE] > 0
  df2 <- df2[rows2keep, cols2keep, drop=FALSE] > 0

  out_df <- data.frame(matrix(NA, nrow=length(cols2keep), ncol=12))
  colnames(out_df) <- c("category", "sample", "acc", "TP", "TN", "FP", "FN", "NPV", "precision", "recall", "fpr", "F1")
  
  row_i = 1
  for(sample in colnames(df1)) {
    
    total_func <- length(rows2keep)
    
    overall_acc <- sum(df1[,sample] == df2[,sample])/total_func
    
    num_true_pos <- length(which(which(df1[,sample]) %in% which(df2[,sample])))
    num_true_neg <- length(which(which(! df1[,sample]) %in% which(! df2[,sample])))
    
    num_false_pos <- length(which(which(! df1[,sample]) %in% which(df2[,sample])))
    num_false_neg <- length(which(which(df1[,sample]) %in% which(! df2[,sample])))
    
    npv <- num_true_neg/(num_true_neg + num_false_neg)
    precision <- num_true_pos/(num_true_pos + num_false_pos)
    recall <- num_true_pos/(num_true_pos + num_false_neg)
    fpr <- num_false_pos/(num_false_pos + num_true_neg)
    F1 <-  2 * ((precision * recall)/(precision + recall))
    
    out_df[row_i, ] <- c(NA, NA, overall_acc, num_true_pos, num_true_neg, num_false_pos, num_false_neg,
                         npv, precision, recall, fpr, F1)
    
    row_i = row_i + 1
  }

  out_df$category <- category  
  out_df$sample <- cols2keep

  return(out_df)
}

# Get basic summaries of trait database (# category, sparsity, mean, max, and sd of all values).
characterize_db <- function(in_tab) {
 
  in_tab <- as.matrix(in_tab)
  
  dim_out <- dim(in_tab)
  
  num_rows <- dim_out[1]
  num_cols <- dim_out[2]
  
  sparsity <- length(which(in_tab == 0)) / (num_rows *num_cols)
  
  tab_mean <- mean(in_tab)

  tab_max <- max(in_tab)
  
  tab_sd <- sd(in_tab)
  
  output <- c(num_rows, num_cols, sparsity, tab_mean, tab_max, tab_sd)
  names(output) <- c("nrows", "ncols", "sparsity", "mean", "max", "sd")
  
  return(output)
  
}

# Function to calculate trait depth using castor package for all traits in dataframe.
calc_trait_depth <- function(in_tree, in_tab, nproc=1) {
  
  # Reorder rows to be order of tip labels.
  in_tab_reordered <- in_tab[in_tree$tip.label,, drop=FALSE]

  # Recode any values > 0 to be 1
  in_tab_reordered[in_tab_reordered > 0] <- 1
  
  return(mclapply(in_tab_reordered, get_trait_depth, tree=in_tree, count_singletons=FALSE, weighted=TRUE))
  
}


rep_func_subset_cor <- function(func, rep, n) {
  
  cor_out <- c()
  
  for(i in 1:rep){
    
    cor_out <- c(cor_out, cor.test(ran_func_subset_sum(func, n), ran_func_subset_sum(func, n), method="spearman")$estimate)
    
  }
  
  return(cor_out)
  
}

# Simpler approach for getting null simply based on taking mean # of gene families across DB
# for subset overlapping in input table and return dataframe of these means with column names
# as samples.
generate_null_mean_db_funcs <- function(db, tab) {
  n_sample <- ncol(tab)
  tab <- tab[rownames(tab)[which(rownames(tab) %in% colnames(db))],]
  db_subset <- db[, rownames(tab), drop=FALSE]
  
  # Create empty null df
  null_df <- data.frame(matrix(NA,nrow=nrow(tab), ncol=ncol(tab)))
  colnames(null_df) <- colnames(tab)
  rownames(null_df) <- rownames(tab)
  
  db_mean_funcs <- colMeans(db_subset)
  
  for (i in 1:n_sample){
    null_df[,colnames(null_df)[i]] <- db_mean_funcs
  }
  
  return(null_df)
}


generate_random_table <- function(db, tab, ngenome=100) {
  n_sample <- ncol(tab)
  tab <- tab[rownames(tab)[which(rownames(tab) %in% colnames(db))],]
  db_subset <- db[,rownames(tab)]
  
  # Create simulated null df
  null_df <- data.frame(matrix(NA,nrow=nrow(tab), ncol=ncol(tab)))
  colnames(null_df) <- colnames(tab)
  rownames(null_df) <- rownames(tab)
  
  for (i in 1:n_sample){
    null_df[,colnames(null_df)[i]] <- ran_func_subset_sum(db_subset, ngenome)
  }
  
  return(null_df)
}


rand_sample_vs_func_table <- function(db, tab, return_df=TRUE, ngenome=100, metric="spearman") {
  n_sample <- ncol(tab)
  tab <- tab[rownames(tab)[which(rownames(tab) %in% colnames(db))],]
  db_subset <- db[,rownames(tab)]
  
  # Create simulated null df
  null_df <- data.frame(matrix(NA,nrow=nrow(tab), ncol=ncol(tab)))
  colnames(null_df) <- colnames(tab)
  rownames(null_df) <- rownames(tab)
  
  for (i in 1:n_sample){
    null_df[,colnames(null_df)[i]] <- ran_func_subset_sum(db_subset, ngenome)
  }
  
  
  metric_out <- get_df_concordance(null_df, tab, metric)
  
  if(return_df){
    return(data.frame(metric=metric_out, cat="null", metric_type=metric, sample_names=colnames(tab)))
  } else{
    return(metric_out)
  }
}

cor_all_cols <- function(tab1, tab2, return_df=TRUE, cat_string="actual", metric="spearman", subset_rows=TRUE) {
  
  cols2keep <- colnames(tab1)[which(colnames(tab1) %in% colnames(tab2))]
  
  if (subset_rows) {
    rows2keep <- rownames(tab1)[which(rownames(tab1) %in% rownames(tab2))]
    tab1_subset <- tab1[rows2keep, cols2keep, drop=FALSE]
    tab2_subset <- tab2[rows2keep, cols2keep, drop=FALSE]
  } else {
    tab1_subset <- tab1[ , cols2keep, drop=FALSE]
    tab2_subset <- tab2[ , cols2keep, drop=FALSE]
  }
  
  metric_out <- get_df_concordance(tab1_subset, tab2_subset, metric)
  
  if(return_df){
    return(data.frame(metric=metric_out, cat=cat_string, metric_type=metric, sample_names=cols2keep))
  } else{
    return(metric_out)
  }
}

# Define function that calculates pairwise concordance 
get_df_concordance <- function(df1, df2, concordance_method) {
  
  coeff <- c()
  
  for (column in colnames(df1)) {
    
    coeff <- c(coeff, compare_vecs(as.numeric(df1[, column]),
                                   as.numeric(df2[, column]),
                                   concordance_method))
  }
  
  return(coeff)
}


cosine_sim <- function(vec1, vec2) {
  
  ### Reads in 2 vectors and calculates the cosine similarity
  ### between them.
  dis_numerator <- sum(vec1 * vec2)
  dis_denominator <- sqrt(sum((vec1**2))) * sqrt(sum(vec2**2))
  
  return(dis_numerator / dis_denominator)
  
}

ran_func_subset_sum <- function(func, n, nbinom_size=10, nbinom_prob=0.7) {
  # Function to subsample table down to n genomes randomly, give each
  # genome an abundance sampled from a negative binomial dist, and
  # sum the abundances of all functions across these genomes.
  
  # Subset table to n genomes.
  func_subset <- func[sample(rownames(func), n),, drop=FALSE]
  
  # Get simulated abundances for all genomes.
  sim <- rnbinom(n, size=nbinom_size, prob=nbinom_prob)
  
  # Multiple the abundances of genes in each genome by the simulated abundance of the genome.
  func_subset <- func_subset * sim
  
  # Return the sum of each function in this simulated community.
  return(colSums(func_subset))
  
}

compare_vecs <- function(vec1, vec2, method) {
  ### Will compare 2 vectors and return correlation coefficient or similarity.
  ### "method" needs to be one of "spearman", "pearson", "jaccard", or "cosine".
  
  if(! method %in% c("spearman", "pearson", "jaccard", "cosine")) {
    stop(paste("method, ", method, " is not one of allowed inputs."))
  }
  
  if(method == "cosine") {
    
    return(cosine_sim(vec1, vec2)) 
    

  } else if(method == "spearman") {
    
    return(cor.test(vec1, vec2,method=method, exact=FALSE)$estimate)
    
  } else if(method == "pearson") {
    
    return(cor.test(vec1, vec2,method=method)$estimate)
    
  } else if(method == "jaccard") {
    
    return(cluster_similarity(convert_to_binary(vec1), convert_to_binary(vec2), 
                              similarity = "jaccard", method = "independence"))
    
  }
  
}

rep_func_subset_compare <- function(func, rep, n, method) {
  ### Function to calculate the similarity or correlation between
  ### random subsets from a functional abundance table (where each
  ### genome per subset is given abundance from negative binomial dist.).
  ### "rep" specifies the number of random subsets pairs to compare and
  ### "n" specifies the number of genomes in each subset.
  ### "method" needs to be one of the permitted methods passed to "compare_vecs".
  
  results <- c()
  
  for(i in 1:rep){
    
    results <- c(results, compare_vecs(ran_func_subset_sum(func, n),
                                       ran_func_subset_sum(func, n),
                                       method=method))
  }
  
  return(results)
  
}



plot_func_compare_reps <- function(func, method) {
  
  func_n1 <- rep_func_subset_compare(func, 10, 1, method)
  func_n2 <- rep_func_subset_compare(func, 10, 2, method)
  func_n3 <- rep_func_subset_compare(func, 10, 3, method)
  func_n4 <- rep_func_subset_compare(func, 10, 4, method)
  func_n5 <- rep_func_subset_compare(func, 10, 5, method)
  func_n6 <- rep_func_subset_compare(func, 10, 6, method)
  func_n7 <- rep_func_subset_compare(func, 10, 7, method)
  func_n8 <- rep_func_subset_compare(func, 10, 8, method)
  func_n9 <- rep_func_subset_compare(func, 10, 9, method)
  func_n10 <- rep_func_subset_compare(func, 10, 10, method)
  func_n25 <- rep_func_subset_compare(func, 10, 25, method)
  func_n50 <- rep_func_subset_compare(func, 10, 50, method)
  func_n75 <- rep_func_subset_compare(func, 10, 75, method)
  func_n100 <- rep_func_subset_compare(func, 10, 100, method)
  
  func_subsample_compare <- data.frame(
    "n=1"=func_n1,
    "n=2"=func_n2,
    "n=3"=func_n3,
    "n=4"=func_n4,
    "n=5"=func_n5,
    "n=6"=func_n6,
    "n=7"=func_n7,
    "n=8"=func_n8,
    "n=9"=func_n9,
    "n=10"=func_n10,
    "n=25"=func_n25,
    "n=50"=func_n50,
    "n=75"=func_n75,
    "n=100"=func_n100)
  
  boxplot(func_subsample_compare,
          names=sub("n.", "",colnames(func_subsample_compare)),
          ylab=method, xlab="Number of genomes subsampled",
          main="", ylim=c(0,1))
  
  return(func_subsample_compare)
}

clean_raw_humann2_out <- function(filename, outfile, col_str_to_remove, first_col, rm_descrip=FALSE, strat=FALSE, str2add="",
                                  replace_pf=FALSE, old_sample=NULL, new_sample=NULL) {
  
  intab <- read.table(filename, header=T, sep="\t", stringsAsFactors = FALSE, quote="", comment.char="", row.names=1, check.names = FALSE)
  colnames(intab) <- gsub(col_str_to_remove, "", colnames(intab))
  orig_cols <- colnames(intab)
  
  if(rm_descrip) {
    intab$first_col <- gsub(":.+$", "", rownames(intab))
  } else {
    intab$first_col <- sapply(rownames(intab), function(x) { strsplit(x, "\\|")[[1]][1] } )
  }
  
  if(strat) {
    intab$taxa <- sapply(rownames(intab), function(x) { strsplit(x, "\\|")[[1]][2] } )
    intab <- intab[, c("first_col", "taxa", orig_cols)]
  } else {
    intab <- intab[, c("first_col", orig_cols)]
  }
  
  if(length(str2add) > 0) {
    intab$first_col <- paste(str2add, intab$first_col, sep="")
  }
  
  if(replace_pf) {
    intab$first_col <- gsub("^PF", "pfam", intab$first_col)
  }
  
  # Extra things to replace in column names:
  colnames(intab) <- gsub("_R1_R2_cat", "", colnames(intab))
  colnames(intab) <- gsub("_kneaddata", "", colnames(intab))
  
  colnames(intab)[1] <- first_col
  
  # Swap sample names if mapping vectors given.
  if(! is.null(old_sample) & ! is.null(new_sample)) {
    
    samples_to_keep <- which(old_sample %in% colnames(intab))
    old_sample <- old_sample[samples_to_keep]
    new_sample <- new_sample[samples_to_keep]
    
    if(! strat) {
      intab <- intab[, c(first_col, old_sample)]
      colnames(intab) <- c(first_col, new_sample)
    } else {
      intab <- intab[, c(first_col, "taxa", old_sample)]
      colnames(intab) <- c(first_col, "taxa", new_sample)
    }
  }
  
  write.table(x = intab, file = outfile,
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  
}

# Function to try trimming end of genome name to see if that can be classified.
get_taxa_from_str <- function(str_in) {
  
  # Get number of spaces in input string and add 1.
  num_loop <- length(grep(" ", unlist(strsplit(str_in, "")))) + 1
  
  for(i in 1:num_loop) {
    taxa_out <- tryCatch(classification(c(str_in), db="ncbi", out_type="summary"), error=function(err) NA)
    # If got matching taxa then return.
    if(! is.na(taxa_out)) {
      return(taxa_out)
    }
    
    str_in <- gsub(" [^ ]*$", "", str_in)
    
  }
  
  return(NA)
  
}

# Function to split up single Taxon column output by QIIME2 into separate columns.
# Note that "unclassified" values, e.g. "s__", will be added even in cases where
# the classifier didn't reach lower levels.
add_tax_cols <- function(taxa_table) {
  
  # Make empty columns for taxa levels.
  taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for(taxa_level in taxa_levels) { 
    taxa_table[, taxa_level] <- NA 
  }
  
  unclassified_labels <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Loop over all ASVs and parse the taxa information into these individual columns.
  for(i in 1:nrow(taxa_table)) {
    
    # Get that ASV's full taxonomy.
    # This string will be altered at each loop of a taxonomic level below.
    seq_taxonomy <- taxa_table[i, "Taxon"]
    
    # Some taxonomy assignments wont go down to species so need to skip 
    # those by checking how many levels are classified (will do this based on length below).
    exp_length <- 7
    
    for(tax_level in rev(taxa_levels)) {
      
      # Go to next iteration if the lowest taxonomic level is higher than current tax_level.
      exp_length_diff <- exp_length - length(stri_split(seq_taxonomy, regex = "; ")[[1]])
      if(exp_length_diff != 0) {
        taxa_table[i, tax_level] <- paste(c(seq_taxonomy, 
                                            unclassified_labels[(exp_length - exp_length_diff + 1):exp_length]), 
                                          collapse="; ")
        exp_length <- exp_length - 1
        next
      }
      
      # Set taxa level for this sequence.
      taxa_table[i, tax_level] <- seq_taxonomy
      
      # Remove the last level from the taxonomy.
      seq_taxonomy <- sub("; [^;]*$", "", seq_taxonomy)
      
      exp_length <- exp_length - 1
    }
  }
  
  return(taxa_table)
}



# Function that reads in table outputted by add_tax_cols and BIOM Table and 
# calculates abundance tables at each taxonomic level. Will return named list.
abun_by_taxa_level <- function(taxa_table_w_levels, in_biom) {
  
  taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  summed_tables <- list()
  
  for(taxa_level in taxa_levels) {
    
    in_biom_tax <- in_biom
    
    in_biom_tax[, taxa_level] <- factor(taxa_table_w_levels[rownames(in_biom), taxa_level])
    
    summed_tab <- aggregate(x=in_biom[, 1:(ncol(in_biom_tax)-1)], by=list(in_biom_tax[, ncol(in_biom_tax)]), sum)
    
    # Set rownames to be first column and reverse this column before adding to list.
    rownames(summed_tab) <- summed_tab[, 1]
    summed_tables[[taxa_level]] <- summed_tab[, -1]
  }
  
  return(summed_tables)
}


### Function to compare ratios of functions contributed by defined ASVs between 2 sample groupings.
compare_func_ratios_by_asv_groups <- function(strat_table, ASVs_of_interest, sample_group1, sample_group2) {
  
  # Set first column name to be "func" (even if it is "pathway") and restrict to samples to those of interest.
  orig_colnames <- colnames(strat_table)
  colnames(strat_table) <- c("func", orig_colnames[2:length(orig_colnames)])
  strat_table <- strat_table[,c("func", "sequence", sample_group1, sample_group2)]
  
  # Remove rows that are all 0s.
  zero_rows <- which(rowSums(strat_table[, 3:ncol(strat_table)]) == 0)
  if(length(zero_rows) > 0) {
    strat_table <- strat_table[-zero_rows,]
  }
  
  # Split table by ASV groups.
  strat_table_asv_interest <- strat_table[which(strat_table$sequence %in% ASVs_of_interest),]
  strat_table_asv_other <- strat_table[which(! strat_table$sequence %in% ASVs_of_interest),]
  
  # Remove sequence column.
  strat_table_asv_interest <- strat_table_asv_interest[, -which(colnames(strat_table_asv_interest) == "sequence")]
  strat_table_asv_other <- strat_table_asv_other[, -which(colnames(strat_table_asv_other) == "sequence")]
  
  # Get summed function counts by function.
  strat_table_asv_interest_summed <- aggregate(. ~ func, data=strat_table_asv_interest, FUN=sum)
  strat_table_asv_other_summed <- aggregate(. ~ func, data=strat_table_asv_other, FUN=sum)
  
  # Set functions as rownames and remove column.
  rownames(strat_table_asv_interest_summed) <- strat_table_asv_interest_summed$func
  rownames(strat_table_asv_other_summed) <- strat_table_asv_other_summed$func
  strat_table_asv_interest_summed <- strat_table_asv_interest_summed[, -1]
  strat_table_asv_other_summed <- strat_table_asv_other_summed[, -1]
  
  # Add in pathways that might be missing from each table.
  all_funcs <- unique(c(rownames(strat_table_asv_interest_summed), rownames(strat_table_asv_other_summed)))
  strat_table_asv_interest_summed_NoMiss <- add_missing_funcs(strat_table_asv_interest_summed, all_funcs)
  strat_table_asv_other_summed_NoMiss <- add_missing_funcs(strat_table_asv_other_summed, all_funcs)  
  
  # Sort dataframes to be in same order and add pseudocount.
  strat_table_asv_interest_summed_NoMiss <- strat_table_asv_interest_summed_NoMiss + 1
  strat_table_asv_other_summed_NoMiss <- strat_table_asv_other_summed_NoMiss[rownames(strat_table_asv_interest_summed_NoMiss),] + 1
  
  # Get table of ratios between the ASV groups.
  out_ratio <- strat_table_asv_interest_summed_NoMiss/strat_table_asv_other_summed_NoMiss
  
  # Test for different ratios between each funcs between the sample groups.
  ratio_wilcoxon_raw_p <- c()
  ratio_wilcoxon_raw_w <- c()
  path_mean_diff <- c()
  
  for(func in rownames(out_ratio)) {
    
    func_wilcox <- wilcox.test(as.numeric(out_ratio[func, sample_group1]),
                               as.numeric(out_ratio[func, sample_group2]))
    
    ratio_wilcoxon_raw_p <- c(ratio_wilcoxon_raw_p, func_wilcox$p.value)
    ratio_wilcoxon_raw_w <- c(ratio_wilcoxon_raw_w, func_wilcox$statistic)
    
    path_mean_diff <- c(path_mean_diff, mean(as.numeric(out_ratio[func, sample_group1])) - mean(as.numeric(out_ratio[func, sample_group2])))
    
  }
  
  ratio_wilcoxon_raw <- data.frame(feature=rownames(out_ratio),
                                   wilcox_p=ratio_wilcoxon_raw_p,
                                   wilcox_w=ratio_wilcoxon_raw_w,
                                   mean_diff=path_mean_diff)
  
  return(list(interest=strat_table_asv_interest_summed_NoMiss,
              other=strat_table_asv_other_summed_NoMiss,
              ratio=out_ratio,
              wilcox=ratio_wilcoxon_raw))
  
}


# First compare numbers of genera contributing pathways (excluding "unclassified" in 16S biopsies and mgs stool).
num_contrib_genera <- function(pred_df, mgs_df) {
  
  # Remove unclassified rows.
  pred_df <- pred_df[-which(pred_df$genus == "unclassified"), ]
  mgs_df <- mgs_df[-which(mgs_df$genus == "unclassified"), ]
  
  all_path <- unique(c(pred_df$pathway, mgs_df$pathway))
  
  num_contrib_df <- data.frame(matrix(NA, nrow=length(all_path), ncol=4))
  colnames(num_contrib_df) <- c("picrust2_mean", "picrust2_sem", "mgs_mean", "mgs_sem")
  rownames(num_contrib_df) <- all_path
  
  for(pathway in all_path) {
    
    if(pathway %in% pred_df$pathway) {
      pred_df_subset <- pred_df[which(pred_df$pathway == pathway), -which(colnames(pred_df) %in% c("genus", "pathway"))]
      
      num_contrib_df[pathway, "picrust2_mean"] <- mean(colSums(pred_df_subset > 0))
      num_contrib_df[pathway, "picrust2_sem"] <- sd(colSums(pred_df_subset > 0)) / sqrt(ncol(pred_df_subset))
    }
    
    if(pathway %in% mgs_df$pathway) {
      mgs_df_subset <- mgs_df[which(mgs_df$pathway == pathway), ]
      
      num_contrib_df[pathway, "mgs_mean"] <- mean(colSums(mgs_df_subset > 0))
      num_contrib_df[pathway, "mgs_sem"] <- sd(colSums(mgs_df_subset > 0)) / sqrt(ncol(mgs_df_subset))
    }
    
  }
  
  return(num_contrib_df)
}

pathway_mean_sem <- function(in_df, pathway, col_str) {
  
  pathway_subset_abun <- in_df[which(in_df$pathway == pathway), ]
  
  pathway_subset_abun <- pathway_subset_abun[-which(pathway_subset_abun$genus == "unclassified"), ]
  
  rownames(pathway_subset_abun) <- pathway_subset_abun$genus
  
  pathway_subset_abun <- pathway_subset_abun[, -which(colnames(pathway_subset_abun) %in% c("genus", "pathway"))]
  
  pathway_subset_abun <- data.frame(sweep(pathway_subset_abun, 2, colSums(pathway_subset_abun), '/'), check.names = FALSE) * 100
  
  if(length(which(is.na(pathway_subset_abun))) > 0) {
    pathway_subset_abun[is.na(pathway_subset_abun)] <- 0
  }
  pathway_subset_abun_genus_mean <- rowMeans(pathway_subset_abun)
  
  pathway_subset_abun_genus_sd <- apply(pathway_subset_abun, 1, sd)
  
  pathway_subset_abun_genus_sem <- pathway_subset_abun_genus_sd / sqrt(ncol(pathway_subset_abun))
  
  pathway_subset_abun_genus_upper <- pathway_subset_abun_genus_mean + pathway_subset_abun_genus_sem * 2
  pathway_subset_abun_genus_lower <- pathway_subset_abun_genus_mean - pathway_subset_abun_genus_sem * 2
  
  if(length(which(pathway_subset_abun_genus_lower < 0)) > 0) {
    pathway_subset_abun_genus_lower[which(pathway_subset_abun_genus_lower < 0)] <- 0
  }
  
  out_df <- data.frame(mean=pathway_subset_abun_genus_mean,
                       sem=pathway_subset_abun_genus_sem,
                       lower=pathway_subset_abun_genus_lower,
                       upper=pathway_subset_abun_genus_upper)
  rownames(out_df) <- names(pathway_subset_abun_genus_mean) 
  colnames(out_df) <- paste(col_str, colnames(out_df), sep="_")
  
  return(out_df)
}

breakdown_mean_genera_contrib <- function(pred_df, mgs_df, pathway) {
  
  pred_summary <- pathway_mean_sem(pred_df, pathway, "picrust2")
  mgs_summary <- pathway_mean_sem(mgs_df, pathway, "mgs")
  
  merged_out <- merge(pred_summary, mgs_summary, by="row.names", all.x=TRUE)
  
  if(length(which(is.na(merged_out))) > 0 ) {
    merged_out[is.na(merged_out)] <- 0
  }
  return(merged_out)
}


run_wilcoxon_relab_tests <- function(table, metadata, dataset_name="unknown") {
  
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  print(paste("Comparison for", dataset_name, "limited to", as.character(length(group1_subset)),
              "vs", as.character(length(group2_subset)), "samples", sep=" "))
  
  table_relab <-  data.frame(sweep(table, 2, colSums(table), `/`)) * 100
  
  wilcox_out_df <- data.frame(matrix(NA, nrow=nrow(table_relab), ncol=3))
  rownames(wilcox_out_df) <- rownames(table)
  colnames(wilcox_out_df) <- c("mean_group1", "mean_group2", "wilcox_p")
  
  for(f in rownames(table_relab)) {
    wilcox_out_df[f, "mean_group1"] <- mean(as.numeric(table_relab[f, group1_subset]))
    wilcox_out_df[f, "mean_group2"] <- mean(as.numeric(table_relab[f, group2_subset]))
    wilcox_out_df[f, "wilcox_p"] <- wilcox.test(as.numeric(table_relab[f, group1_subset]), as.numeric(table_relab[f, group2_subset]))$p.value
  }
  
  wilcox_out_df$wilcox_BH <- p.adjust(wilcox_out_df$wilcox_p, "BH")
  
  return(wilcox_out_df)
}

run_aldex2_two_groups <- function(table, metadata, dataset_name="unknown",
                                  mc.samples.set=128, test.set="t", effect.set=TRUE,
                                  include.sample.summary.set=FALSE, denom.set="all",
                                  verbose.set=FALSE) {
  
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  print(paste("Comparison for", dataset_name, "limited to", as.character(length(group1_subset)),
              "vs", as.character(length(group2_subset)), "samples", sep=" "))
  
  aldex2_out <- aldex(reads = floor(table[, c(group1_subset, group2_subset)]),
                      conditions=c(rep("group1", length(group1_subset)),
                                   rep("group2", length(group2_subset))),
                      mc.samples.set=128,
                      test.set="t",
                      effect.set=TRUE,
                      include.sample.summary.set=FALSE,
                      denom.set="all",
                      verbose.set=FALSE)
  
  # Add in missing functions (which may have been excluded due to being 0 in all samples).
  if(nrow(aldex2_out) < nrow(table)) {
    missing_funcs <- rownames(table)[which(! rownames(table) %in% rownames(aldex2_out))]
    missing_df <- data.frame(matrix(NA, nrow=length(missing_funcs), ncol=ncol(aldex2_out)))
    colnames(missing_df) <- colnames(aldex2_out)
    rownames(missing_df) <- missing_funcs
    orig <- aldex2_out
    aldex2_out <- rbind(aldex2_out, missing_df)
  }
  
  return(aldex2_out)
}


run_default_deseq2_two_groups <- function(table, metadata, alpha.set, dataset_name="unknown") {
  
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  print(paste("Comparison for", dataset_name, "limited to", as.character(length(group1_subset)),
              "vs", as.character(length(group2_subset)), "samples", sep=" "))
  
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  # Create metadata df.
  metadata_tab <- data.frame(matrix(NA, nrow=(length(group1_subset) + length(group2_subset)), ncol=1))
  rownames(metadata_tab) <- c(group1_subset, group2_subset)
  colnames(metadata_tab) <- c("group")
  metadata_tab[group1_subset, "group"] <- "group1"
  metadata_tab[group2_subset, "group"] <- "group2"
  metadata_tab$group <- as.factor(metadata_tab$group)
  
  # Input count df needs to have columns in same order as metadata rows.
  table_subset <- floor(table[, rownames(metadata_tab)])
  
  dds <- DESeqDataSetFromMatrix(countData = table_subset,
                                colData = metadata_tab,
                                design = ~ group)
  default_deseq2 <- DESeq(dds)
  default_deseq2_results <- results(default_deseq2, alpha=alpha.set)
  
  return(default_deseq2_results)
}


run_GMmean_deseq2_two_groups <- function(table, metadata, alpha.set, dataset_name="unknown") {
  
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  print(paste("Comparison for", dataset_name, "limited to", as.character(length(group1_subset)),
              "vs", as.character(length(group2_subset)), "samples", sep=" "))
  
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  # Create metadata df.
  metadata_tab <- data.frame(matrix(NA, nrow=(length(group1_subset) + length(group2_subset)), ncol=1))
  rownames(metadata_tab) <- c(group1_subset, group2_subset)
  colnames(metadata_tab) <- c("group")
  metadata_tab[group1_subset, "group"] <- "group1"
  metadata_tab[group2_subset, "group"] <- "group2"
  metadata_tab$group <- as.factor(metadata_tab$group)
  
  # Input count df needs to have columns in same order as metadata rows.
  table_subset <- floor(table[, rownames(metadata_tab)])
  
  dds <- DESeqDataSetFromMatrix(countData = table_subset,
                                colData = metadata_tab,
                                design = ~ group)
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds = DESeq(dds, fitType="local")
  
  return(results(dds, alpha=alpha.set))
}


summarize_test_performance <- function(table_list, sig_column, reference_name, sig_cutoff=0.05) {
  
  for(n in names(table_list)) {
    table_list[[n]][which(is.na(table_list[[n]][, sig_column])), sig_column] <- 1 
  }
  
  test_categories <- names(table_list)[which(! names(table_list) %in% reference_name)]
  
  perf_df <- data.frame(matrix(NA, nrow=9, ncol=(length(names(table_list)) - 1)))
  colnames(perf_df) <- test_categories
  rownames(perf_df) <- c("false_pos", "true_pos", "false_neg", "true_neg", "precision",
                         "recall", "f1", "specificity", "fpr")
  
  reference_pos <- rownames(table_list[[reference_name]])[which(table_list[[reference_name]][, sig_column] < sig_cutoff)]
  reference_neg <- rownames(table_list[[reference_name]])[which(table_list[[reference_name]][, sig_column] >= sig_cutoff)]
  
  for(category in test_categories) {
    category_pos <- rownames(table_list[[reference_name]])[which(table_list[[category]][, sig_column] < sig_cutoff)]
    category_neg <- rownames(table_list[[reference_name]])[which(table_list[[category]][, sig_column] >= sig_cutoff)]
    
    category_false_pos <- category_pos[which(category_pos %in% reference_neg)]
    category_true_pos <- category_pos[which(category_pos %in% reference_pos)]
    category_false_neg <- category_neg[which(category_neg %in% reference_pos)]
    category_true_neg <- category_neg[which(category_neg %in% reference_neg)]
    
    perf_df["false_pos", category] <- length(category_false_pos)
    perf_df["true_pos", category] <- length(category_true_pos)
    perf_df["false_neg", category] <- length(category_false_neg)
    perf_df["true_neg", category] <- length(category_true_neg)
    perf_df["precision", category] <- length(category_true_pos) / (length(category_true_pos) + length(category_false_pos))
    perf_df["recall", category] <- length(category_true_pos) / (length(category_true_pos) + length(category_false_neg))
    perf_df["f1", category] <- 2 * ((perf_df["precision", category] * perf_df["recall", category]) / (perf_df["precision", category] + perf_df["recall", category]))
    perf_df["specificity", category] <- length(category_true_neg) / length(reference_neg)
    perf_df["fpr", category] <- length(category_false_pos) / length(reference_neg)
  }
  
  return(perf_df)
}


differential_prevalence <- function(table, metadata, dataset_name="unknown") {
  group1_subset <- metadata$group1[which(metadata$group1 %in% colnames(table))]
  group2_subset <- metadata$group2[which(metadata$group2 %in% colnames(table))]
  
  print(paste("Comparison for", dataset_name, "limited to", as.character(length(group1_subset)),
              "vs", as.character(length(group2_subset)), "samples", sep=" "))
  
  prevalence_out_df <- data.frame(matrix(NA, nrow=nrow(table), ncol=6))
  rownames(prevalence_out_df) <- rownames(table)
  colnames(prevalence_out_df) <- c("present_group1", "absent_group1", "present_group2", "absent_group2", "fishers_ratio", "fishers_p")
  
  for(f in rownames(table)) {
    
    present_group1 <- length(which(as.numeric(table[f, group1_subset]) > 0))
    absent_group1 <- length(which(as.numeric(table[f, group1_subset]) == 0))
    
    present_group2 <- length(which(as.numeric(table[f, group2_subset]) > 0))
    absent_group2 <- length(which(as.numeric(table[f, group2_subset]) == 0))
    
    fishers_exact_out <- fisher.test(matrix(c(present_group1, absent_group1, present_group2, absent_group2), nrow=2, ncol=2))
    
    prevalence_out_df[f, "present_group1"] <- present_group1
    prevalence_out_df[f, "present_group2"] <- present_group2
    prevalence_out_df[f, "absent_group1"] <- absent_group1
    prevalence_out_df[f, "absent_group2"] <- absent_group2
    prevalence_out_df[f, "fishers_ratio"] <- fishers_exact_out$estimate
    prevalence_out_df[f, "fishers_p"] <- fishers_exact_out$p.value
    
  }
  
  prevalence_out_df$fishers_BH <- p.adjust(prevalence_out_df$fishers_p, "BH")
  
  return(prevalence_out_df)
}
