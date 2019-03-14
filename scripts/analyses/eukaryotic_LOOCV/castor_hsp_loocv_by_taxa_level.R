#!/usr/bin/env Rscript

# Read in quietly to avoid outputting "Loading required package: Rcpp" to stderr.
library(castor, quietly = TRUE)
library(ape)
# Load parallel package to run over multiple cores.
library(parallel)

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

Args <- commandArgs(TRUE)

# Read in command-line arguments.
full_tree <- read.tree(file=Args[1])
trait_values <- read.delim(Args[2], check.names=FALSE, row.names=1)
taxa_levels <- read.table(Args[3], sep="\t", header=TRUE, stringsAsFactors=FALSE)
tax_level <- Args[4]
expect_outfile <- Args[5]
predict_outfile <- Args[6]
metrics_outfile <- Args[7]
num_cores <- as.numeric(Args[8])

trait_values_ordered <- trait_values[full_tree$tip.label, , drop=FALSE]
trait_values_ordered <- as.data.frame(trait_values_ordered, drop=FALSE) + 1

num_tip <- nrow(trait_values_ordered)

# Sort taxa levels df to be in same order as others.
rownames(taxa_levels) <- taxa_levels$assembly
taxa_levels <- taxa_levels[full_tree$tip.label, , drop=FALSE]

# Get unique ids (excluding NAs) for this level.
unique_tax <- unique(na.omit(taxa_levels[,tax_level]))

# Prep predicted out df.
predict_out <- data.frame(matrix(NA, nrow=0, ncol=(ncol(trait_values_ordered) + 2)))
colnames(predict_out) <- c(colnames(trait_values_ordered), "taxon", "nsti")

exp_out <- data.frame(matrix(NA, nrow=0, ncol=(ncol(trait_values_ordered) + 1)))
colnames(exp_out) <- colnames(trait_values_ordered)

# Prep metric outfile.
metrics_out <- data.frame(matrix(NA, nrow=length(unique_tax), ncol=4))
rownames(metrics_out) <- unique_tax
colnames(metrics_out) <- c("rho", "recall", "precision", "nsti")

# Loop through each unique id and identify rows corresponding to this id.
for(tax in unique_tax) {
  
  matching_rows <- which(taxa_levels[,tax_level] == tax)
  
  rows2look <- matching_rows
  
  # Did not run the below commented out code.
  # Subset to 10 rows randomly if more than 10 match this taxa label (to speed up).
  #if(length(rows2look) > 10) {
  #  rows2look <- sample(rows2look, 10)
  #}

  # Loop through each matching row and get predictions for that tip while leaving the rest out.
  for(row_i in rows2look) {
    
    row_i_orig_index <- which(matching_rows == row_i)
    
    rows2remove <- matching_rows[-row_i_orig_index]
    
    trait_values_ordered_subset <- trait_values_ordered
    full_tree_subset <- full_tree
    
    exp_out <- rbind(exp_out, trait_values_ordered_subset[row_i,])
    
    trait_values_ordered_subset[row_i,] <- NA
    
    if(length(rows2remove) > 0) {
      # Remove these tips from trait table and tree.
      trait_values_ordered_subset <- trait_values_ordered_subset[-rows2remove,]
      full_tree_subset <- drop.tip(full_tree_subset, rows2remove)
    }
    
    hsp_out_models <- mclapply(trait_values_ordered_subset,
                               hsp_max_parsimony,
                               tree = full_tree_subset,
                               transition_costs = "proportional",
                               weight_by_scenarios = TRUE,
                               edge_exponent=0.5,
                               mc.cores = num_cores)
    
    rep_out <- unlist(mclapply(hsp_out_models, function(x) { max.col(x$likelihoods[row_i, , drop=FALSE]) - 1 },
                               mc.cores = num_cores))
    
    rep_df <- data.frame(t(matrix(rep_out)))
    
    colnames(rep_df) <- names(rep_out)
    
    rep_df$taxon <- tax
    rep_df$nsti <- find_nearest_tips(full_tree_subset, target_tips=full_tree_subset$tip.label[-row_i])$nearest_distance_per_tip[row_i]
    predict_out <- rbind(predict_out, rep_df)
    
  }
}

metrics_out_df <- predict_out[, c("taxon", "nsti"), drop=FALSE]
rownames(metrics_out_df) <- rownames(exp_out)
metrics_out_df$spearman <- NA

for(row_i in 1:nrow(exp_out)) {
  metrics_out_df[row_i, "spearman"] <- compare_vecs(as.numeric(exp_out[row_i,]), as.numeric(predict_out[row_i, 1:(ncol(predict_out)-2)]), "spearman")
  
  exp_binary <- as.numeric(exp_out[row_i,]) > 0
  
  predict_binary <- as.numeric(predict_out[row_i, 1:(ncol(predict_out)-2)]) > 0
  
  TP = which(which(exp_binary) %in% which(predict_binary))
  FP = which(which(! exp_binary) %in% which(predict_binary))
  FN = which(which(exp_binary) %in% which(! predict_binary))
  
  metrics_out_df[row_i, "precision"] <- TP/(TP + FP) 

  metrics_out_df[row_i, "recall"] <- TP/(TP + FN)

}

metrics_out_df_mean <- aggregate(.~taxon, data=metrics_out_df, mean)

rownames(predict_out) <- rownames(exp_out)

write.table(x = metrics_out_df_mean,
            file = metrics_outfile,
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE,
            sep="\t")

write.table(x = predict_out,
            file = predict_outfile,
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep="\t")

write.table(x = exp_out,
            file = expect_outfile,
            quote=FALSE,
            row.names=TRUE,
            col.names=TRUE,
            sep="\t")
