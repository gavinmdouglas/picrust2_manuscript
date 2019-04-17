
rm(list=ls())

setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

read_in_blueberry_18S_picrust2 <- function(infile, possible_funcs, funcs2keep, samples2keep) {
  in_table <- read.table(infile, header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
  in_table <- add_missing_funcs(in_table, possible_funcs)
  return(in_table[funcs2keep, samples2keep])
}

ec_18S <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_18S_counts.txt",
                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt", stringsAsFactors = FALSE)$V1
possible_18S_ecs <- colnames(ec_18S)
possible_ec <- possible_mgs_ecs[which(possible_mgs_ecs %in% possible_18S_ecs)]


pathabun_18S <- data.frame(t(read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_18S_pathway/path_abun_unstrat.tsv",
                                        row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

all_possible_18S_paths <- read.table("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt", stringsAsFactors = FALSE, sep="\t")$V1

pathabun_18S_no_miss <- data.frame(t(add_missing_funcs(t(pathabun_18S), all_possible_18S_paths)), check.names = FALSE) 

# Read in MGS.
ec_blueberry_mgs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/blueberry/humann2_ec_unstrat.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blueberry_mgs) <- gsub("^Bact", "", colnames(ec_blueberry_mgs))
ec_blueberry_mgs <- ec_blueberry_mgs[-which(! rownames(ec_blueberry_mgs) %in% possible_ec),]


# Read in MGS.
pathabun_blueberry_mgs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/blueberry/humann2_pathabun_unstrat.tsv",
                               header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blueberry_mgs) <- gsub("^Bact", "", colnames(pathabun_blueberry_mgs))
pathabun_blueberry_mgs <- pathabun_blueberry_mgs[-which(! rownames(pathabun_blueberry_mgs) %in% all_possible_18S_paths),]

# Read in 18S samples.
samples_18S <- colnames(read.table("EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE))

overlapping_samples <- colnames(ec_blueberry_mgs)[which(colnames(ec_blueberry_mgs) %in% samples_18S)]

# Subset to ECs and samples overlapping and calculate spearman and accuracy metrics based on EC numbers.
ec_blueberry_mgs_nomiss <- add_missing_funcs(ec_blueberry_mgs, possible_ec)
ec_blueberry_mgs_nomiss_subset <- ec_blueberry_mgs_nomiss[possible_ec, overlapping_samples]

ec_blueberry_18S_nsti2 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)
ec_blueberry_18S_nsti1.5 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)
ec_blueberry_18S_nsti1 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_1.0/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)
ec_blueberry_18S_nsti0.5 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)
ec_blueberry_18S_nsti0.25 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)
ec_blueberry_18S_nsti0.1 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_0.1/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)
ec_blueberry_18S_nsti0.05 <- read_in_blueberry_18S_picrust2("EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv", possible_18S_ecs, possible_ec, overlapping_samples)

blueberry_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_blueberry_mgs_nomiss_subset)
blueberry_ec_mgs_null_df_round <- round(blueberry_ec_mgs_null_df  - 0.00000001)

blueberry_ec_mgs_null_scc <- cor_all_cols(tab1 = blueberry_ec_mgs_null_df, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="Null", metric="spearman")
blueberry_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, blueberry_ec_mgs_null_df_round, category="Null")

blueberry_ec_mgs_nsti2_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti2, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
blueberry_ec_mgs_nsti2_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti2, category="NSTI=2")
rownames(blueberry_ec_mgs_nsti2_metrics) <- blueberry_ec_mgs_nsti2_metrics$sample
rownames(blueberry_ec_mgs_nsti2_scc) <- blueberry_ec_mgs_nsti2_scc$sample_names
     
blueberry_ec_mgs_nsti1.5_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti1.5, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
blueberry_ec_mgs_nsti1.5_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti1.5, category="NSTI=1.5")
rownames(blueberry_ec_mgs_nsti1.5_metrics) <- blueberry_ec_mgs_nsti1.5_metrics$sample
rownames(blueberry_ec_mgs_nsti1.5_scc) <- blueberry_ec_mgs_nsti1.5_scc$sample_names

blueberry_ec_mgs_nsti1_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti1, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
blueberry_ec_mgs_nsti1_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti1, category="NSTI=1")
rownames(blueberry_ec_mgs_nsti1_metrics) <- blueberry_ec_mgs_nsti1_metrics$sample
rownames(blueberry_ec_mgs_nsti1_scc) <- blueberry_ec_mgs_nsti1_scc$sample_names

blueberry_ec_mgs_nsti0.5_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti0.5, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
blueberry_ec_mgs_nsti0.5_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti0.5, category="NSTI=0.5")
rownames(blueberry_ec_mgs_nsti0.5_metrics) <- blueberry_ec_mgs_nsti0.5_metrics$sample
rownames(blueberry_ec_mgs_nsti0.5_scc) <- blueberry_ec_mgs_nsti0.5_scc$sample_names

blueberry_ec_mgs_nsti0.25_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti0.25, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
blueberry_ec_mgs_nsti0.25_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti0.25, category="NSTI=0.25")
rownames(blueberry_ec_mgs_nsti0.25_metrics) <- blueberry_ec_mgs_nsti0.25_metrics$sample
rownames(blueberry_ec_mgs_nsti0.25_scc) <- blueberry_ec_mgs_nsti0.25_scc$sample_names

blueberry_ec_mgs_nsti0.1_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti0.1, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
blueberry_ec_mgs_nsti0.1_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti0.1, category="NSTI=0.1")
rownames(blueberry_ec_mgs_nsti0.1_metrics) <- blueberry_ec_mgs_nsti0.1_metrics$sample
rownames(blueberry_ec_mgs_nsti0.1_scc) <- blueberry_ec_mgs_nsti0.1_scc$sample_names

blueberry_ec_mgs_nsti0.05_scc <- cor_all_cols(tab1 = ec_blueberry_18S_nsti0.05, tab2 = ec_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
blueberry_ec_mgs_nsti0.05_metrics <- calc_accuracy_metrics(ec_blueberry_mgs_nomiss_subset, ec_blueberry_18S_nsti0.05, category="NSTI=0.05")
rownames(blueberry_ec_mgs_nsti0.05_metrics) <- blueberry_ec_mgs_nsti0.05_metrics$sample
rownames(blueberry_ec_mgs_nsti0.05_scc) <- blueberry_ec_mgs_nsti0.05_scc$sample_names


# Subset to pathways and samples overlapping and calculate spearman and accuracy metrics based on pathway abundances.
pathabun_blueberry_mgs_nomiss <- add_missing_funcs(pathabun_blueberry_mgs, all_possible_18S_paths)
pathabun_blueberry_mgs_nomiss_subset <- pathabun_blueberry_mgs_nomiss[, overlapping_samples]

pathabun_blueberry_18S_nsti2 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_2.0/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)
pathabun_blueberry_18S_nsti1.5 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_1.5/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)
pathabun_blueberry_18S_nsti1 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_1.0/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)
pathabun_blueberry_18S_nsti0.5 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_0.5/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)
pathabun_blueberry_18S_nsti0.25 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_0.25/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)
pathabun_blueberry_18S_nsti0.1 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_0.1/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)
pathabun_blueberry_18S_nsti0.05 <- read_in_blueberry_18S_picrust2("pathways_out_nsti_0.05/path_abun_unstrat.tsv", all_possible_18S_paths, all_possible_18S_paths, overlapping_samples)

blueberry_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun_18S_no_miss, tab = pathabun_blueberry_mgs_nomiss_subset)
blueberry_pathabun_mgs_null_df_round <- round(blueberry_pathabun_mgs_null_df - 0.00000001)

blueberry_pathabun_mgs_null_scc <- cor_all_cols(tab1 = blueberry_pathabun_mgs_null_df, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="Null", metric="spearman")
blueberry_pathabun_mgs_null_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, blueberry_pathabun_mgs_null_df_round, category="Null")

blueberry_pathabun_mgs_nsti2_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti2, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
blueberry_pathabun_mgs_nsti2_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti2, category="NSTI=2")
rownames(blueberry_pathabun_mgs_nsti2_metrics) <- blueberry_pathabun_mgs_nsti2_metrics$sample
rownames(blueberry_pathabun_mgs_nsti2_scc) <- blueberry_pathabun_mgs_nsti2_scc$sample_names

blueberry_pathabun_mgs_nsti1.5_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti1.5, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
blueberry_pathabun_mgs_nsti1.5_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti1.5, category="NSTI=1.5")
rownames(blueberry_pathabun_mgs_nsti1.5_metrics) <- blueberry_pathabun_mgs_nsti1.5_metrics$sample
rownames(blueberry_pathabun_mgs_nsti1.5_scc) <- blueberry_pathabun_mgs_nsti1.5_scc$sample_names

blueberry_pathabun_mgs_nsti1_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti1, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
blueberry_pathabun_mgs_nsti1_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti1, category="NSTI=1")
rownames(blueberry_pathabun_mgs_nsti1_metrics) <- blueberry_pathabun_mgs_nsti1_metrics$sample
rownames(blueberry_pathabun_mgs_nsti1_scc) <- blueberry_pathabun_mgs_nsti1_scc$sample_names

blueberry_pathabun_mgs_nsti0.5_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti0.5, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
blueberry_pathabun_mgs_nsti0.5_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti0.5, category="NSTI=0.5")
rownames(blueberry_pathabun_mgs_nsti0.5_metrics) <- blueberry_pathabun_mgs_nsti0.5_metrics$sample
rownames(blueberry_pathabun_mgs_nsti0.5_scc) <- blueberry_pathabun_mgs_nsti0.5_scc$sample_names

blueberry_pathabun_mgs_nsti0.25_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti0.25, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
blueberry_pathabun_mgs_nsti0.25_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti0.25, category="NSTI=0.25")
rownames(blueberry_pathabun_mgs_nsti0.25_metrics) <- blueberry_pathabun_mgs_nsti0.25_metrics$sample
rownames(blueberry_pathabun_mgs_nsti0.25_scc) <- blueberry_pathabun_mgs_nsti0.25_scc$sample_names

blueberry_pathabun_mgs_nsti0.1_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti0.1, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
blueberry_pathabun_mgs_nsti0.1_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti0.1, category="NSTI=0.1")
rownames(blueberry_pathabun_mgs_nsti0.1_metrics) <- blueberry_pathabun_mgs_nsti0.1_metrics$sample
rownames(blueberry_pathabun_mgs_nsti0.1_scc) <- blueberry_pathabun_mgs_nsti0.1_scc$sample_names

blueberry_pathabun_mgs_nsti0.05_scc <- cor_all_cols(tab1 = pathabun_blueberry_18S_nsti0.05, tab2 = pathabun_blueberry_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
blueberry_pathabun_mgs_nsti0.05_metrics <- calc_accuracy_metrics(pathabun_blueberry_mgs_nomiss_subset, pathabun_blueberry_18S_nsti0.05, category="NSTI=0.05")
rownames(blueberry_pathabun_mgs_nsti0.05_metrics) <- blueberry_pathabun_mgs_nsti0.05_metrics$sample
rownames(blueberry_pathabun_mgs_nsti0.05_scc) <- blueberry_pathabun_mgs_nsti0.05_scc$sample_names


blueberry_ec_mgs_scc_df <- rbind(blueberry_ec_mgs_null_scc,
                            blueberry_ec_mgs_nsti2_scc,
                            blueberry_ec_mgs_nsti1.5_scc,
                            blueberry_ec_mgs_nsti1_scc,
                            blueberry_ec_mgs_nsti0.5_scc,
                            blueberry_ec_mgs_nsti0.25_scc,
                            blueberry_ec_mgs_nsti0.1_scc,
                            blueberry_ec_mgs_nsti0.05_scc)

blueberry_ec_mgs_acc_df <- rbind(blueberry_ec_mgs_null_metrics,
                            blueberry_ec_mgs_nsti2_metrics,
                            blueberry_ec_mgs_nsti1.5_metrics,
                            blueberry_ec_mgs_nsti1_metrics,
                            blueberry_ec_mgs_nsti0.5_metrics,
                            blueberry_ec_mgs_nsti0.25_metrics,
                            blueberry_ec_mgs_nsti0.1_metrics,
                            blueberry_ec_mgs_nsti0.05_metrics)

blueberry_pathabun_mgs_scc_df <- rbind(blueberry_pathabun_mgs_null_scc,
                                  blueberry_pathabun_mgs_nsti2_scc,
                                  blueberry_pathabun_mgs_nsti1.5_scc,
                                  blueberry_pathabun_mgs_nsti1_scc,
                                  blueberry_pathabun_mgs_nsti0.5_scc,
                                  blueberry_pathabun_mgs_nsti0.25_scc,
                                  blueberry_pathabun_mgs_nsti0.1_scc,
                                  blueberry_pathabun_mgs_nsti0.05_scc)

blueberry_pathabun_mgs_acc_df <- rbind(blueberry_pathabun_mgs_null_metrics,
                                  blueberry_pathabun_mgs_nsti2_metrics,
                                  blueberry_pathabun_mgs_nsti1.5_metrics,
                                  blueberry_pathabun_mgs_nsti1_metrics,
                                  blueberry_pathabun_mgs_nsti0.5_metrics,
                                  blueberry_pathabun_mgs_nsti0.25_metrics,
                                  blueberry_pathabun_mgs_nsti0.1_metrics,
                                  blueberry_pathabun_mgs_nsti0.05_metrics)

saveRDS(object = blueberry_ec_mgs_scc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_ec_scc_metrics.rds")

saveRDS(object = blueberry_ec_mgs_acc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_ec_acc_metrics.rds")

saveRDS(object = blueberry_pathabun_mgs_scc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_pathabun_scc_metrics.rds")

saveRDS(object = blueberry_pathabun_mgs_acc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/blueberry_pathabun_acc_metrics.rds")
