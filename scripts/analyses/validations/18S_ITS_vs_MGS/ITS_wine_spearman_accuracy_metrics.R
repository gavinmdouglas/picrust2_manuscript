
rm(list=ls())

setwd("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

read_in_wine_ITS_picrust2 <- function(infile, possible_ecs, ecs2keep, samples2keep) {
  in_table <- read.table(infile, header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
  colnames(in_table) <- gsub("-", ".", colnames(in_table))
  in_table <- add_missing_funcs(in_table, possible_ecs)
  return(in_table[ecs2keep, samples2keep])
}

ec_ITS <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_ITS_counts.txt",
                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt", stringsAsFactors = FALSE)$V1
possible_ITS_ecs <- colnames(ec_ITS)
overlapping_ec <- possible_mgs_ecs[which(possible_mgs_ecs %in% possible_ITS_ecs)]


pathabun_ITS <- data.frame(t(read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_ITS_pathway/path_abun_unstrat.tsv",
                                        row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)
# Note that the "all" set of pathways includes all possibly predicted pathways - since the genome database is so small not all
# possible pathways were identified in individual genomes so these needed to be filled in as 0s.
all_possible_ITS_paths <- read.table("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt", stringsAsFactors = FALSE, sep="\t")$V1
pathabun_ITS_no_miss <- data.frame(t(add_missing_funcs(t(pathabun_ITS), all_possible_ITS_paths)), check.names = FALSE)

wine_id_map <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/wine_id_mapping.txt", header=T, sep="\t", stringsAsFactors = FALSE)
rownames(wine_id_map) <- wine_id_map$MGS_run

# Read in MGS.
ec_wine_mgs <- read.table("../../../mgs/humann2_final_out/humann2_level4ec_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")
rownames(ec_wine_mgs) <- gsub("^", "EC:", rownames(ec_wine_mgs))
colnames(ec_wine_mgs) <- gsub("_kneaddata_Abundance.RPKs", "", colnames(ec_wine_mgs))
ec_wine_mgs <- ec_wine_mgs[-which(rownames(ec_wine_mgs) %in% c("EC:UNMAPPED", "EC:UNGROUPED")),]
colnames(ec_wine_mgs) <- wine_id_map[colnames(ec_wine_mgs), "X16S_name"]

# Read in MGS.
pathabun_wine_mgs <- read.table("../../../mgs/humann2_final_out/humann2_pathabundance_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")
rownames(pathabun_wine_mgs) <- gsub(":.*$", "", rownames(pathabun_wine_mgs))
colnames(pathabun_wine_mgs) <- gsub("_kneaddata_Abundance", "", colnames(pathabun_wine_mgs))
pathabun_wine_mgs <- pathabun_wine_mgs[-which(! rownames(pathabun_wine_mgs) %in% all_possible_ITS_paths),]
colnames(pathabun_wine_mgs) <- wine_id_map[colnames(pathabun_wine_mgs), "X16S_name"]

# Read in ITS samples.
samples_ITS <- colnames(read.table("EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE))
samples_ITS <- gsub("-", ".", samples_ITS)
overlapping_samples <- colnames(ec_wine_mgs)[which(colnames(ec_wine_mgs) %in% samples_ITS)]

# Remove mock samples:
overlapping_samples <- overlapping_samples[-grep("control", overlapping_samples)]

# Restrict D1 samples to replicates A and D3 samples to replicates B (so that at least they aren't the same biological sample).
overlapping_samples <- overlapping_samples[-grep("D1.B", overlapping_samples)]
overlapping_samples <- overlapping_samples[-grep("D3.A", overlapping_samples)]

# Subset to ECs and samples overlapping and calculate spearman and accuracy metrics based on EC numbers.
ec_wine_mgs_nomiss <- add_missing_funcs(ec_wine_mgs, possible_mgs_ecs)
ec_wine_mgs_nomiss_subset <- ec_wine_mgs_nomiss[overlapping_ec, overlapping_samples]

ec_wine_ITS_nsti2 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)
ec_wine_ITS_nsti1.5 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)
ec_wine_ITS_nsti1 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_1.0/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)
ec_wine_ITS_nsti0.5 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)
ec_wine_ITS_nsti0.25 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)
ec_wine_ITS_nsti0.1 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.1/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)
ec_wine_ITS_nsti0.05 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv", possible_ITS_ecs, overlapping_ec, overlapping_samples)

wine_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_wine_mgs_nomiss_subset)
wine_ec_mgs_null_df_round <- round(wine_ec_mgs_null_df  - 0.00000001)

wine_ec_mgs_null_scc <- cor_all_cols(tab1 = wine_ec_mgs_null_df, tab2 = ec_wine_mgs_nomiss_subset, cat_string="Null", metric="spearman")
wine_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, wine_ec_mgs_null_df_round, category="Null")

wine_ec_mgs_nsti2_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti2, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
wine_ec_mgs_nsti2_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti2, category="NSTI=2")
rownames(wine_ec_mgs_nsti2_metrics) <- wine_ec_mgs_nsti2_metrics$sample
rownames(wine_ec_mgs_nsti2_scc) <- wine_ec_mgs_nsti2_scc$sample_names

wine_ec_mgs_nsti1.5_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti1.5, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
wine_ec_mgs_nsti1.5_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti1.5, category="NSTI=1.5")
rownames(wine_ec_mgs_nsti1.5_metrics) <- wine_ec_mgs_nsti1.5_metrics$sample
rownames(wine_ec_mgs_nsti1.5_scc) <- wine_ec_mgs_nsti1.5_scc$sample_names

wine_ec_mgs_nsti1_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti1, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
wine_ec_mgs_nsti1_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti1, category="NSTI=1")
rownames(wine_ec_mgs_nsti1_metrics) <- wine_ec_mgs_nsti1_metrics$sample
rownames(wine_ec_mgs_nsti1_scc) <- wine_ec_mgs_nsti1_scc$sample_names

wine_ec_mgs_nsti0.5_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.5, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
wine_ec_mgs_nsti0.5_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.5, category="NSTI=0.5")
rownames(wine_ec_mgs_nsti0.5_metrics) <- wine_ec_mgs_nsti0.5_metrics$sample
rownames(wine_ec_mgs_nsti0.5_scc) <- wine_ec_mgs_nsti0.5_scc$sample_names

wine_ec_mgs_nsti0.25_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.25, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
wine_ec_mgs_nsti0.25_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.25, category="NSTI=0.25")
rownames(wine_ec_mgs_nsti0.25_metrics) <- wine_ec_mgs_nsti0.25_metrics$sample
rownames(wine_ec_mgs_nsti0.25_scc) <- wine_ec_mgs_nsti0.25_scc$sample_names

wine_ec_mgs_nsti0.1_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.1, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
wine_ec_mgs_nsti0.1_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.1, category="NSTI=0.1")
rownames(wine_ec_mgs_nsti0.1_metrics) <- wine_ec_mgs_nsti0.1_metrics$sample
rownames(wine_ec_mgs_nsti0.1_scc) <- wine_ec_mgs_nsti0.1_scc$sample_names

wine_ec_mgs_nsti0.05_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.05, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
wine_ec_mgs_nsti0.05_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.05, category="NSTI=0.05")
rownames(wine_ec_mgs_nsti0.05_metrics) <- wine_ec_mgs_nsti0.05_metrics$sample
rownames(wine_ec_mgs_nsti0.05_scc) <- wine_ec_mgs_nsti0.05_scc$sample_names


# Subset to pathways and samples overlapping and calculate spearman and accuracy metrics based on pathway abundances.
pathabun_wine_mgs_nomiss <- add_missing_funcs(pathabun_wine_mgs, all_possible_ITS_paths)
pathabun_wine_mgs_nomiss_subset <- pathabun_wine_mgs_nomiss[all_possible_ITS_paths, overlapping_samples]

pathabun_wine_ITS_nsti2 <- read_in_wine_ITS_picrust2("pathways_out_nsti_2.0/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)
pathabun_wine_ITS_nsti1.5 <- read_in_wine_ITS_picrust2("pathways_out_nsti_1.5/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)
pathabun_wine_ITS_nsti1 <- read_in_wine_ITS_picrust2("pathways_out_nsti_1.0/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)
pathabun_wine_ITS_nsti0.5 <- read_in_wine_ITS_picrust2("pathways_out_nsti_0.5/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)
pathabun_wine_ITS_nsti0.25 <- read_in_wine_ITS_picrust2("pathways_out_nsti_0.25/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)
pathabun_wine_ITS_nsti0.1 <- read_in_wine_ITS_picrust2("pathways_out_nsti_0.1/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)
pathabun_wine_ITS_nsti0.05 <- read_in_wine_ITS_picrust2("pathways_out_nsti_0.05/path_abun_unstrat.tsv", all_possible_ITS_paths, all_possible_ITS_paths, overlapping_samples)

wine_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun_ITS_no_miss, tab = pathabun_wine_mgs_nomiss_subset)
wine_pathabun_mgs_null_df_round <- round(wine_pathabun_mgs_null_df  - 0.00000001)

wine_pathabun_mgs_null_scc <- cor_all_cols(tab1 = wine_pathabun_mgs_null_df, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="Null", metric="spearman")
wine_pathabun_mgs_null_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, wine_pathabun_mgs_null_df_round, category="Null")

wine_pathabun_mgs_nsti2_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti2, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
wine_pathabun_mgs_nsti2_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti2, category="NSTI=2")
rownames(wine_pathabun_mgs_nsti2_metrics) <- wine_pathabun_mgs_nsti2_metrics$sample
rownames(wine_pathabun_mgs_nsti2_scc) <- wine_pathabun_mgs_nsti2_scc$sample_names

wine_pathabun_mgs_nsti1.5_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti1.5, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
wine_pathabun_mgs_nsti1.5_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti1.5, category="NSTI=1.5")
rownames(wine_pathabun_mgs_nsti1.5_metrics) <- wine_pathabun_mgs_nsti1.5_metrics$sample
rownames(wine_pathabun_mgs_nsti1.5_scc) <- wine_pathabun_mgs_nsti1.5_scc$sample_names

wine_pathabun_mgs_nsti1_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti1, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
wine_pathabun_mgs_nsti1_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti1, category="NSTI=1")
rownames(wine_pathabun_mgs_nsti1_metrics) <- wine_pathabun_mgs_nsti1_metrics$sample
rownames(wine_pathabun_mgs_nsti1_scc) <- wine_pathabun_mgs_nsti1_scc$sample_names

wine_pathabun_mgs_nsti0.5_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti0.5, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
wine_pathabun_mgs_nsti0.5_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti0.5, category="NSTI=0.5")
rownames(wine_pathabun_mgs_nsti0.5_metrics) <- wine_pathabun_mgs_nsti0.5_metrics$sample
rownames(wine_pathabun_mgs_nsti0.5_scc) <- wine_pathabun_mgs_nsti0.5_scc$sample_names

wine_pathabun_mgs_nsti0.25_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti0.25, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
wine_pathabun_mgs_nsti0.25_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti0.25, category="NSTI=0.25")
rownames(wine_pathabun_mgs_nsti0.25_metrics) <- wine_pathabun_mgs_nsti0.25_metrics$sample
rownames(wine_pathabun_mgs_nsti0.25_scc) <- wine_pathabun_mgs_nsti0.25_scc$sample_names

wine_pathabun_mgs_nsti0.1_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti0.1, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
wine_pathabun_mgs_nsti0.1_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti0.1, category="NSTI=0.1")
rownames(wine_pathabun_mgs_nsti0.1_metrics) <- wine_pathabun_mgs_nsti0.1_metrics$sample
rownames(wine_pathabun_mgs_nsti0.1_scc) <- wine_pathabun_mgs_nsti0.1_scc$sample_names

wine_pathabun_mgs_nsti0.05_scc <- cor_all_cols(tab1 = pathabun_wine_ITS_nsti0.05, tab2 = pathabun_wine_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
wine_pathabun_mgs_nsti0.05_metrics <- calc_accuracy_metrics(pathabun_wine_mgs_nomiss_subset, pathabun_wine_ITS_nsti0.05, category="NSTI=0.05")
rownames(wine_pathabun_mgs_nsti0.05_metrics) <- wine_pathabun_mgs_nsti0.05_metrics$sample
rownames(wine_pathabun_mgs_nsti0.05_scc) <- wine_pathabun_mgs_nsti0.05_scc$sample_names


wine_ec_mgs_scc_df <- rbind(wine_ec_mgs_null_scc,
                            wine_ec_mgs_nsti2_scc,
                            wine_ec_mgs_nsti1.5_scc,
                            wine_ec_mgs_nsti1_scc,
                            wine_ec_mgs_nsti0.5_scc,
                            wine_ec_mgs_nsti0.25_scc,
                            wine_ec_mgs_nsti0.1_scc,
                            wine_ec_mgs_nsti0.05_scc)

wine_ec_mgs_acc_df <- rbind(wine_ec_mgs_null_metrics,
                            wine_ec_mgs_nsti2_metrics,
                            wine_ec_mgs_nsti1.5_metrics,
                            wine_ec_mgs_nsti1_metrics,
                            wine_ec_mgs_nsti0.5_metrics,
                            wine_ec_mgs_nsti0.25_metrics,
                            wine_ec_mgs_nsti0.1_metrics,
                            wine_ec_mgs_nsti0.05_metrics)



wine_pathabun_mgs_scc_df <- rbind(wine_pathabun_mgs_null_scc,
                                  wine_pathabun_mgs_nsti2_scc,
                                  wine_pathabun_mgs_nsti1.5_scc,
                                  wine_pathabun_mgs_nsti1_scc,
                                  wine_pathabun_mgs_nsti0.5_scc,
                                  wine_pathabun_mgs_nsti0.25_scc,
                                  wine_pathabun_mgs_nsti0.1_scc,
                                  wine_pathabun_mgs_nsti0.05_scc)

wine_pathabun_mgs_acc_df <- rbind(wine_pathabun_mgs_null_metrics,
                                  wine_pathabun_mgs_nsti2_metrics,
                                  wine_pathabun_mgs_nsti1.5_metrics,
                                  wine_pathabun_mgs_nsti1_metrics,
                                  wine_pathabun_mgs_nsti0.5_metrics,
                                  wine_pathabun_mgs_nsti0.25_metrics,
                                  wine_pathabun_mgs_nsti0.1_metrics,
                                  wine_pathabun_mgs_nsti0.05_metrics)

saveRDS(object = wine_ec_mgs_scc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/wine_ec_scc_metrics.rds")

saveRDS(object = wine_ec_mgs_acc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/wine_ec_acc_metrics.rds")

saveRDS(object = wine_pathabun_mgs_scc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/wine_pathabun_scc_metrics.rds")

saveRDS(object = wine_pathabun_mgs_acc_df,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/18S_ITS_vs_MGS_metrics/wine_pathabun_acc_metrics.rds")
