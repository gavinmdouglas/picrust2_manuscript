setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

# Read in counts of each kingdom for each sample.
kingdom_count <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/sortmerna_out/kingdom_counts.tsv",
                            header=F, sep="\t", row.names=1)
colnames(kingdom_count) <- c("archaea", "bacteria", "eukaryota", "mixed")
rownames(kingdom_count) <- gsub("_16S_sortmerna_out.fastq.blast", "", rownames(kingdom_count))

kingdom_percent <- sweep(kingdom_count, 1, rowSums(kingdom_count), '/') * 100

# Identify possible ECs.
possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt")$V1
possible_18S_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/picrust2_ECs_18S.txt")$V1

ec_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_18S.txt",
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

ec_blue_mgs <- read.table("mgs/humann2_output_final/humann2_level4ec_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")
rownames(ec_blue_mgs) <- gsub("^", "EC:", rownames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("_Abundance.RPKs", "", colnames(ec_blue_mgs))
ec_blue_mgs <- ec_blue_mgs[-which(rownames(ec_blue_mgs) %in% c("EC:UNMAPPED", "EC:UNGROUPED")),]
ec_blue_mgs_nomiss <- add_missing_funcs(ec_blue_mgs, possible_mgs_ecs)

ec_blue_18S_nsti0.05 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti0.05) <- gsub("^", "BB", colnames(ec_blue_18S_nsti0.05))
ec_blue_18S_nsti0.05_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.05, possible_18S_ecs)

ec_blue_18S_nsti0.1 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_0.1/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti0.1) <- gsub("^", "BB", colnames(ec_blue_18S_nsti0.1))
ec_blue_18S_nsti0.1_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.1, possible_18S_ecs)

ec_blue_18S_nsti0.25 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti0.25) <- gsub("^", "BB", colnames(ec_blue_18S_nsti0.25))
ec_blue_18S_nsti0.25_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.25, possible_18S_ecs)

ec_blue_18S_nsti0.5 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti0.5) <- gsub("^", "BB", colnames(ec_blue_18S_nsti0.5))
ec_blue_18S_nsti0.5_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.5, possible_18S_ecs)

ec_blue_18S_nsti1 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_1/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti1) <- gsub("^", "BB", colnames(ec_blue_18S_nsti1))
ec_blue_18S_nsti1_nomiss <- add_missing_funcs(ec_blue_18S_nsti1, possible_18S_ecs)

ec_blue_18S_nsti1.5 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti1.5) <- gsub("^", "BB", colnames(ec_blue_18S_nsti1.5))
ec_blue_18S_nsti1.5_nomiss <- add_missing_funcs(ec_blue_18S_nsti1.5, possible_18S_ecs)

ec_blue_18S_nsti2 <- read.table("18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(ec_blue_18S_nsti2) <- gsub("^", "BB", colnames(ec_blue_18S_nsti2))
ec_blue_18S_nsti2_nomiss <- add_missing_funcs(ec_blue_18S_nsti2, possible_18S_ecs)

overlapping_ec <- rownames(ec_blue_18S_nsti2_nomiss)[which(rownames(ec_blue_18S_nsti2_nomiss) %in% rownames(ec_blue_mgs_nomiss))]
overlapping_col <- colnames(ec_blue_18S_nsti0.05_nomiss)[which(colnames(ec_blue_18S_nsti0.05_nomiss) %in% colnames(ec_blue_mgs_nomiss))]

ec_blue_mgs_nomiss_subset <- ec_blue_mgs_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti0.05_nomiss_subset <- ec_blue_18S_nsti0.05_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti0.1_nomiss_subset <- ec_blue_18S_nsti0.1_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti0.25_nomiss_subset <- ec_blue_18S_nsti0.25_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti0.5_nomiss_subset <- ec_blue_18S_nsti0.5_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti1_nomiss_subset <- ec_blue_18S_nsti1_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti1.5_nomiss_subset <- ec_blue_18S_nsti1.5_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti2_nomiss_subset <- ec_blue_18S_nsti2_nomiss[overlapping_ec, overlapping_col]

# Test:
blue_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_blue_mgs_nomiss_subset)
blue_ec_mgs_null_df_round <- round(blue_ec_mgs_null_df  - 0.00000001)

blue_ec_mgs_null_scc <- cor_all_cols(tab1 = blue_ec_mgs_null_df, tab2 = ec_blue_mgs_nomiss_subset, cat_string="Null", metric="spearman")
blue_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, blue_ec_mgs_null_df, category="Null")

blue_ec_mgs_nsti0.05_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti0.05_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
blue_ec_mgs_nsti0.05_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti0.05_nomiss_subset, category="NSTI=0.05")

blue_ec_mgs_nsti0.1_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti0.1_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
blue_ec_mgs_nsti0.1_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti0.1_nomiss_subset, category="NSTI=0.1")

blue_ec_mgs_nsti0.25_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti0.25_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
blue_ec_mgs_nsti0.25_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti0.25_nomiss_subset, category="NSTI=0.25")

blue_ec_mgs_nsti0.5_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti0.5_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
blue_ec_mgs_nsti0.5_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti0.5_nomiss_subset, category="NSTI=0.5")

blue_ec_mgs_nsti1_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti1_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
blue_ec_mgs_nsti1_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti1_nomiss_subset, category="NSTI=1")

blue_ec_mgs_nsti1.5_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti1.5_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
blue_ec_mgs_nsti1.5_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti1.5_nomiss_subset, category="NSTI=1.5")

blue_ec_mgs_nsti2_scc <- cor_all_cols(tab1 = ec_blue_18S_nsti2_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
blue_ec_mgs_nsti2_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nsti2_nomiss_subset, category="NSTI=2")




possible_mgs_pathways <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_metacyc_pathways/humann2_pathways.txt")$V1
possible_18S_pathways <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_metacyc_pathways/picrust2_eukaryotic_pathways.txt")$V1

pathabun_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ref_genome_metacyc_18S/path_abun_unstrat.tsv",
                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)
pathabun_18S_nomiss <- data.frame(t(add_missing_funcs(pathabun_18S, possible_18S_pathways)), check.names = FALSE)
  
pathabun_blue_mgs <- read.table("mgs/humann2_output_final/humann2_pathabundance_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")
rownames(pathabun_blue_mgs) <- gsub(":.*$", "", rownames(pathabun_blue_mgs))
colnames(pathabun_blue_mgs) <- gsub("_Abundance", "", colnames(pathabun_blue_mgs))
pathabun_blue_mgs <- pathabun_blue_mgs[-which(rownames(pathabun_blue_mgs) %in% c("UNMAPPED", "UNGROUPED", "UNINTEGRATED")),]
pathabun_blue_mgs_nomiss <- add_missing_funcs(pathabun_blue_mgs, possible_mgs_pathways)

pathabun_blue_18S_nsti0.05 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.05/path_abun_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti0.05) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti0.05))
pathabun_blue_18S_nsti0.05_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti0.05, possible_18S_pathways)

pathabun_blue_18S_nsti0.1 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.1/path_abun_unstrat.tsv",
                                  header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti0.1) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti0.1))
pathabun_blue_18S_nsti0.1_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti0.1, possible_18S_pathways)

pathabun_blue_18S_nsti0.25 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.25/path_abun_unstrat.tsv",
                                   header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti0.25) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti0.25))
pathabun_blue_18S_nsti0.25_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti0.25, possible_18S_pathways)

pathabun_blue_18S_nsti0.5 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.5/path_abun_unstrat.tsv",
                                  header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti0.5) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti0.5))
pathabun_blue_18S_nsti0.5_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti0.5, possible_18S_pathways)

pathabun_blue_18S_nsti1 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_1/path_abun_unstrat.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti1) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti1))
pathabun_blue_18S_nsti1_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti1, possible_18S_pathways)

pathabun_blue_18S_nsti1.5 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_1.5/path_abun_unstrat.tsv",
                                  header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti1.5) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti1.5))
pathabun_blue_18S_nsti1.5_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti1.5, possible_18S_pathways)

pathabun_blue_18S_nsti2 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_2/path_abun_unstrat.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathabun_blue_18S_nsti2) <- gsub("^", "BB", colnames(pathabun_blue_18S_nsti2))
pathabun_blue_18S_nsti2_nomiss <- add_missing_funcs(pathabun_blue_18S_nsti2, possible_18S_pathways)

overlapping_pathabun <- rownames(pathabun_blue_18S_nsti2_nomiss)[which(rownames(pathabun_blue_18S_nsti2_nomiss) %in% rownames(pathabun_blue_mgs_nomiss))]
overlapping_col <- colnames(pathabun_blue_18S_nsti0.05_nomiss)[which(colnames(pathabun_blue_18S_nsti0.05_nomiss) %in% colnames(pathabun_blue_mgs_nomiss))]

pathabun_blue_mgs_nomiss_subset <- pathabun_blue_mgs_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti0.05_nomiss_subset <- pathabun_blue_18S_nsti0.05_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti0.1_nomiss_subset <- pathabun_blue_18S_nsti0.1_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti0.25_nomiss_subset <- pathabun_blue_18S_nsti0.25_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti0.5_nomiss_subset <- pathabun_blue_18S_nsti0.5_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti1_nomiss_subset <- pathabun_blue_18S_nsti1_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti1.5_nomiss_subset <- pathabun_blue_18S_nsti1.5_nomiss[overlapping_pathabun, overlapping_col]
pathabun_blue_18S_nsti2_nomiss_subset <- pathabun_blue_18S_nsti2_nomiss[overlapping_pathabun, overlapping_col]

blue_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun_18S_nomiss, tab = pathabun_blue_mgs_nomiss_subset)
blue_pathabun_mgs_null_df_round <- round(blue_pathabun_mgs_null_df  - 0.00000001)

blue_pathabun_mgs_null_scc <- cor_all_cols(tab1 = blue_pathabun_mgs_null_df, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="Null", metric="spearman")
blue_pathabun_mgs_null_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, blue_pathabun_mgs_null_df, category="Null")

blue_pathabun_mgs_nsti0.05_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti0.05_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
blue_pathabun_mgs_nsti0.05_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti0.05_nomiss_subset, category="NSTI=0.05")

blue_pathabun_mgs_nsti0.1_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti0.1_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
blue_pathabun_mgs_nsti0.1_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti0.1_nomiss_subset, category="NSTI=0.1")

blue_pathabun_mgs_nsti0.25_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti0.25_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
blue_pathabun_mgs_nsti0.25_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti0.25_nomiss_subset, category="NSTI=0.25")

blue_pathabun_mgs_nsti0.5_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti0.5_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
blue_pathabun_mgs_nsti0.5_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti0.5_nomiss_subset, category="NSTI=0.5")

blue_pathabun_mgs_nsti1_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti1_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
blue_pathabun_mgs_nsti1_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti1_nomiss_subset, category="NSTI=1")

blue_pathabun_mgs_nsti1.5_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti1.5_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
blue_pathabun_mgs_nsti1.5_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti1.5_nomiss_subset, category="NSTI=1.5")

blue_pathabun_mgs_nsti2_scc <- cor_all_cols(tab1 = pathabun_blue_18S_nsti2_nomiss_subset, tab2 = pathabun_blue_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
blue_pathabun_mgs_nsti2_metrics <- calc_accuracy_metrics(pathabun_blue_mgs_nomiss_subset, pathabun_blue_18S_nsti2_nomiss_subset, category="NSTI=2")


# Path cov:
pathcov_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ref_genome_metacyc_18S/path_cov_unstrat.tsv",
                           row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)
pathcov_18S_nomiss <- data.frame(t(add_missing_funcs(pathcov_18S, possible_18S_pathways)), check.names = FALSE)

pathcov_blue_mgs <- read.table("mgs/humann2_output_final/humann2_pathcoverage_unstratified.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="")
rownames(pathcov_blue_mgs) <- gsub(":.*$", "", rownames(pathcov_blue_mgs))
colnames(pathcov_blue_mgs) <- gsub("_Coverage", "", colnames(pathcov_blue_mgs))
pathcov_blue_mgs <- pathcov_blue_mgs[-which(rownames(pathcov_blue_mgs) %in% c("UNMAPPED", "UNGROUPED", "UNINTEGRATED")),]
pathcov_blue_mgs_nomiss <- add_missing_funcs(pathcov_blue_mgs, possible_mgs_pathways)

pathcov_blue_18S_nsti0.05 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.05/path_cov_unstrat.tsv",
                                         header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti0.05) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti0.05))
pathcov_blue_18S_nsti0.05_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti0.05, possible_18S_pathways)

pathcov_blue_18S_nsti0.1 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.1/path_cov_unstrat.tsv",
                                        header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti0.1) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti0.1))
pathcov_blue_18S_nsti0.1_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti0.1, possible_18S_pathways)

pathcov_blue_18S_nsti0.25 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.25/path_cov_unstrat.tsv",
                                         header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti0.25) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti0.25))
pathcov_blue_18S_nsti0.25_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti0.25, possible_18S_pathways)

pathcov_blue_18S_nsti0.5 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_0.5/path_cov_unstrat.tsv",
                                        header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti0.5) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti0.5))
pathcov_blue_18S_nsti0.5_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti0.5, possible_18S_pathways)

pathcov_blue_18S_nsti1 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_1/path_cov_unstrat.tsv",
                                      header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti1) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti1))
pathcov_blue_18S_nsti1_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti1, possible_18S_pathways)

pathcov_blue_18S_nsti1.5 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_1.5/path_cov_unstrat.tsv",
                                        header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti1.5) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti1.5))
pathcov_blue_18S_nsti1.5_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti1.5, possible_18S_pathways)

pathcov_blue_18S_nsti2 <- read.table("18S/picrust2_full_output_pipeline/pathways_out_nsti_2/path_cov_unstrat.tsv",
                                      header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)
colnames(pathcov_blue_18S_nsti2) <- gsub("^", "BB", colnames(pathcov_blue_18S_nsti2))
pathcov_blue_18S_nsti2_nomiss <- add_missing_funcs(pathcov_blue_18S_nsti2, possible_18S_pathways)

overlapping_pathcov <- rownames(pathcov_blue_18S_nsti2_nomiss)[which(rownames(pathcov_blue_18S_nsti2_nomiss) %in% rownames(pathcov_blue_mgs_nomiss))]
overlapping_col <- colnames(pathcov_blue_18S_nsti0.05_nomiss)[which(colnames(pathcov_blue_18S_nsti0.05_nomiss) %in% colnames(pathcov_blue_mgs_nomiss))]

pathcov_blue_mgs_nomiss_subset <- pathcov_blue_mgs_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti0.05_nomiss_subset <- pathcov_blue_18S_nsti0.05_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti0.1_nomiss_subset <- pathcov_blue_18S_nsti0.1_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti0.25_nomiss_subset <- pathcov_blue_18S_nsti0.25_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti0.5_nomiss_subset <- pathcov_blue_18S_nsti0.5_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti1_nomiss_subset <- pathcov_blue_18S_nsti1_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti1.5_nomiss_subset <- pathcov_blue_18S_nsti1.5_nomiss[overlapping_pathcov, overlapping_col]
pathcov_blue_18S_nsti2_nomiss_subset <- pathcov_blue_18S_nsti2_nomiss[overlapping_pathcov, overlapping_col]

blue_pathcov_mgs_null_df <- generate_null_mean_db_funcs(db = pathcov_18S_nomiss, tab = pathcov_blue_mgs_nomiss_subset)
blue_pathcov_mgs_null_df_round <- round(blue_pathcov_mgs_null_df  - 0.00000001)

blue_pathcov_mgs_null_scc <- cor_all_cols(tab1 = blue_pathcov_mgs_null_df, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="Null", metric="spearman")
blue_pathcov_mgs_null_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, blue_pathcov_mgs_null_df, category="Null")

blue_pathcov_mgs_nsti0.05_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti0.05_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
blue_pathcov_mgs_nsti0.05_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti0.05_nomiss_subset, category="NSTI=0.05")

blue_pathcov_mgs_nsti0.1_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti0.1_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
blue_pathcov_mgs_nsti0.1_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti0.1_nomiss_subset, category="NSTI=0.1")

blue_pathcov_mgs_nsti0.25_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti0.25_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
blue_pathcov_mgs_nsti0.25_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti0.25_nomiss_subset, category="NSTI=0.25")

blue_pathcov_mgs_nsti0.5_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti0.5_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
blue_pathcov_mgs_nsti0.5_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti0.5_nomiss_subset, category="NSTI=0.5")

blue_pathcov_mgs_nsti1_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti1_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
blue_pathcov_mgs_nsti1_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti1_nomiss_subset, category="NSTI=1")

blue_pathcov_mgs_nsti1.5_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti1.5_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
blue_pathcov_mgs_nsti1.5_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti1.5_nomiss_subset, category="NSTI=1.5")

blue_pathcov_mgs_nsti2_scc <- cor_all_cols(tab1 = pathcov_blue_18S_nsti2_nomiss_subset, tab2 = pathcov_blue_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
blue_pathcov_mgs_nsti2_metrics <- calc_accuracy_metrics(pathcov_blue_mgs_nomiss_subset, pathcov_blue_18S_nsti2_nomiss_subset, category="NSTI=2")

ec_blueberry_combined_metrics <- rbind(blue_ec_mgs_null_metrics,
                                       blue_ec_mgs_nsti0.05_metrics,
                                       blue_ec_mgs_nsti0.1_metrics,
                                       blue_ec_mgs_nsti0.25_metrics,
                                       blue_ec_mgs_nsti0.5_metrics,
                                       blue_ec_mgs_nsti1_metrics,
                                       blue_ec_mgs_nsti1.5_metrics,
                                       blue_ec_mgs_nsti2_metrics)

ec_blueberry_combined_scc <- rbind(blue_ec_mgs_null_scc,
                                   blue_ec_mgs_nsti0.05_scc,
                                   blue_ec_mgs_nsti0.1_scc,
                                   blue_ec_mgs_nsti0.25_scc,
                                   blue_ec_mgs_nsti0.5_scc,
                                   blue_ec_mgs_nsti1_scc,
                                   blue_ec_mgs_nsti1.5_scc,
                                   blue_ec_mgs_nsti2_scc)

pathabun_blueberry_combined_metrics <- rbind(blue_pathabun_mgs_null_metrics,
                                       blue_pathabun_mgs_nsti0.05_metrics,
                                       blue_pathabun_mgs_nsti0.1_metrics,
                                       blue_pathabun_mgs_nsti0.25_metrics,
                                       blue_pathabun_mgs_nsti0.5_metrics,
                                       blue_pathabun_mgs_nsti1_metrics,
                                       blue_pathabun_mgs_nsti1.5_metrics,
                                       blue_pathabun_mgs_nsti2_metrics)

pathabun_blueberry_combined_scc <- rbind(blue_pathabun_mgs_null_scc,
                                   blue_pathabun_mgs_nsti0.05_scc,
                                   blue_pathabun_mgs_nsti0.1_scc,
                                   blue_pathabun_mgs_nsti0.25_scc,
                                   blue_pathabun_mgs_nsti0.5_scc,
                                   blue_pathabun_mgs_nsti1_scc,
                                   blue_pathabun_mgs_nsti1.5_scc,
                                   blue_pathabun_mgs_nsti2_scc)

pathcov_blueberry_combined_metrics <- rbind(blue_pathcov_mgs_null_metrics,
                                             blue_pathcov_mgs_nsti0.05_metrics,
                                             blue_pathcov_mgs_nsti0.1_metrics,
                                             blue_pathcov_mgs_nsti0.25_metrics,
                                             blue_pathcov_mgs_nsti0.5_metrics,
                                             blue_pathcov_mgs_nsti1_metrics,
                                             blue_pathcov_mgs_nsti1.5_metrics,
                                             blue_pathcov_mgs_nsti2_metrics)

pathcov_blueberry_combined_scc <- rbind(blue_pathcov_mgs_null_scc,
                                         blue_pathcov_mgs_nsti0.05_scc,
                                         blue_pathcov_mgs_nsti0.1_scc,
                                         blue_pathcov_mgs_nsti0.25_scc,
                                         blue_pathcov_mgs_nsti0.5_scc,
                                         blue_pathcov_mgs_nsti1_scc,
                                         blue_pathcov_mgs_nsti1.5_scc,
                                         blue_pathcov_mgs_nsti2_scc)

saveRDS(object = ec_blueberry_combined_metrics, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_18S_ec_acc_metrics.rds")
saveRDS(object = ec_blueberry_combined_scc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_18S_ec_spearman_df.rds")

saveRDS(object = pathabun_blueberry_combined_metrics, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_18S_pathabun_acc_metrics.rds")
saveRDS(object = pathabun_blueberry_combined_scc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_18S_pathabun_spearman_df.rds")

saveRDS(object = pathcov_blueberry_combined_metrics, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_18S_pathcov_acc_metrics.rds")
saveRDS(object = pathcov_blueberry_combined_scc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_18S_pathcov_spearman_df.rds")
