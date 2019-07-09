### Specifically run pathway validations on PAPRICA MetaCyc pathways separately since there is not much intersect in the pathways identified as present.

rm(list=ls(all=TRUE))

library(ggplot2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

### Read in expected metacyc database (based on reference E.C. number database).
pathabun <- data.frame(t(read.table(gzfile("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_16S_pathway/path_abun_unstrat.tsv"),
                         row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

# Read in possible pathways predicted by PICRUSt2 (which is limited to same set as HUMAnN2 except limited to prokaryotes).
possible_picrust2_pathways <- read.table("possible_metacyc_pathways/picrust2_prokaryotic_pathways.txt", header=F, stringsAsFactors = FALSE)$V1
possible_paprica_pathways <- read.table("possible_metacyc_pathways/paprica_bacteria_pathways.txt", header=F, stringsAsFactors = FALSE)$V1
possible_pathways <- possible_picrust2_pathways[which(possible_picrust2_pathways %in% possible_paprica_pathways)]

length(possible_pathways)
# Only 193 overlapping pathways. The majority of the pathways specific to PAPRICA appear to be encoded by a small number of genes and were filtered out in HUMAnN2.

hmp_paprica_pathabun <- read_table_check_exists("paprica_out/hmp_paprica_pathabun.tsv", header=T, sep="\t", row.names=1)
hmp_paprica_pathabun <- hmp_paprica_pathabun[which(rownames(hmp_paprica_pathabun) %in% possible_pathways), ]
hmp_mgs_pathabun <- read_table_check_exists("../mgs_validation/hmp/humann2_pathabun_unstrat.tsv", header=T, sep="\t", row.names=1)
hmp_mgs_pathabun <- hmp_mgs_pathabun[which(rownames(hmp_mgs_pathabun) %in% possible_pathways), ]
hmp_overlapping_samples <- colnames(hmp_paprica_pathabun)[which(colnames(hmp_paprica_pathabun) %in% colnames(hmp_mgs_pathabun))]

hmp_paprica_pathabun <- hmp_paprica_pathabun[, hmp_overlapping_samples]
hmp_mgs_pathabun <- hmp_mgs_pathabun[, hmp_overlapping_samples]

hmp_paprica_pathabun <- add_missing_funcs(hmp_paprica_pathabun, possible_pathways)
hmp_mgs_pathabun <- add_missing_funcs(hmp_mgs_pathabun, possible_pathways)

hmp_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = hmp_mgs_pathabun)
hmp_pathabun_mgs_null_df_round <- round(hmp_pathabun_mgs_null_df  - 0.00000001)


hmp_pathabun_mgs_null <- cor_all_cols(tab1 = hmp_pathabun_mgs_null_df, tab2 = hmp_mgs_pathabun, cat_string="Null", metric="spearman")
hmp_pathabun_paprica_vs_mgs <- cor_all_cols(tab1 = hmp_paprica_pathabun, tab2 = hmp_mgs_pathabun, cat_string="PAPRICA", metric="spearman")

hmp_null_pathabun_metrics <- calc_accuracy_metrics(hmp_mgs_pathabun, hmp_pathabun_mgs_null_df_round, category="Null")
hmp_paprica_pathabun_metrics <- calc_accuracy_metrics(hmp_paprica_pathabun, hmp_pathabun_mgs_null_df_round, category="PAPRICA")

hmp_pathabun_spearman_df <- rbind(hmp_pathabun_mgs_null, hmp_pathabun_paprica_vs_mgs)
hmp_pathabun_acc_df <- rbind(hmp_null_pathabun_metrics, hmp_paprica_pathabun_metrics)


mammal_paprica_pathabun <- read_table_check_exists("paprica_out/mammal_paprica_pathabun.tsv", header=T, sep="\t", row.names=1)
mammal_paprica_pathabun <- mammal_paprica_pathabun[which(rownames(mammal_paprica_pathabun) %in% possible_pathways), ]
mammal_mgs_pathabun <- read_table_check_exists("../mgs_validation/mammal/humann2_pathabun_unstrat.tsv", header=T, sep="\t", row.names=1)
mammal_mgs_pathabun <- mammal_mgs_pathabun[which(rownames(mammal_mgs_pathabun) %in% possible_pathways), ]
mammal_overlapping_samples <- colnames(mammal_paprica_pathabun)[which(colnames(mammal_paprica_pathabun) %in% colnames(mammal_mgs_pathabun))]

mammal_paprica_pathabun <- mammal_paprica_pathabun[, mammal_overlapping_samples]
mammal_mgs_pathabun <- mammal_mgs_pathabun[, mammal_overlapping_samples]

mammal_paprica_pathabun <- add_missing_funcs(mammal_paprica_pathabun, possible_pathways)
mammal_mgs_pathabun <- add_missing_funcs(mammal_mgs_pathabun, possible_pathways)

mammal_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = mammal_mgs_pathabun)
mammal_pathabun_mgs_null_df_round <- round(mammal_pathabun_mgs_null_df  - 0.00000001)


mammal_pathabun_mgs_null <- cor_all_cols(tab1 = mammal_pathabun_mgs_null_df, tab2 = mammal_mgs_pathabun, cat_string="Null", metric="spearman")
mammal_pathabun_paprica_vs_mgs <- cor_all_cols(tab1 = mammal_paprica_pathabun, tab2 = mammal_mgs_pathabun, cat_string="PAPRICA", metric="spearman")

mammal_null_pathabun_metrics <- calc_accuracy_metrics(mammal_mgs_pathabun, mammal_pathabun_mgs_null_df_round, category="Null")
mammal_paprica_pathabun_metrics <- calc_accuracy_metrics(mammal_paprica_pathabun, mammal_pathabun_mgs_null_df_round, category="PAPRICA")

mammal_pathabun_spearman_df <- rbind(mammal_pathabun_mgs_null, mammal_pathabun_paprica_vs_mgs)
mammal_pathabun_acc_df <- rbind(mammal_null_pathabun_metrics, mammal_paprica_pathabun_metrics)


ocean_paprica_pathabun <- read_table_check_exists("paprica_out/ocean_paprica_pathabun.tsv", header=T, sep="\t", row.names=1)
ocean_paprica_pathabun <- ocean_paprica_pathabun[which(rownames(ocean_paprica_pathabun) %in% possible_pathways), ]
ocean_mgs_pathabun <- read_table_check_exists("../mgs_validation/ocean/humann2_pathabun_unstrat.tsv", header=T, sep="\t", row.names=1)
ocean_mgs_pathabun <- ocean_mgs_pathabun[which(rownames(ocean_mgs_pathabun) %in% possible_pathways), ]
ocean_overlapping_samples <- colnames(ocean_paprica_pathabun)[which(colnames(ocean_paprica_pathabun) %in% colnames(ocean_mgs_pathabun))]

ocean_paprica_pathabun <- ocean_paprica_pathabun[, ocean_overlapping_samples]
ocean_mgs_pathabun <- ocean_mgs_pathabun[, ocean_overlapping_samples]

ocean_paprica_pathabun <- add_missing_funcs(ocean_paprica_pathabun, possible_pathways)
ocean_mgs_pathabun <- add_missing_funcs(ocean_mgs_pathabun, possible_pathways)

ocean_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = ocean_mgs_pathabun)
ocean_pathabun_mgs_null_df_round <- round(ocean_pathabun_mgs_null_df  - 0.00000001)


ocean_pathabun_mgs_null <- cor_all_cols(tab1 = ocean_pathabun_mgs_null_df, tab2 = ocean_mgs_pathabun, cat_string="Null", metric="spearman")
ocean_pathabun_paprica_vs_mgs <- cor_all_cols(tab1 = ocean_paprica_pathabun, tab2 = ocean_mgs_pathabun, cat_string="PAPRICA", metric="spearman")

ocean_null_pathabun_metrics <- calc_accuracy_metrics(ocean_mgs_pathabun, ocean_pathabun_mgs_null_df_round, category="Null")
ocean_paprica_pathabun_metrics <- calc_accuracy_metrics(ocean_paprica_pathabun, ocean_pathabun_mgs_null_df_round, category="PAPRICA")

ocean_pathabun_spearman_df <- rbind(ocean_pathabun_mgs_null, ocean_pathabun_paprica_vs_mgs)
ocean_pathabun_acc_df <- rbind(ocean_null_pathabun_metrics, ocean_paprica_pathabun_metrics)


blueberry_paprica_pathabun <- read_table_check_exists("paprica_out/blueberry_paprica_pathabun.tsv", header=T, sep="\t", row.names=1)
blueberry_paprica_pathabun <- blueberry_paprica_pathabun[which(rownames(blueberry_paprica_pathabun) %in% possible_pathways), ]
blueberry_mgs_pathabun <- read_table_check_exists("../mgs_validation/blueberry/humann2_pathabun_unstrat.tsv", header=T, sep="\t", row.names=1)
blueberry_mgs_pathabun <- blueberry_mgs_pathabun[which(rownames(blueberry_mgs_pathabun) %in% possible_pathways), ]
blueberry_overlapping_samples <- colnames(blueberry_paprica_pathabun)[which(colnames(blueberry_paprica_pathabun) %in% colnames(blueberry_mgs_pathabun))]

blueberry_paprica_pathabun <- blueberry_paprica_pathabun[, blueberry_overlapping_samples]
blueberry_mgs_pathabun <- blueberry_mgs_pathabun[, blueberry_overlapping_samples]

blueberry_paprica_pathabun <- add_missing_funcs(blueberry_paprica_pathabun, possible_pathways)
blueberry_mgs_pathabun <- add_missing_funcs(blueberry_mgs_pathabun, possible_pathways)

blueberry_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = blueberry_mgs_pathabun)
blueberry_pathabun_mgs_null_df_round <- round(blueberry_pathabun_mgs_null_df  - 0.00000001)


blueberry_pathabun_mgs_null <- cor_all_cols(tab1 = blueberry_pathabun_mgs_null_df, tab2 = blueberry_mgs_pathabun, cat_string="Null", metric="spearman")
blueberry_pathabun_paprica_vs_mgs <- cor_all_cols(tab1 = blueberry_paprica_pathabun, tab2 = blueberry_mgs_pathabun, cat_string="PAPRICA", metric="spearman")

blueberry_null_pathabun_metrics <- calc_accuracy_metrics(blueberry_mgs_pathabun, blueberry_pathabun_mgs_null_df_round, category="Null")
blueberry_paprica_pathabun_metrics <- calc_accuracy_metrics(blueberry_paprica_pathabun, blueberry_pathabun_mgs_null_df_round, category="PAPRICA")

blueberry_pathabun_spearman_df <- rbind(blueberry_pathabun_mgs_null, blueberry_pathabun_paprica_vs_mgs)
blueberry_pathabun_acc_df <- rbind(blueberry_null_pathabun_metrics, blueberry_paprica_pathabun_metrics)




saveRDS(object = hmp_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_pathabun_PAPRICA_spearman_df.rds")
saveRDS(object = hmp_pathabun_acc_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_pathabun_PAPRICA_acc_metrics.rds")

saveRDS(object = mammal_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_pathabun_PAPRICA_spearman_df.rds")
saveRDS(object = mammal_pathabun_acc_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_pathabun_PAPRICA_acc_metrics.rds")

saveRDS(object = ocean_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_pathabun_PAPRICA_spearman_df.rds")
saveRDS(object = ocean_pathabun_acc_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_pathabun_PAPRICA_acc_metrics.rds")

saveRDS(object = blueberry_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_pathabun_PAPRICA_spearman_df.rds")
saveRDS(object = blueberry_pathabun_acc_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_pathabun_PAPRICA_acc_metrics.rds")
