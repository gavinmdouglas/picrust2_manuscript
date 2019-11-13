### Specifically run pathway validations on PAPRICA MetaCyc pathways
### separately since there is not much intersect in the pathways identified as present.

rm(list=ls(all=TRUE))

library(ggplot2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

calc_paprica_path_metrics <- function(dataset, possible_path, in_db) {

  paprica_pathabun <- read_table_check_exists(paste("paprica_out/", dataset, "_paprica_path.tsv", sep=""), header=T, sep="\t", row.names=1)
  paprica_pathabun <- paprica_pathabun[which(rownames(paprica_pathabun) %in% possible_path), ]

  mgs_pathabun <- read_table_check_exists(paste("../mgs_validation/", dataset, "/humann2_pathabun_unstrat.tsv", sep=""), header=T, sep="\t", row.names=1)
  mgs_pathabun <- mgs_pathabun[which(rownames(mgs_pathabun) %in% possible_path), ]

  overlapping_samples <- colnames(paprica_pathabun)[which(colnames(paprica_pathabun) %in% colnames(mgs_pathabun))]
  
  paprica_pathabun <- paprica_pathabun[, overlapping_samples]
  mgs_pathabun <- mgs_pathabun[, overlapping_samples]
  
  paprica_pathabun <- add_missing_funcs(paprica_pathabun, possible_path)
  mgs_pathabun <- add_missing_funcs(mgs_pathabun, possible_path)
  
  pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = in_db, tab = mgs_pathabun)
  pathabun_mgs_null_df_round <- round(pathabun_mgs_null_df  - 0.00000001)

  pathabun_mgs_null <- cor_all_cols(tab1 = pathabun_mgs_null_df, tab2 = mgs_pathabun, cat_string="Null", metric="spearman")
  pathabun_paprica_vs_mgs <- cor_all_cols(tab1 = paprica_pathabun, tab2 = mgs_pathabun, cat_string="PAPRICA", metric="spearman")
  
  null_pathabun_metrics <- calc_accuracy_metrics(mgs_pathabun, pathabun_mgs_null_df_round, category="Null")
  paprica_pathabun_metrics <- calc_accuracy_metrics(paprica_pathabun, pathabun_mgs_null_df_round, category="PAPRICA")
  
  return(list("spearman"=rbind(pathabun_mgs_null, pathabun_paprica_vs_mgs),
              "acc"=rbind(null_pathabun_metrics, paprica_pathabun_metrics)))
}

# Read in expected metacyc database (based on reference EC number database).
path_db <- data.frame(t(read.table(gzfile("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_16S_pathway/path_abun_unstrat.tsv"),
                                    row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

# Read in possible pathways output by HUMAnN2 and PAPRICA.
poss_mgs_pathways <- read.table("possible_path/humann2_path.txt", header=FALSE, stringsAsFactors = FALSE)$V1
poss_paprica_pathways <- read.table("possible_path/paprica_path.txt", header=FALSE, stringsAsFactors = FALSE)$V1

poss_pathways <- poss_mgs_pathways[which(poss_mgs_pathways %in% poss_paprica_pathways)]

datasets <- c("cameroon", "hmp", "indian", "primate", "mammal", "ocean", "blueberry")

paprica_path_performance <- list()

for(d in datasets) {
 
  paprica_path_performance[[d]] <- calc_paprica_path_metrics(dataset=d, possible_path=poss_pathways, in_db=path_db)
}
  
saveRDS(object = paprica_path_performance,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/paprica_path_metrics.rds")
