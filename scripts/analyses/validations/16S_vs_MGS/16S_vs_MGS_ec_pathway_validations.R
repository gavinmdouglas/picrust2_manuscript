### Comparing predicted EC numbers based on 16S to MGS "gold standard".
### Also compared to Paprica, which output EC number predictions.
### Calculated spearman correlations between MGS and 16S data and saved as RDS object.
### Also calculated accuracy metrics and saved as RDS.

library(ggplot2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

### Read in database files (used for calculating null distribution).
ec <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ec.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

### Read in expected metacyc database (based on reference E.C. number database).
pathabun <- data.frame(t(read.table(gzfile("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_16S_pathway/path_abun_unstrat.tsv"),
                         row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

pathcov <- data.frame(t(read.table(gzfile("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_16S_pathway/path_cov_unstrat.tsv"),
                                    row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

# Read in all tables, restrict to overlapping samples only, and get subsets with all possible ECs filled in and with
# only ECs overlapping between PICRUSt2 and PAPRICA. Focus analyses on ECs overlapping between both tools, but outputted others as well
# for sanity checks.
hmp_ec_infiles <- read_in_ec_predictions("hmp")
mammal_ec_infiles <- read_in_ec_predictions("mammal")
ocean_ec_infiles <- read_in_ec_predictions("ocean")
blueberry_ec_infiles <- read_in_ec_predictions("blueberry")

# Do same for pathway abundances and coverages (which are only output by PICRUSt2 and HUMAnN2).
hmp_pathway_infiles <- read_in_pathway_predictions("hmp")
mammal_pathway_infiles <- read_in_pathway_predictions("mammal")
ocean_pathway_infiles <- read_in_pathway_predictions("ocean")
blueberry_pathway_infiles <- read_in_pathway_predictions("blueberry")

# Generate random tables based on subsampling database to calculate null distributions.
hmp_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec, tab = hmp_ec_infiles$all_ecs_overlap$mgs_ec)
mammal_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec, tab = mammal_ec_infiles$all_ecs_overlap$mgs_ec)
ocean_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec, tab = ocean_ec_infiles$all_ecs_overlap$mgs_ec)
blueberry_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec, tab = blueberry_ec_infiles$all_ecs_overlap$mgs_ec)

hmp_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = hmp_pathway_infiles$all_pathabun$mgs_pathabun)
mammal_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = mammal_pathway_infiles$all_pathabun$mgs_pathabun)
ocean_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = ocean_pathway_infiles$all_pathabun$mgs_pathabun)
blueberry_pathabun_mgs_null_df <- generate_null_mean_db_funcs(db = pathabun, tab = blueberry_pathway_infiles$all_pathabun$mgs_pathabun)

hmp_pathcov_mgs_null_df <- generate_null_mean_db_funcs(db = pathcov, tab = hmp_pathway_infiles$all_pathcov$mgs_pathcov)
mammal_pathcov_mgs_null_df <- generate_null_mean_db_funcs(db = pathcov, tab = mammal_pathway_infiles$all_pathcov$mgs_pathcov)
ocean_pathcov_mgs_null_df <- generate_null_mean_db_funcs(db = pathcov, tab = ocean_pathway_infiles$all_pathcov$mgs_pathcov)
blueberry_pathcov_mgs_null_df <- generate_null_mean_db_funcs(db = pathcov, tab = blueberry_pathway_infiles$all_pathcov$mgs_pathcov)

hmp_ec_mgs_null_df_round <- round(hmp_ec_mgs_null_df  - 0.00000001)
mammal_ec_mgs_null_df_round <- round(mammal_ec_mgs_null_df  - 0.00000001)
ocean_ec_mgs_null_df_round <- round(ocean_ec_mgs_null_df  - 0.00000001)
blueberry_ec_mgs_null_df_round <- round(blueberry_ec_mgs_null_df  - 0.00000001)

hmp_pathabun_mgs_null_df_round <- round(hmp_pathabun_mgs_null_df  - 0.00000001)
mammal_pathabun_mgs_null_df_round <- round(mammal_pathabun_mgs_null_df  - 0.00000001)
ocean_pathabun_mgs_null_df_round <- round(ocean_pathabun_mgs_null_df  - 0.00000001)
blueberry_pathabun_mgs_null_df_round <- round(blueberry_pathabun_mgs_null_df  - 0.00000001)

hmp_pathcov_mgs_null_df_round <- round(hmp_pathcov_mgs_null_df  - 0.00000001)
mammal_pathcov_mgs_null_df_round <- round(mammal_pathcov_mgs_null_df  - 0.00000001)
ocean_pathcov_mgs_null_df_round <- round(ocean_pathcov_mgs_null_df  - 0.00000001)
blueberry_pathcov_mgs_null_df_round <- round(blueberry_pathcov_mgs_null_df  - 0.00000001)

# Get spearman correlations for each table with MGS:
hmp_ec_mgs_null <- cor_all_cols(tab1 = hmp_ec_mgs_null_df, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="Null", metric="spearman")
hmp_ec_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2 (GG)", metric="spearman")
hmp_ec_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2", metric="spearman")
hmp_ec_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1.5", metric="spearman")
hmp_ec_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1", metric="spearman")
hmp_ec_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.5", metric="spearman")
hmp_ec_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.25", metric="spearman")
hmp_ec_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.1", metric="spearman")
hmp_ec_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.05", metric="spearman")
hmp_ec_paprica_vs_mgs <- cor_all_cols(tab1 = hmp_ec_infiles$all_ecs_overlap$paprica_ec, tab2 = hmp_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="PAPRICA", metric="spearman")

hmp_pathabun_mgs_null <- cor_all_cols(tab1 = hmp_pathabun_mgs_null_df, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="Null", metric="spearman")
hmp_pathabun_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2 (GG)", metric="spearman")
hmp_pathabun_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2", metric="spearman")
hmp_pathabun_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1.5", metric="spearman")
hmp_pathabun_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1", metric="spearman")
hmp_pathabun_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.5", metric="spearman")
hmp_pathabun_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.25", metric="spearman")
hmp_pathabun_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.1", metric="spearman")
hmp_pathabun_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, tab2 = hmp_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.05", metric="spearman")

hmp_pathcov_mgs_null <- cor_all_cols(tab1 = hmp_pathcov_mgs_null_df, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="Null", metric="spearman")
hmp_pathcov_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2 (GG)", metric="spearman")
hmp_pathcov_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2", metric="spearman")
hmp_pathcov_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1.5", metric="spearman")
hmp_pathcov_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1", metric="spearman")
hmp_pathcov_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.5", metric="spearman")
hmp_pathcov_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.25", metric="spearman")
hmp_pathcov_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.1", metric="spearman")
hmp_pathcov_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, tab2 = hmp_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.05", metric="spearman")

# Calculate accuracy metrics:
hmp_picrust2_ec_nsti2_gg_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, category="NSTI=2 (GG)")
hmp_picrust2_ec_nsti2_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, category="NSTI=2")
hmp_picrust2_ec_nsti1.5_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, category="NSTI=1.5")
hmp_picrust2_ec_nsti1_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, category="NSTI=1")
hmp_picrust2_ec_nsti0.5_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, category="NSTI=0.5")
hmp_picrust2_ec_nsti0.25_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, category="NSTI=0.25")
hmp_picrust2_ec_nsti0.1_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, category="NSTI=0.1")
hmp_picrust2_ec_nsti0.05_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, category="NSTI=0.05")
hmp_paprica_ec_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_infiles$all_ecs_overlap$paprica_ec, category="PAPRICA")
hmp_null_ec_metrics <- calc_accuracy_metrics(hmp_ec_infiles$all_ecs_overlap$mgs_ec, hmp_ec_mgs_null_df_round, category="Null")

hmp_picrust2_pathabun_nsti2_gg_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, category="NSTI=2 (GG)")
hmp_picrust2_pathabun_nsti2_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, category="NSTI=2")
hmp_picrust2_pathabun_nsti1.5_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, category="NSTI=1.5")
hmp_picrust2_pathabun_nsti1_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, category="NSTI=1")
hmp_picrust2_pathabun_nsti0.5_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, category="NSTI=0.5")
hmp_picrust2_pathabun_nsti0.25_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, category="NSTI=0.25")
hmp_picrust2_pathabun_nsti0.1_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, category="NSTI=0.1")
hmp_picrust2_pathabun_nsti0.05_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, category="NSTI=0.05")
hmp_null_pathabun_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathabun$mgs_pathabun, hmp_pathabun_mgs_null_df_round, category="Null")

hmp_picrust2_pathcov_nsti2_gg_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, category="NSTI=2 (GG)")
hmp_picrust2_pathcov_nsti2_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, category="NSTI=2")
hmp_picrust2_pathcov_nsti1.5_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, category="NSTI=1.5")
hmp_picrust2_pathcov_nsti1_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, category="NSTI=1")
hmp_picrust2_pathcov_nsti0.5_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, category="NSTI=0.5")
hmp_picrust2_pathcov_nsti0.25_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, category="NSTI=0.25")
hmp_picrust2_pathcov_nsti0.1_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, category="NSTI=0.1")
hmp_picrust2_pathcov_nsti0.05_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, category="NSTI=0.05")
hmp_null_pathcov_metrics <- calc_accuracy_metrics(hmp_pathway_infiles$all_pathcov$mgs_pathcov, hmp_pathcov_mgs_null_df_round, category="Null")

# Create lists of spearman correlation and accuracy metrics.
hmp_ec_spearman_df <- rbind(hmp_ec_mgs_null, hmp_ec_paprica_vs_mgs, hmp_ec_picrust2_nsti2_gg_vs_mgs, hmp_ec_picrust2_nsti2_vs_mgs,
                            hmp_ec_picrust2_nsti1.5_vs_mgs, hmp_ec_picrust2_nsti1_vs_mgs, hmp_ec_picrust2_nsti0.5_vs_mgs,
                            hmp_ec_picrust2_nsti0.25_vs_mgs, hmp_ec_picrust2_nsti0.1_vs_mgs, hmp_ec_picrust2_nsti0.05_vs_mgs)

hmp_pathabun_spearman_df <- rbind(hmp_pathabun_mgs_null, hmp_pathabun_picrust2_nsti2_gg_vs_mgs, hmp_pathabun_picrust2_nsti2_vs_mgs,
                            hmp_pathabun_picrust2_nsti1.5_vs_mgs, hmp_pathabun_picrust2_nsti1_vs_mgs, hmp_pathabun_picrust2_nsti0.5_vs_mgs,
                            hmp_pathabun_picrust2_nsti0.25_vs_mgs, hmp_pathabun_picrust2_nsti0.1_vs_mgs, hmp_pathabun_picrust2_nsti0.05_vs_mgs)

hmp_pathcov_spearman_df <- rbind(hmp_pathcov_mgs_null, hmp_pathcov_picrust2_nsti2_gg_vs_mgs, hmp_pathcov_picrust2_nsti2_vs_mgs,
                            hmp_pathcov_picrust2_nsti1.5_vs_mgs, hmp_pathcov_picrust2_nsti1_vs_mgs, hmp_pathcov_picrust2_nsti0.5_vs_mgs,
                            hmp_pathcov_picrust2_nsti0.25_vs_mgs, hmp_pathcov_picrust2_nsti0.1_vs_mgs, hmp_pathcov_picrust2_nsti0.05_vs_mgs)

hmp_ec_acc <- rbind(hmp_null_ec_metrics,
                 hmp_paprica_ec_metrics,
                 hmp_picrust2_ec_nsti2_gg_metrics,
                 hmp_picrust2_ec_nsti2_metrics,
                 hmp_picrust2_ec_nsti1.5_metrics, 
                 hmp_picrust2_ec_nsti1_metrics,
                 hmp_picrust2_ec_nsti0.5_metrics,
                 hmp_picrust2_ec_nsti0.25_metrics,
                 hmp_picrust2_ec_nsti0.1_metrics,
                 hmp_picrust2_ec_nsti0.05_metrics)

hmp_pathabun_acc <- rbind(hmp_null_pathabun_metrics,
                    hmp_picrust2_pathabun_nsti2_gg_metrics,
                    hmp_picrust2_pathabun_nsti2_metrics,
                    hmp_picrust2_pathabun_nsti1.5_metrics, 
                    hmp_picrust2_pathabun_nsti1_metrics,
                    hmp_picrust2_pathabun_nsti0.5_metrics,
                    hmp_picrust2_pathabun_nsti0.25_metrics,
                    hmp_picrust2_pathabun_nsti0.1_metrics,
                    hmp_picrust2_pathabun_nsti0.05_metrics)

hmp_pathcov_acc <- rbind(hmp_null_pathcov_metrics,
                    hmp_picrust2_pathcov_nsti2_gg_metrics,
                    hmp_picrust2_pathcov_nsti2_metrics,
                    hmp_picrust2_pathcov_nsti1.5_metrics, 
                    hmp_picrust2_pathcov_nsti1_metrics,
                    hmp_picrust2_pathcov_nsti0.5_metrics,
                    hmp_picrust2_pathcov_nsti0.25_metrics,
                    hmp_picrust2_pathcov_nsti0.1_metrics,
                    hmp_picrust2_pathcov_nsti0.05_metrics)

# Get spearman correlations for each table with MGS:
mammal_ec_mgs_null <- cor_all_cols(tab1 = mammal_ec_mgs_null_df, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="Null", metric="spearman")
mammal_ec_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2 (GG)", metric="spearman")
mammal_ec_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2", metric="spearman")
mammal_ec_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1.5", metric="spearman")
mammal_ec_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1", metric="spearman")
mammal_ec_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.5", metric="spearman")
mammal_ec_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.25", metric="spearman")
mammal_ec_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.1", metric="spearman")
mammal_ec_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.05", metric="spearman")
mammal_ec_paprica_vs_mgs <- cor_all_cols(tab1 = mammal_ec_infiles$all_ecs_overlap$paprica_ec, tab2 = mammal_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="PAPRICA", metric="spearman")

mammal_pathabun_mgs_null <- cor_all_cols(tab1 = mammal_pathabun_mgs_null_df, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="Null", metric="spearman")
mammal_pathabun_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2 (GG)", metric="spearman")
mammal_pathabun_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2", metric="spearman")
mammal_pathabun_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1.5", metric="spearman")
mammal_pathabun_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1", metric="spearman")
mammal_pathabun_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.5", metric="spearman")
mammal_pathabun_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.25", metric="spearman")
mammal_pathabun_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.1", metric="spearman")
mammal_pathabun_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, tab2 = mammal_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.05", metric="spearman")

mammal_pathcov_mgs_null <- cor_all_cols(tab1 = mammal_pathcov_mgs_null_df, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="Null", metric="spearman")
mammal_pathcov_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2 (GG)", metric="spearman")
mammal_pathcov_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2", metric="spearman")
mammal_pathcov_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1.5", metric="spearman")
mammal_pathcov_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1", metric="spearman")
mammal_pathcov_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.5", metric="spearman")
mammal_pathcov_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.25", metric="spearman")
mammal_pathcov_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.1", metric="spearman")
mammal_pathcov_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, tab2 = mammal_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.05", metric="spearman")

# Calculate accuracy metrics:
mammal_picrust2_ec_nsti2_gg_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, category="NSTI=2 (GG)")
mammal_picrust2_ec_nsti2_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, category="NSTI=2")
mammal_picrust2_ec_nsti1.5_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, category="NSTI=1.5")
mammal_picrust2_ec_nsti1_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, category="NSTI=1")
mammal_picrust2_ec_nsti0.5_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, category="NSTI=0.5")
mammal_picrust2_ec_nsti0.25_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, category="NSTI=0.25")
mammal_picrust2_ec_nsti0.1_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, category="NSTI=0.1")
mammal_picrust2_ec_nsti0.05_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, category="NSTI=0.05")
mammal_paprica_ec_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_infiles$all_ecs_overlap$paprica_ec, category="PAPRICA")
mammal_null_ec_metrics <- calc_accuracy_metrics(mammal_ec_infiles$all_ecs_overlap$mgs_ec, mammal_ec_mgs_null_df_round, category="Null")

mammal_picrust2_pathabun_nsti2_gg_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, category="NSTI=2 (GG)")
mammal_picrust2_pathabun_nsti2_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, category="NSTI=2")
mammal_picrust2_pathabun_nsti1.5_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, category="NSTI=1.5")
mammal_picrust2_pathabun_nsti1_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, category="NSTI=1")
mammal_picrust2_pathabun_nsti0.5_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, category="NSTI=0.5")
mammal_picrust2_pathabun_nsti0.25_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, category="NSTI=0.25")
mammal_picrust2_pathabun_nsti0.1_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, category="NSTI=0.1")
mammal_picrust2_pathabun_nsti0.05_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, category="NSTI=0.05")
mammal_null_pathabun_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathabun$mgs_pathabun, mammal_pathabun_mgs_null_df_round, category="Null")

mammal_picrust2_pathcov_nsti2_gg_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, category="NSTI=2 (GG)")
mammal_picrust2_pathcov_nsti2_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, category="NSTI=2")
mammal_picrust2_pathcov_nsti1.5_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, category="NSTI=1.5")
mammal_picrust2_pathcov_nsti1_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, category="NSTI=1")
mammal_picrust2_pathcov_nsti0.5_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, category="NSTI=0.5")
mammal_picrust2_pathcov_nsti0.25_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, category="NSTI=0.25")
mammal_picrust2_pathcov_nsti0.1_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, category="NSTI=0.1")
mammal_picrust2_pathcov_nsti0.05_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, category="NSTI=0.05")
mammal_null_pathcov_metrics <- calc_accuracy_metrics(mammal_pathway_infiles$all_pathcov$mgs_pathcov, mammal_pathcov_mgs_null_df_round, category="Null")

# Create lists of spearman correlation and accuracy metrics.
mammal_ec_spearman_df <- rbind(mammal_ec_mgs_null, mammal_ec_paprica_vs_mgs, mammal_ec_picrust2_nsti2_gg_vs_mgs, mammal_ec_picrust2_nsti2_vs_mgs,
                            mammal_ec_picrust2_nsti1.5_vs_mgs, mammal_ec_picrust2_nsti1_vs_mgs, mammal_ec_picrust2_nsti0.5_vs_mgs,
                            mammal_ec_picrust2_nsti0.25_vs_mgs, mammal_ec_picrust2_nsti0.1_vs_mgs, mammal_ec_picrust2_nsti0.05_vs_mgs)

mammal_pathabun_spearman_df <- rbind(mammal_pathabun_mgs_null, mammal_pathabun_picrust2_nsti2_gg_vs_mgs, mammal_pathabun_picrust2_nsti2_vs_mgs,
                                  mammal_pathabun_picrust2_nsti1.5_vs_mgs, mammal_pathabun_picrust2_nsti1_vs_mgs, mammal_pathabun_picrust2_nsti0.5_vs_mgs,
                                  mammal_pathabun_picrust2_nsti0.25_vs_mgs, mammal_pathabun_picrust2_nsti0.1_vs_mgs, mammal_pathabun_picrust2_nsti0.05_vs_mgs)

mammal_pathcov_spearman_df <- rbind(mammal_pathcov_mgs_null, mammal_pathcov_picrust2_nsti2_gg_vs_mgs, mammal_pathcov_picrust2_nsti2_vs_mgs,
                                 mammal_pathcov_picrust2_nsti1.5_vs_mgs, mammal_pathcov_picrust2_nsti1_vs_mgs, mammal_pathcov_picrust2_nsti0.5_vs_mgs,
                                 mammal_pathcov_picrust2_nsti0.25_vs_mgs, mammal_pathcov_picrust2_nsti0.1_vs_mgs, mammal_pathcov_picrust2_nsti0.05_vs_mgs)

mammal_ec_acc <- rbind(mammal_null_ec_metrics,
                    mammal_paprica_ec_metrics,
                    mammal_picrust2_ec_nsti2_gg_metrics,
                    mammal_picrust2_ec_nsti2_metrics,
                    mammal_picrust2_ec_nsti1.5_metrics, 
                    mammal_picrust2_ec_nsti1_metrics,
                    mammal_picrust2_ec_nsti0.5_metrics,
                    mammal_picrust2_ec_nsti0.25_metrics,
                    mammal_picrust2_ec_nsti0.1_metrics,
                    mammal_picrust2_ec_nsti0.05_metrics)

mammal_pathabun_acc <- rbind(mammal_null_pathabun_metrics,
                          mammal_picrust2_pathabun_nsti2_gg_metrics,
                          mammal_picrust2_pathabun_nsti2_metrics,
                          mammal_picrust2_pathabun_nsti1.5_metrics, 
                          mammal_picrust2_pathabun_nsti1_metrics,
                          mammal_picrust2_pathabun_nsti0.5_metrics,
                          mammal_picrust2_pathabun_nsti0.25_metrics,
                          mammal_picrust2_pathabun_nsti0.1_metrics,
                          mammal_picrust2_pathabun_nsti0.05_metrics)

mammal_pathcov_acc <- rbind(mammal_null_pathcov_metrics,
                         mammal_picrust2_pathcov_nsti2_gg_metrics,
                         mammal_picrust2_pathcov_nsti2_metrics,
                         mammal_picrust2_pathcov_nsti1.5_metrics, 
                         mammal_picrust2_pathcov_nsti1_metrics,
                         mammal_picrust2_pathcov_nsti0.5_metrics,
                         mammal_picrust2_pathcov_nsti0.25_metrics,
                         mammal_picrust2_pathcov_nsti0.1_metrics,
                         mammal_picrust2_pathcov_nsti0.05_metrics)



# Ocean
# Get spearman correlations for each table with MGS:
ocean_ec_mgs_null <- cor_all_cols(tab1 = ocean_ec_mgs_null_df, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="Null", metric="spearman")
ocean_ec_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2 (GG)", metric="spearman")
ocean_ec_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2", metric="spearman")
ocean_ec_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1.5", metric="spearman")
ocean_ec_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1", metric="spearman")
ocean_ec_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.5", metric="spearman")
ocean_ec_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.25", metric="spearman")
ocean_ec_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.1", metric="spearman")
ocean_ec_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.05", metric="spearman")
ocean_ec_paprica_vs_mgs <- cor_all_cols(tab1 = ocean_ec_infiles$all_ecs_overlap$paprica_ec, tab2 = ocean_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="PAPRICA", metric="spearman")

ocean_pathabun_mgs_null <- cor_all_cols(tab1 = ocean_pathabun_mgs_null_df, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="Null", metric="spearman")
ocean_pathabun_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2 (GG)", metric="spearman")
ocean_pathabun_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2", metric="spearman")
ocean_pathabun_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1.5", metric="spearman")
ocean_pathabun_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1", metric="spearman")
ocean_pathabun_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.5", metric="spearman")
ocean_pathabun_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.25", metric="spearman")
ocean_pathabun_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.1", metric="spearman")
ocean_pathabun_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, tab2 = ocean_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.05", metric="spearman")

ocean_pathcov_mgs_null <- cor_all_cols(tab1 = ocean_pathcov_mgs_null_df, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="Null", metric="spearman")
ocean_pathcov_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2 (GG)", metric="spearman")
ocean_pathcov_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2", metric="spearman")
ocean_pathcov_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1.5", metric="spearman")
ocean_pathcov_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1", metric="spearman")
ocean_pathcov_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.5", metric="spearman")
ocean_pathcov_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.25", metric="spearman")
ocean_pathcov_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.1", metric="spearman")
ocean_pathcov_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, tab2 = ocean_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.05", metric="spearman")

# Calculate accuracy metrics:
ocean_picrust2_ec_nsti2_gg_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, category="NSTI=2 (GG)")
ocean_picrust2_ec_nsti2_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, category="NSTI=2")
ocean_picrust2_ec_nsti1.5_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, category="NSTI=1.5")
ocean_picrust2_ec_nsti1_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, category="NSTI=1")
ocean_picrust2_ec_nsti0.5_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, category="NSTI=0.5")
ocean_picrust2_ec_nsti0.25_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, category="NSTI=0.25")
ocean_picrust2_ec_nsti0.1_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, category="NSTI=0.1")
ocean_picrust2_ec_nsti0.05_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, category="NSTI=0.05")
ocean_paprica_ec_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_infiles$all_ecs_overlap$paprica_ec, category="PAPRICA")
ocean_null_ec_metrics <- calc_accuracy_metrics(ocean_ec_infiles$all_ecs_overlap$mgs_ec, ocean_ec_mgs_null_df_round, category="Null")

ocean_picrust2_pathabun_nsti2_gg_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, category="NSTI=2 (GG)")
ocean_picrust2_pathabun_nsti2_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, category="NSTI=2")
ocean_picrust2_pathabun_nsti1.5_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, category="NSTI=1.5")
ocean_picrust2_pathabun_nsti1_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, category="NSTI=1")
ocean_picrust2_pathabun_nsti0.5_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, category="NSTI=0.5")
ocean_picrust2_pathabun_nsti0.25_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, category="NSTI=0.25")
ocean_picrust2_pathabun_nsti0.1_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, category="NSTI=0.1")
ocean_picrust2_pathabun_nsti0.05_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, category="NSTI=0.05")
ocean_null_pathabun_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathabun$mgs_pathabun, ocean_pathabun_mgs_null_df_round, category="Null")

ocean_picrust2_pathcov_nsti2_gg_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, category="NSTI=2 (GG)")
ocean_picrust2_pathcov_nsti2_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, category="NSTI=2")
ocean_picrust2_pathcov_nsti1.5_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, category="NSTI=1.5")
ocean_picrust2_pathcov_nsti1_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, category="NSTI=1")
ocean_picrust2_pathcov_nsti0.5_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, category="NSTI=0.5")
ocean_picrust2_pathcov_nsti0.25_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, category="NSTI=0.25")
ocean_picrust2_pathcov_nsti0.1_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, category="NSTI=0.1")
ocean_picrust2_pathcov_nsti0.05_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, category="NSTI=0.05")
ocean_null_pathcov_metrics <- calc_accuracy_metrics(ocean_pathway_infiles$all_pathcov$mgs_pathcov, ocean_pathcov_mgs_null_df_round, category="Null")

# Create lists of spearman correlation and accuracy metrics.
ocean_ec_spearman_df <- rbind(ocean_ec_mgs_null, ocean_ec_paprica_vs_mgs, ocean_ec_picrust2_nsti2_gg_vs_mgs, ocean_ec_picrust2_nsti2_vs_mgs,
                            ocean_ec_picrust2_nsti1.5_vs_mgs, ocean_ec_picrust2_nsti1_vs_mgs, ocean_ec_picrust2_nsti0.5_vs_mgs,
                            ocean_ec_picrust2_nsti0.25_vs_mgs, ocean_ec_picrust2_nsti0.1_vs_mgs, ocean_ec_picrust2_nsti0.05_vs_mgs)

ocean_pathabun_spearman_df <- rbind(ocean_pathabun_mgs_null, ocean_pathabun_picrust2_nsti2_gg_vs_mgs, ocean_pathabun_picrust2_nsti2_vs_mgs,
                                  ocean_pathabun_picrust2_nsti1.5_vs_mgs, ocean_pathabun_picrust2_nsti1_vs_mgs, ocean_pathabun_picrust2_nsti0.5_vs_mgs,
                                  ocean_pathabun_picrust2_nsti0.25_vs_mgs, ocean_pathabun_picrust2_nsti0.1_vs_mgs, ocean_pathabun_picrust2_nsti0.05_vs_mgs)

ocean_pathcov_spearman_df <- rbind(ocean_pathcov_mgs_null, ocean_pathcov_picrust2_nsti2_gg_vs_mgs, ocean_pathcov_picrust2_nsti2_vs_mgs,
                                 ocean_pathcov_picrust2_nsti1.5_vs_mgs, ocean_pathcov_picrust2_nsti1_vs_mgs, ocean_pathcov_picrust2_nsti0.5_vs_mgs,
                                 ocean_pathcov_picrust2_nsti0.25_vs_mgs, ocean_pathcov_picrust2_nsti0.1_vs_mgs, ocean_pathcov_picrust2_nsti0.05_vs_mgs)

ocean_ec_acc <- rbind(ocean_null_ec_metrics,
                    ocean_paprica_ec_metrics,
                    ocean_picrust2_ec_nsti2_gg_metrics,
                    ocean_picrust2_ec_nsti2_metrics,
                    ocean_picrust2_ec_nsti1.5_metrics, 
                    ocean_picrust2_ec_nsti1_metrics,
                    ocean_picrust2_ec_nsti0.5_metrics,
                    ocean_picrust2_ec_nsti0.25_metrics,
                    ocean_picrust2_ec_nsti0.1_metrics,
                    ocean_picrust2_ec_nsti0.05_metrics)

ocean_pathabun_acc <- rbind(ocean_null_pathabun_metrics,
                          ocean_picrust2_pathabun_nsti2_gg_metrics,
                          ocean_picrust2_pathabun_nsti2_metrics,
                          ocean_picrust2_pathabun_nsti1.5_metrics, 
                          ocean_picrust2_pathabun_nsti1_metrics,
                          ocean_picrust2_pathabun_nsti0.5_metrics,
                          ocean_picrust2_pathabun_nsti0.25_metrics,
                          ocean_picrust2_pathabun_nsti0.1_metrics,
                          ocean_picrust2_pathabun_nsti0.05_metrics)

ocean_pathcov_acc <- rbind(ocean_null_pathcov_metrics,
                         ocean_picrust2_pathcov_nsti2_gg_metrics,
                         ocean_picrust2_pathcov_nsti2_metrics,
                         ocean_picrust2_pathcov_nsti1.5_metrics, 
                         ocean_picrust2_pathcov_nsti1_metrics,
                         ocean_picrust2_pathcov_nsti0.5_metrics,
                         ocean_picrust2_pathcov_nsti0.25_metrics,
                         ocean_picrust2_pathcov_nsti0.1_metrics,
                         ocean_picrust2_pathcov_nsti0.05_metrics)


# blueberry
# Get spearman correlations for each table with MGS:
blueberry_ec_mgs_null <- cor_all_cols(tab1 = blueberry_ec_mgs_null_df, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="Null", metric="spearman")
blueberry_ec_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2 (GG)", metric="spearman")
blueberry_ec_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=2", metric="spearman")
blueberry_ec_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1.5", metric="spearman")
blueberry_ec_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=1", metric="spearman")
blueberry_ec_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.5", metric="spearman")
blueberry_ec_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.25", metric="spearman")
blueberry_ec_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.1", metric="spearman")
blueberry_ec_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="NSTI=0.05", metric="spearman")
blueberry_ec_paprica_vs_mgs <- cor_all_cols(tab1 = blueberry_ec_infiles$all_ecs_overlap$paprica_ec, tab2 = blueberry_ec_infiles$all_ecs_overlap$mgs_ec, cat_string="PAPRICA", metric="spearman")

blueberry_pathabun_mgs_null <- cor_all_cols(tab1 = blueberry_pathabun_mgs_null_df, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="Null", metric="spearman")
blueberry_pathabun_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2 (GG)", metric="spearman")
blueberry_pathabun_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=2", metric="spearman")
blueberry_pathabun_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1.5", metric="spearman")
blueberry_pathabun_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=1", metric="spearman")
blueberry_pathabun_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.5", metric="spearman")
blueberry_pathabun_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.25", metric="spearman")
blueberry_pathabun_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.1", metric="spearman")
blueberry_pathabun_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, tab2 = blueberry_pathway_infiles$all_pathabun$mgs_pathabun, cat_string="NSTI=0.05", metric="spearman")

blueberry_pathcov_mgs_null <- cor_all_cols(tab1 = blueberry_pathcov_mgs_null_df, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="Null", metric="spearman")
blueberry_pathcov_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2 (GG)", metric="spearman")
blueberry_pathcov_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=2", metric="spearman")
blueberry_pathcov_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1.5", metric="spearman")
blueberry_pathcov_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=1", metric="spearman")
blueberry_pathcov_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.5", metric="spearman")
blueberry_pathcov_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.25", metric="spearman")
blueberry_pathcov_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.1", metric="spearman")
blueberry_pathcov_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, tab2 = blueberry_pathway_infiles$all_pathcov$mgs_pathcov, cat_string="NSTI=0.05", metric="spearman")

# Calculate accuracy metrics:
blueberry_picrust2_ec_nsti2_gg_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2_gg, category="NSTI=2 (GG)")
blueberry_picrust2_ec_nsti2_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti2, category="NSTI=2")
blueberry_picrust2_ec_nsti1.5_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1.5, category="NSTI=1.5")
blueberry_picrust2_ec_nsti1_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti1, category="NSTI=1")
blueberry_picrust2_ec_nsti0.5_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.5, category="NSTI=0.5")
blueberry_picrust2_ec_nsti0.25_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.25, category="NSTI=0.25")
blueberry_picrust2_ec_nsti0.1_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.1, category="NSTI=0.1")
blueberry_picrust2_ec_nsti0.05_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$picrust2_ec_nsti0.05, category="NSTI=0.05")
blueberry_paprica_ec_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_infiles$all_ecs_overlap$paprica_ec, category="PAPRICA")
blueberry_null_ec_metrics <- calc_accuracy_metrics(blueberry_ec_infiles$all_ecs_overlap$mgs_ec, blueberry_ec_mgs_null_df_round, category="Null")

blueberry_picrust2_pathabun_nsti2_gg_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2_gg, category="NSTI=2 (GG)")
blueberry_picrust2_pathabun_nsti2_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti2, category="NSTI=2")
blueberry_picrust2_pathabun_nsti1.5_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1.5, category="NSTI=1.5")
blueberry_picrust2_pathabun_nsti1_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti1, category="NSTI=1")
blueberry_picrust2_pathabun_nsti0.5_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.5, category="NSTI=0.5")
blueberry_picrust2_pathabun_nsti0.25_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.25, category="NSTI=0.25")
blueberry_picrust2_pathabun_nsti0.1_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.1, category="NSTI=0.1")
blueberry_picrust2_pathabun_nsti0.05_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathway_infiles$all_pathabun$picrust2_pathabun_nsti0.05, category="NSTI=0.05")
blueberry_null_pathabun_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathabun$mgs_pathabun, blueberry_pathabun_mgs_null_df_round, category="Null")

blueberry_picrust2_pathcov_nsti2_gg_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2_gg, category="NSTI=2 (GG)")
blueberry_picrust2_pathcov_nsti2_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti2, category="NSTI=2")
blueberry_picrust2_pathcov_nsti1.5_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1.5, category="NSTI=1.5")
blueberry_picrust2_pathcov_nsti1_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti1, category="NSTI=1")
blueberry_picrust2_pathcov_nsti0.5_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.5, category="NSTI=0.5")
blueberry_picrust2_pathcov_nsti0.25_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.25, category="NSTI=0.25")
blueberry_picrust2_pathcov_nsti0.1_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.1, category="NSTI=0.1")
blueberry_picrust2_pathcov_nsti0.05_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathway_infiles$all_pathcov$picrust2_pathcov_nsti0.05, category="NSTI=0.05")
blueberry_null_pathcov_metrics <- calc_accuracy_metrics(blueberry_pathway_infiles$all_pathcov$mgs_pathcov, blueberry_pathcov_mgs_null_df_round, category="Null")

# Create lists of spearman correlation and accuracy metrics.
blueberry_ec_spearman_df <- rbind(blueberry_ec_mgs_null, blueberry_ec_paprica_vs_mgs, blueberry_ec_picrust2_nsti2_gg_vs_mgs, blueberry_ec_picrust2_nsti2_vs_mgs,
                            blueberry_ec_picrust2_nsti1.5_vs_mgs, blueberry_ec_picrust2_nsti1_vs_mgs, blueberry_ec_picrust2_nsti0.5_vs_mgs,
                            blueberry_ec_picrust2_nsti0.25_vs_mgs, blueberry_ec_picrust2_nsti0.1_vs_mgs, blueberry_ec_picrust2_nsti0.05_vs_mgs)

blueberry_pathabun_spearman_df <- rbind(blueberry_pathabun_mgs_null, blueberry_pathabun_picrust2_nsti2_gg_vs_mgs, blueberry_pathabun_picrust2_nsti2_vs_mgs,
                                  blueberry_pathabun_picrust2_nsti1.5_vs_mgs, blueberry_pathabun_picrust2_nsti1_vs_mgs, blueberry_pathabun_picrust2_nsti0.5_vs_mgs,
                                  blueberry_pathabun_picrust2_nsti0.25_vs_mgs, blueberry_pathabun_picrust2_nsti0.1_vs_mgs, blueberry_pathabun_picrust2_nsti0.05_vs_mgs)

blueberry_pathcov_spearman_df <- rbind(blueberry_pathcov_mgs_null, blueberry_pathcov_picrust2_nsti2_gg_vs_mgs, blueberry_pathcov_picrust2_nsti2_vs_mgs,
                                 blueberry_pathcov_picrust2_nsti1.5_vs_mgs, blueberry_pathcov_picrust2_nsti1_vs_mgs, blueberry_pathcov_picrust2_nsti0.5_vs_mgs,
                                 blueberry_pathcov_picrust2_nsti0.25_vs_mgs, blueberry_pathcov_picrust2_nsti0.1_vs_mgs, blueberry_pathcov_picrust2_nsti0.05_vs_mgs)

blueberry_ec_acc <- rbind(blueberry_null_ec_metrics,
                    blueberry_paprica_ec_metrics,
                    blueberry_picrust2_ec_nsti2_gg_metrics,
                    blueberry_picrust2_ec_nsti2_metrics,
                    blueberry_picrust2_ec_nsti1.5_metrics, 
                    blueberry_picrust2_ec_nsti1_metrics,
                    blueberry_picrust2_ec_nsti0.5_metrics,
                    blueberry_picrust2_ec_nsti0.25_metrics,
                    blueberry_picrust2_ec_nsti0.1_metrics,
                    blueberry_picrust2_ec_nsti0.05_metrics)

blueberry_pathabun_acc <- rbind(blueberry_null_pathabun_metrics,
                          blueberry_picrust2_pathabun_nsti2_gg_metrics,
                          blueberry_picrust2_pathabun_nsti2_metrics,
                          blueberry_picrust2_pathabun_nsti1.5_metrics, 
                          blueberry_picrust2_pathabun_nsti1_metrics,
                          blueberry_picrust2_pathabun_nsti0.5_metrics,
                          blueberry_picrust2_pathabun_nsti0.25_metrics,
                          blueberry_picrust2_pathabun_nsti0.1_metrics,
                          blueberry_picrust2_pathabun_nsti0.05_metrics)

blueberry_pathcov_acc <- rbind(blueberry_null_pathcov_metrics,
                         blueberry_picrust2_pathcov_nsti2_gg_metrics,
                         blueberry_picrust2_pathcov_nsti2_metrics,
                         blueberry_picrust2_pathcov_nsti1.5_metrics, 
                         blueberry_picrust2_pathcov_nsti1_metrics,
                         blueberry_picrust2_pathcov_nsti0.5_metrics,
                         blueberry_picrust2_pathcov_nsti0.25_metrics,
                         blueberry_picrust2_pathcov_nsti0.1_metrics,
                         blueberry_picrust2_pathcov_nsti0.05_metrics)




# Save RDS objects:
saveRDS(object = hmp_ec_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_ec_spearman_df.rds")
saveRDS(object = hmp_ec_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_ec_acc_metrics.rds")
saveRDS(object = hmp_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_pathabun_spearman_df.rds")
saveRDS(object = hmp_pathabun_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_pathabun_acc_metrics.rds")
saveRDS(object = hmp_pathcov_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_pathcov_spearman_df.rds")
saveRDS(object = hmp_pathcov_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/hmp_pathcov_acc_metrics.rds")

saveRDS(object = mammal_ec_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_ec_spearman_df.rds")
saveRDS(object = mammal_ec_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_ec_acc_metrics.rds")
saveRDS(object = mammal_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_pathabun_spearman_df.rds")
saveRDS(object = mammal_pathabun_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_pathabun_acc_metrics.rds")
saveRDS(object = mammal_pathcov_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_pathcov_spearman_df.rds")
saveRDS(object = mammal_pathcov_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/mammal_pathcov_acc_metrics.rds")

saveRDS(object = ocean_ec_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_ec_spearman_df.rds")
saveRDS(object = ocean_ec_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_ec_acc_metrics.rds")
saveRDS(object = ocean_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_pathabun_spearman_df.rds")
saveRDS(object = ocean_pathabun_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_pathabun_acc_metrics.rds")
saveRDS(object = ocean_pathcov_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_pathcov_spearman_df.rds")
saveRDS(object = ocean_pathcov_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/ocean_pathcov_acc_metrics.rds")

saveRDS(object = blueberry_ec_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_ec_spearman_df.rds")
saveRDS(object = blueberry_ec_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_ec_acc_metrics.rds")
saveRDS(object = blueberry_pathabun_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_pathabun_spearman_df.rds")
saveRDS(object = blueberry_pathabun_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_pathabun_acc_metrics.rds")
saveRDS(object = blueberry_pathcov_spearman_df, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_pathcov_spearman_df.rds")
saveRDS(object = blueberry_pathcov_acc, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_metrics/blueberry_pathcov_acc_metrics.rds")
