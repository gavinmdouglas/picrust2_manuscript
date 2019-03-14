### Comparing predicted functions based on 16S to MGS "gold standard".
### Comparisons made between KEGG orthologs predicted by each tool to MGS.
### RDS files saved for spearman correlations and accuracy metrics.

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

### Read in database files (used for calculating null distribution).
ko <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)

# Read in all tables, restrict to overlapping samples only, and get subsets with all possible KOs filled in and with
# only KOs overlapping between all tools. Focus analyses on KOs overlapping between all tools, but outputted others as well
# for sanity checks.
hmp_infiles <- read_in_ko_predictions("hmp")
mammal_infiles <- read_in_ko_predictions("mammal")
ocean_infiles <- read_in_ko_predictions("ocean")
blueberry_infiles <- read_in_ko_predictions("blueberry")

# Generate random tables based on subsampling database to calculate null distributions.
hmp_ko_mgs_null_df <- generate_null_mean_db_funcs(db = ko, tab = hmp_infiles$all_kos_overlap$mgs_ko)
mammal_ko_mgs_null_df <- generate_null_mean_db_funcs(db = ko, tab = mammal_infiles$all_kos_overlap$mgs_ko)
ocean_ko_mgs_null_df <- generate_null_mean_db_funcs(db = ko, tab = ocean_infiles$all_kos_overlap$mgs_ko)
blueberry_ko_mgs_null_df <- generate_null_mean_db_funcs(db = ko, tab = blueberry_infiles$all_kos_overlap$mgs_ko)

hmp_ko_mgs_null_df_round <- round(hmp_ko_mgs_null_df - 0.00000001)
mammal_ko_mgs_null_df_round <- round(mammal_ko_mgs_null_df - 0.00000001)
ocean_ko_mgs_null_df_round <- round(ocean_ko_mgs_null_df - 0.00000001)
blueberry_ko_mgs_null_df_round <- round(blueberry_ko_mgs_null_df - 0.00000001)

# Get spearman correlations for each table with MGS:
hmp_ko_mgs_null <- cor_all_cols(tab1 = hmp_ko_mgs_null_df, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="Null", metric="spearman")
hmp_ko_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2 (GG)", metric="spearman")
hmp_ko_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti2, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2", metric="spearman")
hmp_ko_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti1.5, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1.5", metric="spearman")
hmp_ko_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti1, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1", metric="spearman")
hmp_ko_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.5, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.5", metric="spearman")
hmp_ko_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.25, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.25", metric="spearman")
hmp_ko_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.1, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.1", metric="spearman")
hmp_ko_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.05, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.05", metric="spearman")
hmp_ko_picrust1_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$picrust1_ko, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="PICRUSt1", metric="spearman")
hmp_ko_panfp_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$panfp_ko, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="PanFP", metric="spearman")
hmp_ko_piphillin_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$piphillin_ko, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="Piphillin", metric="spearman")
hmp_ko_tax4fun_vs_mgs <- cor_all_cols(tab1 = hmp_infiles$all_kos_overlap$tax4fun_ko, tab2 = hmp_infiles$all_kos_overlap$mgs_ko, cat_string="Tax4Fun", metric="spearman")

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


ocean_ko_mgs_null <- cor_all_cols(tab1 = ocean_ko_mgs_null_df, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="Null", metric="spearman")
ocean_ko_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2 (GG)", metric="spearman")
ocean_ko_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti2, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2", metric="spearman")
ocean_ko_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti1.5, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1.5", metric="spearman")
ocean_ko_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti1, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1", metric="spearman")
ocean_ko_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.5, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.5", metric="spearman")
ocean_ko_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.25, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.25", metric="spearman")
ocean_ko_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.1, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.1", metric="spearman")
ocean_ko_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.05, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.05", metric="spearman")
ocean_ko_picrust1_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$picrust1_ko, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="PICRUSt1", metric="spearman")
ocean_ko_panfp_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$panfp_ko, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="PanFP", metric="spearman")
ocean_ko_piphillin_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$piphillin_ko, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="Piphillin", metric="spearman")
ocean_ko_tax4fun_vs_mgs <- cor_all_cols(tab1 = ocean_infiles$all_kos_overlap$tax4fun_ko, tab2 = ocean_infiles$all_kos_overlap$mgs_ko, cat_string="Tax4Fun", metric="spearman")


blueberry_ko_mgs_null <- cor_all_cols(tab1 = blueberry_ko_mgs_null_df, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="Null", metric="spearman")
blueberry_ko_picrust2_nsti2_gg_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2 (GG)", metric="spearman")
blueberry_ko_picrust2_nsti2_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti2, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=2", metric="spearman")
blueberry_ko_picrust2_nsti1.5_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti1.5, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1.5", metric="spearman")
blueberry_ko_picrust2_nsti1_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti1, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=1", metric="spearman")
blueberry_ko_picrust2_nsti0.5_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.5, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.5", metric="spearman")
blueberry_ko_picrust2_nsti0.25_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.25, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.25", metric="spearman")
blueberry_ko_picrust2_nsti0.1_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.1, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.1", metric="spearman")
blueberry_ko_picrust2_nsti0.05_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.05, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="NSTI=0.05", metric="spearman")
blueberry_ko_picrust1_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$picrust1_ko, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="PICRUSt1", metric="spearman")
blueberry_ko_panfp_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$panfp_ko, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="PanFP", metric="spearman")
blueberry_ko_piphillin_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$piphillin_ko, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="Piphillin", metric="spearman")
blueberry_ko_tax4fun_vs_mgs <- cor_all_cols(tab1 = blueberry_infiles$all_kos_overlap$tax4fun_ko, tab2 = blueberry_infiles$all_kos_overlap$mgs_ko, cat_string="Tax4Fun", metric="spearman")

# Make combined dfs of spearman correlation coefficient subsets of interest:
hmp_ko_spearman_df <- rbind(hmp_ko_mgs_null, hmp_ko_tax4fun_vs_mgs, hmp_ko_panfp_vs_mgs, hmp_ko_piphillin_vs_mgs,
                            hmp_ko_picrust1_vs_mgs, hmp_ko_picrust2_nsti2_gg_vs_mgs, hmp_ko_picrust2_nsti2_vs_mgs,
                            hmp_ko_picrust2_nsti1.5_vs_mgs, hmp_ko_picrust2_nsti1_vs_mgs, hmp_ko_picrust2_nsti0.5_vs_mgs,
                            hmp_ko_picrust2_nsti0.25_vs_mgs, hmp_ko_picrust2_nsti0.1_vs_mgs, hmp_ko_picrust2_nsti0.05_vs_mgs)

mammal_ko_spearman_df <- rbind(mammal_ko_mgs_null, mammal_ko_tax4fun_vs_mgs, mammal_ko_panfp_vs_mgs, mammal_ko_piphillin_vs_mgs,
                            mammal_ko_picrust1_vs_mgs, mammal_ko_picrust2_nsti2_gg_vs_mgs, mammal_ko_picrust2_nsti2_vs_mgs,
                            mammal_ko_picrust2_nsti1.5_vs_mgs, mammal_ko_picrust2_nsti1_vs_mgs, mammal_ko_picrust2_nsti0.5_vs_mgs,
                            mammal_ko_picrust2_nsti0.25_vs_mgs, mammal_ko_picrust2_nsti0.1_vs_mgs, mammal_ko_picrust2_nsti0.05_vs_mgs)

ocean_ko_spearman_df <- rbind(ocean_ko_mgs_null, ocean_ko_tax4fun_vs_mgs, ocean_ko_panfp_vs_mgs, ocean_ko_piphillin_vs_mgs,
                            ocean_ko_picrust1_vs_mgs, ocean_ko_picrust2_nsti2_gg_vs_mgs, ocean_ko_picrust2_nsti2_vs_mgs,
                            ocean_ko_picrust2_nsti1.5_vs_mgs, ocean_ko_picrust2_nsti1_vs_mgs, ocean_ko_picrust2_nsti0.5_vs_mgs,
                            ocean_ko_picrust2_nsti0.25_vs_mgs, ocean_ko_picrust2_nsti0.1_vs_mgs, ocean_ko_picrust2_nsti0.05_vs_mgs)


blueberry_ko_spearman_df <- rbind(blueberry_ko_mgs_null, blueberry_ko_tax4fun_vs_mgs, blueberry_ko_panfp_vs_mgs, blueberry_ko_piphillin_vs_mgs,
                            blueberry_ko_picrust1_vs_mgs, blueberry_ko_picrust2_nsti2_gg_vs_mgs, blueberry_ko_picrust2_nsti2_vs_mgs,
                            blueberry_ko_picrust2_nsti1.5_vs_mgs, blueberry_ko_picrust2_nsti1_vs_mgs, blueberry_ko_picrust2_nsti0.5_vs_mgs,
                            blueberry_ko_picrust2_nsti0.25_vs_mgs, blueberry_ko_picrust2_nsti0.1_vs_mgs, blueberry_ko_picrust2_nsti0.05_vs_mgs)


hmp_picrust2_ko_nsti2_gg_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, category="NSTI=2 (GG)")
hmp_picrust2_ko_nsti2_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti2, category="NSTI=2")
hmp_picrust2_ko_nsti1.5_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti1.5, category="NSTI=1.5")
hmp_picrust2_ko_nsti1_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti1, category="NSTI=1")
hmp_picrust2_ko_nsti0.5_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.5, category="NSTI=0.5")
hmp_picrust2_ko_nsti0.25_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.25, category="NSTI=0.25")
hmp_picrust2_ko_nsti0.1_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.1, category="NSTI=0.1")
hmp_picrust2_ko_nsti0.05_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust2_ko_nsti0.05, category="NSTI=0.05")
hmp_picrust1_ko_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$picrust1_ko, category="PICRUSt1")
hmp_piphillin_ko_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$piphillin, category="Piphillin")
hmp_panfp_ko_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$panfp, category="PanFP")
hmp_tax4fun_ko_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_infiles$all_kos_overlap$tax4fun, category="Tax4Fun")
hmp_null_ko_metrics <- calc_accuracy_metrics(hmp_infiles$all_kos_overlap$mgs_ko, hmp_ko_mgs_null_df_round, category="Null")

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


ocean_picrust2_ko_nsti2_gg_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, category="NSTI=2 (GG)")
ocean_picrust2_ko_nsti2_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti2, category="NSTI=2")
ocean_picrust2_ko_nsti1.5_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti1.5, category="NSTI=1.5")
ocean_picrust2_ko_nsti1_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti1, category="NSTI=1")
ocean_picrust2_ko_nsti0.5_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.5, category="NSTI=0.5")
ocean_picrust2_ko_nsti0.25_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.25, category="NSTI=0.25")
ocean_picrust2_ko_nsti0.1_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.1, category="NSTI=0.1")
ocean_picrust2_ko_nsti0.05_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust2_ko_nsti0.05, category="NSTI=0.05")
ocean_picrust1_ko_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$picrust1_ko, category="PICRUSt1")
ocean_piphillin_ko_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$piphillin, category="Piphillin")
ocean_panfp_ko_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$panfp, category="PanFP")
ocean_tax4fun_ko_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_infiles$all_kos_overlap$tax4fun, category="Tax4Fun")
ocean_null_ko_metrics <- calc_accuracy_metrics(ocean_infiles$all_kos_overlap$mgs_ko, ocean_ko_mgs_null_df_round, category="Null")


blueberry_picrust2_ko_nsti2_gg_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti2_gg, category="NSTI=2 (GG)")
blueberry_picrust2_ko_nsti2_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti2, category="NSTI=2")
blueberry_picrust2_ko_nsti1.5_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti1.5, category="NSTI=1.5")
blueberry_picrust2_ko_nsti1_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti1, category="NSTI=1")
blueberry_picrust2_ko_nsti0.5_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.5, category="NSTI=0.5")
blueberry_picrust2_ko_nsti0.25_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.25, category="NSTI=0.25")
blueberry_picrust2_ko_nsti0.1_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.1, category="NSTI=0.1")
blueberry_picrust2_ko_nsti0.05_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust2_ko_nsti0.05, category="NSTI=0.05")
blueberry_picrust1_ko_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$picrust1_ko, category="PICRUSt1")
blueberry_piphillin_ko_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$piphillin, category="Piphillin")
blueberry_panfp_ko_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$panfp, category="PanFP")
blueberry_tax4fun_ko_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_infiles$all_kos_overlap$tax4fun, category="Tax4Fun")
blueberry_null_ko_metrics <- calc_accuracy_metrics(blueberry_infiles$all_kos_overlap$mgs_ko, blueberry_ko_mgs_null_df_round, category="Null")

hmp_acc <- rbind(hmp_null_ko_metrics,
                 hmp_tax4fun_ko_metrics,
                 hmp_panfp_ko_metrics,
                 hmp_piphillin_ko_metrics,
                 hmp_picrust1_ko_metrics,
                 hmp_picrust2_ko_nsti2_gg_metrics,
                 hmp_picrust2_ko_nsti2_metrics,
                 hmp_picrust2_ko_nsti1.5_metrics, 
                 hmp_picrust2_ko_nsti1_metrics,
                 hmp_picrust2_ko_nsti0.5_metrics,
                 hmp_picrust2_ko_nsti0.25_metrics,
                 hmp_picrust2_ko_nsti0.1_metrics,
                 hmp_picrust2_ko_nsti0.05_metrics)

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

ocean_acc <- rbind(ocean_null_ko_metrics,
                 ocean_tax4fun_ko_metrics,
                 ocean_panfp_ko_metrics,
                 ocean_piphillin_ko_metrics,
                 ocean_picrust1_ko_metrics,
                 ocean_picrust2_ko_nsti2_gg_metrics,
                 ocean_picrust2_ko_nsti2_metrics,
                 ocean_picrust2_ko_nsti1.5_metrics, 
                 ocean_picrust2_ko_nsti1_metrics,
                 ocean_picrust2_ko_nsti0.5_metrics,
                 ocean_picrust2_ko_nsti0.25_metrics,
                 ocean_picrust2_ko_nsti0.1_metrics,
                 ocean_picrust2_ko_nsti0.05_metrics)

blueberry_acc <- rbind(blueberry_null_ko_metrics,
                 blueberry_tax4fun_ko_metrics,
                 blueberry_panfp_ko_metrics,
                 blueberry_piphillin_ko_metrics,
                 blueberry_picrust1_ko_metrics,
                 blueberry_picrust2_ko_nsti2_gg_metrics,
                 blueberry_picrust2_ko_nsti2_metrics,
                 blueberry_picrust2_ko_nsti1.5_metrics, 
                 blueberry_picrust2_ko_nsti1_metrics,
                 blueberry_picrust2_ko_nsti0.5_metrics,
                 blueberry_picrust2_ko_nsti0.25_metrics,
                 blueberry_picrust2_ko_nsti0.1_metrics,
                 blueberry_picrust2_ko_nsti0.05_metrics)

hmp_acc$category <- factor(hmp_acc$category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

mammal_acc$category <- factor(mammal_acc$category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

ocean_acc$category <- factor(ocean_acc$category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

blueberry_acc$category <- factor(blueberry_acc$category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2",
                                                      "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

saveRDS(object = hmp_ko_spearman_df, file = "../../data/saved_RDS/16S_vs_MGS_metrics/hmp_ko_spearman_df.rds")
saveRDS(object = hmp_acc, file = "../../data/saved_RDS/16S_vs_MGS_metrics/hmp_ko_acc_metrics.rds")

saveRDS(object = mammal_ko_spearman_df, file = "../../data/saved_RDS/16S_vs_MGS_metrics/mammal_ko_spearman_df.rds")
saveRDS(object = mammal_acc, file = "../../data/saved_RDS/16S_vs_MGS_metrics/mammal_ko_acc_metrics.rds")

saveRDS(object = ocean_ko_spearman_df, file = "../../data/saved_RDS/16S_vs_MGS_metrics/ocean_ko_spearman_df.rds")
saveRDS(object = ocean_acc, file = "../../data/saved_RDS/16S_vs_MGS_metrics/ocean_ko_acc_metrics.rds")

saveRDS(object = blueberry_ko_spearman_df, file = "../../data/saved_RDS/16S_vs_MGS_metrics/blueberry_ko_spearman_df.rds")
saveRDS(object = blueberry_acc, file = "../../data/saved_RDS/16S_vs_MGS_metrics/blueberry_ko_acc_metrics.rds")
