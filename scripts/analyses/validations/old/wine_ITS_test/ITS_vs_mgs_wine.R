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

wine_id_map <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/wine_id_mapping.txt", header=T, sep="\t", stringsAsFactors = FALSE)
rownames(wine_id_map) <- wine_id_map$MGS_run

# Read in MGS.
ec_wine_mgs <- read.table("../../../mgs/humann2_final_out/humann2_level4ec_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")
rownames(ec_wine_mgs) <- gsub("^", "EC:", rownames(ec_wine_mgs))
colnames(ec_wine_mgs) <- gsub("_kneaddata_Abundance.RPKs", "", colnames(ec_wine_mgs))
ec_wine_mgs <- ec_wine_mgs[-which(rownames(ec_wine_mgs) %in% c("EC:UNMAPPED", "EC:UNGROUPED")),]
colnames(ec_wine_mgs) <- wine_id_map[colnames(ec_wine_mgs), "X16S_name"]

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

# Subset to ECs and samples overlapping..
ec_wine_mgs_nomiss <- add_missing_funcs(ec_wine_mgs, possible_mgs_ecs)
ec_wine_mgs_nomiss_subset <- ec_wine_mgs_nomiss[overlapping_ec, overlapping_samples]

ec_wine_ITS_nsti2 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv",
                                                    possible_ITS_ecs, overlapping_ec, overlapping_samples)

ec_wine_ITS_nsti1.5 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv",
                                                      possible_ITS_ecs, overlapping_ec, overlapping_samples)

ec_wine_ITS_nsti1 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_1.0/pred_metagenome_unstrat.tsv",
                                                    possible_ITS_ecs, overlapping_ec, overlapping_samples)

ec_wine_ITS_nsti0.5 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv",
                                                      possible_ITS_ecs, overlapping_ec, overlapping_samples)

ec_wine_ITS_nsti0.25 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv",
                                                      possible_ITS_ecs, overlapping_ec, overlapping_samples)

ec_wine_ITS_nsti0.1 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.1/pred_metagenome_unstrat.tsv",
                                                      possible_ITS_ecs, overlapping_ec, overlapping_samples)

ec_wine_ITS_nsti0.05 <- read_in_wine_ITS_picrust2("EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv",
                                                       possible_ITS_ecs, overlapping_ec, overlapping_samples)


# Read in % of each kingdom identified in MGS.
wine_kingdom_counts <- read.table("mgs/wine_metaxa2_kingdom_summary.txt",
                                  header=T, sep="\t", stringsAsFactors = FALSE, row.names = 1)

# Get % eukaryota.
wine_kingdom_percent <- data.frame(t(sweep(wine_kingdom_counts, 2, colSums(wine_kingdom_counts), '/')))*100



wine_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_wine_mgs_nomiss_subset)
wine_ec_mgs_null_df_round <- round(wine_ec_mgs_null_df  - 0.00000001)

wine_ec_mgs_null_scc <- cor_all_cols(tab1 = wine_ec_mgs_null_df, tab2 = ec_wine_mgs_nomiss_subset, cat_string="Null", metric="spearman")
wine_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, wine_ec_mgs_null_df, category="Null")

wine_ec_mgs_nsti2_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti2, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=2", metric="spearman")
wine_ec_mgs_nsti2_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti2, category="NSTI=2")
rownames(wine_ec_mgs_nsti2_metrics) <- wine_ec_mgs_nsti2_metrics$sample
wine_ec_mgs_nsti2_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti2_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti2_scc) <- wine_ec_mgs_nsti2_scc$sample_names
wine_ec_mgs_nsti2_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti2_scc$sample_names, "Eukaryota"]

wine_ec_mgs_nsti1.5_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti1.5, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=1.5", metric="spearman")
wine_ec_mgs_nsti1.5_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti1.5, category="NSTI=1.5")
rownames(wine_ec_mgs_nsti1.5_metrics) <- wine_ec_mgs_nsti1.5_metrics$sample
wine_ec_mgs_nsti1.5_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti1.5_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti1.5_scc) <- wine_ec_mgs_nsti1.5_scc$sample_names
wine_ec_mgs_nsti1.5_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti1.5_scc$sample_names, "Eukaryota"]

wine_ec_mgs_nsti1_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti1, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=1", metric="spearman")
wine_ec_mgs_nsti1_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti1, category="NSTI=1")
rownames(wine_ec_mgs_nsti1_metrics) <- wine_ec_mgs_nsti1_metrics$sample
wine_ec_mgs_nsti1_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti1_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti1_scc) <- wine_ec_mgs_nsti1_scc$sample_names
wine_ec_mgs_nsti1_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti1_scc$sample_names, "Eukaryota"]

wine_ec_mgs_nsti0.5_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.5, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.5", metric="spearman")
wine_ec_mgs_nsti0.5_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.5, category="NSTI=0.5")
rownames(wine_ec_mgs_nsti0.5_metrics) <- wine_ec_mgs_nsti0.5_metrics$sample
wine_ec_mgs_nsti0.5_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.5_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti0.5_scc) <- wine_ec_mgs_nsti0.5_scc$sample_names
wine_ec_mgs_nsti0.5_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.5_scc$sample_names, "Eukaryota"]

wine_ec_mgs_nsti0.25_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.25, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.25", metric="spearman")
wine_ec_mgs_nsti0.25_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.25, category="NSTI=0.25")
rownames(wine_ec_mgs_nsti0.25_metrics) <- wine_ec_mgs_nsti0.25_metrics$sample
wine_ec_mgs_nsti0.25_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.25_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti0.25_scc) <- wine_ec_mgs_nsti0.25_scc$sample_names
wine_ec_mgs_nsti0.25_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.25_scc$sample_names, "Eukaryota"]

wine_ec_mgs_nsti0.1_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.1, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.1", metric="spearman")
wine_ec_mgs_nsti0.1_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.1, category="NSTI=0.1")
rownames(wine_ec_mgs_nsti0.1_metrics) <- wine_ec_mgs_nsti0.1_metrics$sample
wine_ec_mgs_nsti0.1_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.1_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti0.1_scc) <- wine_ec_mgs_nsti0.1_scc$sample_names
wine_ec_mgs_nsti0.1_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.1_scc$sample_names, "Eukaryota"]

wine_ec_mgs_nsti0.05_scc <- cor_all_cols(tab1 = ec_wine_ITS_nsti0.05, tab2 = ec_wine_mgs_nomiss_subset, cat_string="NSTI=0.05", metric="spearman")
wine_ec_mgs_nsti0.05_metrics <- calc_accuracy_metrics(ec_wine_mgs_nomiss_subset, ec_wine_ITS_nsti0.05, category="NSTI=0.05")
rownames(wine_ec_mgs_nsti0.05_metrics) <- wine_ec_mgs_nsti0.05_metrics$sample
wine_ec_mgs_nsti0.05_metrics$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.05_metrics$sample, "Eukaryota"]
rownames(wine_ec_mgs_nsti0.05_scc) <- wine_ec_mgs_nsti0.05_scc$sample_names
wine_ec_mgs_nsti0.05_scc$eukaryota <- wine_kingdom_percent[wine_ec_mgs_nsti0.05_scc$sample_names, "Eukaryota"]


par(mfrow=c(1,1))
plot(wine_ec_mgs_nsti2_scc$eukaryota, wine_ec_mgs_nsti2_scc$metric, xlab="% Eukaryotic", ylab="Spearman ITS vs MGS", main="")


par(mfrow=c(1,1))
boxplot(wine_ec_mgs_null_scc$metric,
        wine_ec_mgs_nsti2_scc$metric,
        wine_ec_mgs_nsti1.5_scc$metric,
        wine_ec_mgs_nsti1_scc$metric,
        wine_ec_mgs_nsti0.5_scc$metric,
        wine_ec_mgs_nsti0.25_scc$metric,
        wine_ec_mgs_nsti0.1_scc$metric,
        wine_ec_mgs_nsti0.05_scc$metric,
        xlab="", ylab="Spearman ITS vs MGS", main="",
        names=c("Null", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"),
        col=c("light grey", rep(x = "#00BFC4", 7)),
        ylim=c(0.2, 0.7))

boxplot(wine_ec_mgs_null_metrics$precision,
        wine_ec_mgs_nsti2_metrics$precision,
        wine_ec_mgs_nsti1.5_metrics$precision,
        wine_ec_mgs_nsti1_metrics$precision,
        wine_ec_mgs_nsti0.5_metrics$precision,
        wine_ec_mgs_nsti0.25_metrics$precision,
        wine_ec_mgs_nsti0.1_metrics$precision,
        wine_ec_mgs_nsti0.05_metrics$precision,
        xlab="", ylab="Precision ITS vs MGS", main="",
        names=c("Null", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

boxplot(wine_ec_mgs_null_metrics$recall,
        wine_ec_mgs_nsti2_metrics$recall,
        wine_ec_mgs_nsti1.5_metrics$recall,
        wine_ec_mgs_nsti1_metrics$recall,
        wine_ec_mgs_nsti0.5_metrics$recall,
        wine_ec_mgs_nsti0.25_metrics$recall,
        wine_ec_mgs_nsti0.1_metrics$recall,
        wine_ec_mgs_nsti0.05_metrics$recall,
        xlab="", ylab="Recall ITS vs MGS", main="",
        names=c("Null", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
