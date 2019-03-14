source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

ec_euk <- read.table(gzfile("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_ITS.txt"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

humann2_ec <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt",
                         header=F, stringsAsFactors = FALSE)$V1

euk_overlapping_ec <- colnames(ec_euk)[which(colnames(ec_euk) %in% humann2_ec)]

wine_mapfile <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/wine_id_mapping.txt",
                           header=T, sep="\t", stringsAsFactors = FALSE, quote="", comment.char="")
rownames(wine_mapfile) <- wine_mapfile$MGS_run

wine_predicted_ec_nsti2 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                                header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_predicted_ec_nsti1.5 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_predicted_ec_nsti1 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_1/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_predicted_ec_nsti0.5 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_predicted_ec_nsti0.25 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_predicted_ec_nsti0.1 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_0.1/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_predicted_ec_nsti0.05 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_ITS_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_mgs_ec <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/mgs/humann2_final_out/humann2_level4ec_unstratified.tsv",
                                header=T, sep="\t", stringsAsFactors = FALSE, comment.char = "", row.names=1)

cor.test(wine_ec_null_df[rownames(wine_mgs_ec), "control_ferment1a"],
         wine_mgs_ec[rownames(wine_mgs_ec), "control_ferment1a"],
         method = "spearman")

cor.test(wine_predicted_ec[rownames(wine_mgs_ec), "control_ferment1a"],
         wine_mgs_ec[rownames(wine_mgs_ec), "control_ferment1a"],
         method = "spearman")


colnames(wine_mgs_ec) <- gsub("_kneaddata_Abundance.RPKs", "", colnames(wine_mgs_ec))
rownames(wine_mgs_ec) <- paste("EC:", rownames(wine_mgs_ec), sep="")
colnames(wine_mgs_ec) <- wine_mapfile[colnames(wine_mgs_ec), "X16S_name"]

wine_ec_null_mean_df <- generate_null_mean_db_funcs(db = ec_euk,
                                                    tab = wine_mgs_ec_all)

wine_ec_null_mean_df_rounded <- round(wine_ec_null_mean_df - 0.00000001)

ec_null_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_ec_null_mean_df, metric = "spearman")
ec_nsti2_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti2, metric = "spearman")
ec_nsti1.5_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti1.5, metric = "spearman")
ec_nsti1_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti1, metric = "spearman")
ec_nsti0.5_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.5, metric = "spearman")
ec_nsti0.25_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.25, metric = "spearman")
ec_nsti0.1_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.1, metric = "spearman")
ec_nsti0.05_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.05, metric = "spearman")

ec_nsti0.05_all <- add_missing_funcs(wine_predicted_ec_nsti0.05, colnames(ec_euk))
ec_nsti0.1_all <- add_missing_funcs(wine_predicted_ec_nsti0.1, colnames(ec_euk))
ec_nsti0.25_all <- add_missing_funcs(wine_predicted_ec_nsti0.25, colnames(ec_euk))
ec_nsti0.5_all <- add_missing_funcs(wine_predicted_ec_nsti0.5, colnames(ec_euk))
ec_nsti1_all <- add_missing_funcs(wine_predicted_ec_nsti1, colnames(ec_euk))
ec_nsti1.5_all <- add_missing_funcs(wine_predicted_ec_nsti1.5, colnames(ec_euk))
ec_nsti2_all <- add_missing_funcs(wine_predicted_ec_nsti2, colnames(ec_euk))

mgs_rows_to_remove <- which(rownames(wine_mgs_ec) %in% c("EC:UNMAPPED", "EC:UNGROUPED"))
wine_mgs_ec <- wine_mgs_ec[-mgs_rows_to_remove, ]
wine_mgs_ec_all <- add_missing_funcs(wine_mgs_ec, humann2_ec)


ec_null_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_ec_null_mean_df_rounded, category="Null")
ec_nsti2_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti2_all, category="NSTI=2")
ec_nsti1.5_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti1.5_all, category="NSTI=1.5")
ec_nsti1_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti1_all, category="NSTI=1")
ec_nsti0.5_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti0.5_all, category="NSTI=0.5")
ec_nsti0.25_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti0.25_all, category="NSTI=0.25")
ec_nsti0.1_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti0.1_all, category="NSTI=0.1")
ec_nsti0.05_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, ec_nsti0.05_all, category="NSTI=0.05")

ec_euk_metrics <- rbind(ec_null_metrics,
                        ec_nsti2_metrics,
                        ec_nsti1.5_metrics,
                        ec_nsti1_metrics,
                        ec_nsti0.5_metrics,
                        ec_nsti0.25_metrics,
                        ec_nsti0.1_metrics,
                        ec_nsti0.05_metrics)

ec_euk_metrics_subset <- ec_euk_metrics[,c("sample", "precision", "recall", "fpr", "category")]

colnames(ec_euk_metrics_subset) <- c("Sample", "Precision", "Recall", "FPR", "Category")

ec_euk_metrics_subset$Category <- factor(ec_euk_metrics_subset$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

ec_euk_metrics_subset_melt <- melt(ec_euk_metrics_subset)
ggplot(ec_euk_metrics_subset_melt, aes(x=Category, y=value)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.0, 1.0)) + ylab(c("")) + ggtitle("Wine") + geom_boxplot(fill="#00BFC4")



ec_nsti2_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti2, metric = "spearman")
ec_nsti1.5_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti1.5, metric = "spearman")
ec_nsti1_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti1, metric = "spearman")
ec_nsti0.5_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.5, metric = "spearman")
ec_nsti0.25_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.25, metric = "spearman")
ec_nsti0.1_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.1, metric = "spearman")
ec_nsti0.05_cor <- cor_all_cols(tab1 = wine_mgs_ec, tab2 = wine_predicted_ec_nsti0.05, metric = "spearman")




write.table(x = colnames(wine_predicted_ec), file = "/home/gavin/tmp/wine_columns.txt",
            sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
