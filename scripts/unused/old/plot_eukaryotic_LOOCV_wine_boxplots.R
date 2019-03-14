library(ggplot2)
library(cowplot)

setwd("/home/gavin/projects/picrust2_manuscript/data/taxa_LOOCV/")

source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

# Calculate null distribution for each E.C. database.
ec_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_18S.txt",
                     header=T, row.names=1, stringsAsFactors = FALSE, sep="\t", check.names=FALSE)

ec_ITS <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_ITS.txt",
                     header=T, row.names=1, stringsAsFactors = FALSE, sep="\t", check.names=FALSE)

ec_Class_18S <- read.table("18S_Class_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Order_18S <- read.table("18S_Order_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Family_18S <- read.table("18S_Family_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Genus_18S <- read.table("18S_Genus_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Species_18S <- read.table("18S_Species_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_assembly_18S <- read.table("18S_assembly_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)

ec_Class_ITS <- read.table("ITS_Class_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Order_ITS <- read.table("ITS_Order_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Family_ITS <- read.table("ITS_Family_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Genus_ITS <- read.table("ITS_Genus_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_Species_ITS <- read.table("ITS_Species_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_assembly_ITS <- read.table("ITS_assembly_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)

ec_assembly_18S_df <- data.frame(metric=ec_assembly_18S$spearman, cat="Genome",
                                 metric_type="spearman", sample_names=ec_assembly_18S$taxon)

ec_assembly_ITS_df <- data.frame(metric=ec_assembly_ITS$spearman, cat="Genome",
                                 metric_type="spearman", sample_names=ec_assembly_ITS$taxon)

ec_assembly_18S_expected <- t(read.table("18S_assembly_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_18S_spearman_assembly <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_assembly_18S_expected)
ec_assembly_18S_null <- cor_all_cols(tab1 = null_18S_spearman_assembly, tab2 = ec_assembly_18S_expected, cat_string="Genome", metric="spearman")

ec_assembly_ITS_expected <- t(read.table("ITS_assembly_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_ITS_spearman_assembly <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_assembly_ITS_expected)
ec_assembly_ITS_null <- cor_all_cols(tab1 = null_ITS_spearman_assembly, tab2 = ec_assembly_ITS_expected, cat_string="Genome", metric="spearman")


ec_Species_18S_df <- data.frame(metric=ec_Species_18S$spearman, cat="Species",
                                 metric_type="spearman", sample_names=ec_Species_18S$taxon)

ec_Species_ITS_df <- data.frame(metric=ec_Species_ITS$spearman, cat="Species",
                                metric_type="spearman", sample_names=ec_Species_ITS$taxon)

ec_Species_18S_expected <- t(read.table("18S_Species_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_18S_spearman_Species <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_Species_18S_expected)
ec_Species_18S_null <- cor_all_cols(tab1 = null_18S_spearman_Species, tab2 = ec_Species_18S_expected, cat_string="Species", metric="spearman")

ec_Species_ITS_expected <- t(read.table("ITS_Species_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_ITS_spearman_Species <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_Species_ITS_expected)
ec_Species_ITS_null <- cor_all_cols(tab1 = null_ITS_spearman_Species, tab2 = ec_Species_ITS_expected, cat_string="Species", metric="spearman")


ec_Genus_18S_df <- data.frame(metric=ec_Genus_18S$spearman, cat="Genus",
                                 metric_type="spearman", sample_names=ec_Genus_18S$taxon)

ec_Genus_ITS_df <- data.frame(metric=ec_Genus_ITS$spearman, cat="Genus",
                              metric_type="spearman", sample_names=ec_Genus_ITS$taxon)

ec_Genus_18S_expected <- t(read.table("18S_Genus_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_18S_spearman_Genus <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_Genus_18S_expected)
ec_Genus_18S_null <- cor_all_cols(tab1 = null_18S_spearman_Genus, tab2 = ec_Genus_18S_expected, cat_string="Genus", metric="spearman")

ec_Genus_ITS_expected <- t(read.table("ITS_Genus_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_ITS_spearman_Genus <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_Genus_ITS_expected)
ec_Genus_ITS_null <- cor_all_cols(tab1 = null_ITS_spearman_Genus, tab2 = ec_Genus_ITS_expected, cat_string="Genus", metric="spearman")



ec_Family_18S_df <- data.frame(metric=ec_Family_18S$spearman, cat="Family",
                                 metric_type="spearman", sample_names=ec_Family_18S$taxon)

ec_Family_ITS_df <- data.frame(metric=ec_Family_ITS$spearman, cat="Family",
                               metric_type="spearman", sample_names=ec_Family_ITS$taxon)

ec_Family_18S_expected <- t(read.table("18S_Family_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_18S_spearman_Family <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_Family_18S_expected)
ec_Family_18S_null <- cor_all_cols(tab1 = null_18S_spearman_Family, tab2 = ec_Family_18S_expected, cat_string="Family", metric="spearman")

ec_Family_ITS_expected <- t(read.table("ITS_Family_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_ITS_spearman_Family <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_Family_ITS_expected)
ec_Family_ITS_null <- cor_all_cols(tab1 = null_ITS_spearman_Family, tab2 = ec_Family_ITS_expected, cat_string="Family", metric="spearman")


ec_Order_18S_df <- data.frame(metric=ec_Order_18S$spearman, cat="Order",
                                 metric_type="spearman", sample_names=ec_Order_18S$taxon)

ec_Order_ITS_df <- data.frame(metric=ec_Order_ITS$spearman, cat="Order",
                              metric_type="spearman", sample_names=ec_Order_ITS$taxon)

ec_Order_18S_expected <- t(read.table("18S_Order_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_18S_spearman_Order <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_Order_18S_expected)
ec_Order_18S_null <- cor_all_cols(tab1 = null_18S_spearman_Order, tab2 = ec_Order_18S_expected, cat_string="Order", metric="spearman")

ec_Order_ITS_expected <- t(read.table("ITS_Order_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_ITS_spearman_Order <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_Order_ITS_expected)
ec_Order_ITS_null <- cor_all_cols(tab1 = null_ITS_spearman_Order, tab2 = ec_Order_ITS_expected, cat_string="Order", metric="spearman")


ec_Class_18S_df <- data.frame(metric=ec_Class_18S$spearman, cat="Class",
                                 metric_type="spearman", sample_names=ec_Class_18S$taxon)

ec_Class_ITS_df <- data.frame(metric=ec_Class_ITS$spearman, cat="Class",
                              metric_type="spearman", sample_names=ec_Class_ITS$taxon)

ec_Class_18S_expected <- t(read.table("18S_Class_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_18S_spearman_Class <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_Class_18S_expected)
ec_Class_18S_null <- cor_all_cols(tab1 = null_18S_spearman_Class, tab2 = ec_Class_18S_expected, cat_string="Class", metric="spearman")

ec_Class_ITS_expected <- t(read.table("ITS_Class_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE))
null_ITS_spearman_Class <- generate_null_mean_db_funcs(db = ec_ITS, tab = ec_Class_ITS_expected)
ec_Class_ITS_null <- cor_all_cols(tab1 = null_ITS_spearman_Class, tab2 = ec_Class_ITS_expected, cat_string="Class", metric="spearman")


ec_18S_null_df <- rbind(ec_Class_18S_null, ec_Order_18S_null, ec_Family_18S_null, ec_Genus_18S_null, ec_Species_18S_null, ec_assembly_18S_null)
ec_ITS_null_df <- rbind(ec_Class_ITS_null, ec_Order_ITS_null, ec_Family_ITS_null, ec_Genus_ITS_null, ec_Species_ITS_null, ec_assembly_ITS_null)
ec_18S_df <- rbind(ec_Class_18S_df, ec_Order_18S_df, ec_Family_18S_df, ec_Genus_18S_df, ec_Species_18S_df, ec_assembly_18S_df)
ec_ITS_df <- rbind(ec_Class_ITS_df, ec_Order_ITS_df, ec_Family_ITS_df, ec_Genus_ITS_df, ec_Species_ITS_df, ec_assembly_ITS_df)

ec_18S_null_df$Type <- "Null"
ec_ITS_null_df$Type <- "Null"
ec_18S_df$Type <- "Predicted"
ec_ITS_df$Type <- "Predicted" 

ec_18S_combined_df <- rbind(ec_18S_null_df, ec_18S_df)
ec_ITS_combined_df <- rbind(ec_ITS_null_df, ec_ITS_df)

ec_18S_combined_df$cat <- factor(ec_18S_combined_df$cat, levels=c("Genome", "Species", "Genus", "Family", "Order", "Class"))
ec_ITS_combined_df$cat <- factor(ec_ITS_combined_df$cat, levels=c("Genome", "Species", "Genus", "Family", "Order", "Class"))

ec_18S_combined_df_plot <- ggplot(ec_18S_combined_df, aes(x=Type, y=metric, fill=Type)) + geom_boxplot() +
  facet_grid(. ~ cat, scales = "free", space = "free", switch="x") +
  ylim(c(0.0, 1.0)) + ylab(c("")) + ggtitle("") +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  ylab("Spearman correlation coefficient") + xlab("") +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

ec_ITS_combined_df_plot <- ggplot(ec_ITS_combined_df, aes(x=Type, y=metric, fill=Type)) + geom_boxplot() +
  facet_grid(. ~ cat, scales = "free", space = "free", switch="x") +
  ylim(c(0.0, 1.0)) + ylab(c("")) + ggtitle("") +
  scale_fill_manual(values=c("light grey", "#00BFC4")) +
  ylab("Spearman correlation coefficient") + xlab("") +  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

plot_grid(ec_18S_combined_df_plot, ec_ITS_combined_df_plot, labels=c("A", "B"))


wilcox.test(ec_Class_18S_null$metric, ec_Class_18S_df$metric)
wilcox.test(ec_Order_18S_null$metric, ec_Order_18S_df$metric)
wilcox.test(ec_Family_18S_null$metric, ec_Family_18S_df$metric)
wilcox.test(ec_Genus_18S_null$metric, ec_Genus_18S_df$metric)
wilcox.test(ec_Species_18S_null$metric, ec_Species_18S_df$metric)
wilcox.test(ec_assembly_18S_null$metric, ec_assembly_18S_df$metric)

wilcox.test(ec_Class_ITS_null$metric, ec_Class_ITS_df$metric)
wilcox.test(ec_Order_ITS_null$metric, ec_Order_ITS_df$metric)
wilcox.test(ec_Family_ITS_null$metric, ec_Family_ITS_df$metric)
wilcox.test(ec_Genus_ITS_null$metric, ec_Genus_ITS_df$metric)
wilcox.test(ec_Species_ITS_null$metric, ec_Species_ITS_df$metric)
wilcox.test(ec_assembly_ITS_null$metric, ec_assembly_ITS_df$metric)

setwd("/home/gavin/projects/picrust2_manuscript/data/taxa_LOOCV/")
ec_ass_18S_expected <- read.table("18S_assembly_exp.tsv", header=T, sep="\t", row.names=1, check.names=FALSE)
ec_ass_18S_pred <- read.table("18S_assembly_predict.tsv", header=T, sep="\t", row.names=1, check.names=FALSE)
ec_ass_18S_pred <- ec_ass_18S_pred[,-which(colnames(ec_ass_18S_pred) == "nsti")]
ec_ass_18S_pred <- ec_ass_18S_pred[,-which(colnames(ec_ass_18S_pred) == "taxon")]

ec_ass_18S_metrics <- calc_accuracy_metrics(t(ec_ass_18S_expected), t(ec_ass_18S_pred), "test")

aggregate(. ~ func, data=strat_table_asv_interest, FUN=sum)




#### Run validation on wine dataset
humann2_ec <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt",
                         header=F, stringsAsFactors = FALSE)$V1
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

colnames(wine_mgs_ec) <- gsub("_kneaddata_Abundance.RPKs", "", colnames(wine_mgs_ec))
rownames(wine_mgs_ec) <- paste("EC:", rownames(wine_mgs_ec), sep="")
colnames(wine_mgs_ec) <- wine_mapfile[colnames(wine_mgs_ec), "X16S_name"]

wine_predicted_ec_nsti0.05_all <- add_missing_funcs(wine_predicted_ec_nsti0.05, colnames(ec_ITS))
wine_predicted_ec_nsti0.1_all <- add_missing_funcs(wine_predicted_ec_nsti0.1, colnames(ec_ITS))
wine_predicted_ec_nsti0.25_all <- add_missing_funcs(wine_predicted_ec_nsti0.25, colnames(ec_ITS))
wine_predicted_ec_nsti0.5_all <- add_missing_funcs(wine_predicted_ec_nsti0.5, colnames(ec_ITS))
wine_predicted_ec_nsti1_all <- add_missing_funcs(wine_predicted_ec_nsti1, colnames(ec_ITS))
wine_predicted_ec_nsti1.5_all <- add_missing_funcs(wine_predicted_ec_nsti1.5, colnames(ec_ITS))
wine_predicted_ec_nsti2_all <- add_missing_funcs(wine_predicted_ec_nsti2, colnames(ec_ITS))

mgs_rows_to_remove <- which(rownames(wine_mgs_ec) %in% c("EC:UNMAPPED", "EC:UNGROUPED"))
wine_mgs_ec <- wine_mgs_ec[-mgs_rows_to_remove, ]
wine_mgs_ec_all <- add_missing_funcs(wine_mgs_ec, humann2_ec)

wine_ITS_overlapping_ecs <- humann2_ec[which(humann2_ec %in% colnames(ec_ITS))]

wine_mgs_ec_all <- wine_mgs_ec_all[wine_ITS_overlapping_ecs,]

wine_ec_null_df <- generate_null_mean_db_funcs(db = ec_ITS, tab = wine_mgs_ec_all)
wine_ec_null_df_round <- round(wine_ec_null_df  - 0.00000001)

wine_ec_null_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_ec_null_df, metric = "spearman", cat_string = "Null")
wine_ec_nsti2_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti2_all, metric = "spearman", cat_string = "NSTI=2")
wine_ec_nsti1.5_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti1.5_all, metric = "spearman", cat_string = "NSTI=1.5")
wine_ec_nsti1_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti1_all, metric = "spearman", cat_string = "NSTI=1")
wine_ec_nsti0.5_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti0.5_all, metric = "spearman", cat_string = "NSTI=0.5")
wine_ec_nsti0.25_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti0.25_all, metric = "spearman", cat_string = "NSTI=0.25")
wine_ec_nsti0.1_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti0.1_all, metric = "spearman", cat_string = "NSTI=0.1")
wine_ec_nsti0.05_cor <- cor_all_cols(tab1 = wine_mgs_ec_all, tab2 = wine_predicted_ec_nsti0.05_all, metric = "spearman", cat_string = "NSTI=0.05")

wine_ec_null_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_ec_null_df_round, category="Null")
wine_ec_nsti2_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti2_all, category="NSTI=2")
wine_ec_nsti1.5_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti1.5_all, category="NSTI=1.5")
wine_ec_nsti1_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti1_all, category="NSTI=1")
wine_ec_nsti0.5_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti0.5_all, category="NSTI=0.5")
wine_ec_nsti0.25_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti0.25_all, category="NSTI=0.25")
wine_ec_nsti0.1_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti0.1_all, category="NSTI=0.1")
wine_ec_nsti0.05_metrics <- calc_accuracy_metrics(wine_mgs_ec_all, wine_predicted_ec_nsti0.05_all, category="NSTI=0.05")

wine_ec_euk_cor <- rbind(wine_ec_null_cor,
                             wine_ec_nsti2_cor,
                             wine_ec_nsti1.5_cor,
                             wine_ec_nsti1_cor,
                             wine_ec_nsti0.5_cor,
                             wine_ec_nsti0.25_cor,
                             wine_ec_nsti0.1_cor,
                             wine_ec_nsti0.05_cor)

wine_ec_euk_metrics <- rbind(wine_ec_null_metrics,
                        wine_ec_nsti2_metrics,
                        wine_ec_nsti1.5_metrics,
                        wine_ec_nsti1_metrics,
                        wine_ec_nsti0.5_metrics,
                        wine_ec_nsti0.25_metrics,
                        wine_ec_nsti0.1_metrics,
                        wine_ec_nsti0.05_metrics)

wine_ec_euk_metrics_subset <- wine_ec_euk_metrics[,c("sample", "precision", "recall", "fpr", "category")]

colnames(wine_ec_euk_metrics_subset) <- c("Sample", "Precision", "Recall", "FPR", "Category")

wine_ec_euk_metrics_subset$Category <- factor(wine_ec_euk_metrics_subset$Category, levels=c("Null", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

wine_ec_euk_metrics_subset_melt <- melt(wine_ec_euk_metrics_subset)
wine_ec_euk_metrics_boxplot <- ggplot(wine_ec_euk_metrics_subset_melt, aes(x=Category, y=value)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.0, 1.0)) + ylab(c("")) + ggtitle("") + geom_boxplot(fill="#00BFC4")

wine_ec_euk_cor$cat <- factor(wine_ec_euk_cor$cat, levels=c("Null","NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
wine_ec_euk_cor_melt <- melt(wine_ec_euk_cor)
wine_ec_euk_cor_boxplot <- ggplot(wine_ec_euk_cor_melt, aes(x=cat, y=value)) + geom_boxplot() +
  coord_flip() + xlab(c("Spearman correlation coefficient")) +
  ylim(c(0.0, 1.0)) + ylab(c("")) + ggtitle("") + geom_boxplot(fill="#00BFC4")

plot_grid(ec_18S_combined_df_plot, ec_ITS_combined_df_plot,
          #wine_ec_euk_cor_boxplot, wine_ec_euk_metrics_boxplot,
          labels=c("A", "B"))
