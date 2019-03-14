### Analyses to identify EC numbers specifically enriched (or present only in eukaryotes)
### in Eukaryotic databases compared to the prokaryotic database.

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

# First read in all 3 databases.
prokaryotic_ec <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ec.txt.gz"),
                             row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

eukaryotic_ec_uniref50 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_18S_uniref50.txt",
                                            row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

eukaryotic_ec_uniref50_newer_humann2 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_18S_uniref50_newer_humann2_test.txt",
                                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)


eukaryotic_ec_uniref90 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_18S.txt",
                                     row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt")$V1

# First identify EC numbers that are only found in eukaryotic databases.
eukaryotic_ec_uniref50_no_pro <- colnames(eukaryotic_ec_uniref50)[which(! colnames(eukaryotic_ec_uniref50) %in% colnames(prokaryotic_ec))]
eukaryotic_ec_uniref50_newer_humann2_no_pro <- colnames(eukaryotic_ec_uniref50_newer_humann2)[which(! colnames(eukaryotic_ec_uniref50_newer_humann2) %in% colnames(prokaryotic_ec))]
eukaryotic_ec_uniref90_no_pro <- colnames(eukaryotic_ec_uniref90)[which(! colnames(eukaryotic_ec_uniref90) %in% colnames(prokaryotic_ec))]

# Compare how common the unique ECs to eukaryotes are compared to all others.
eukaryotic_ec_uniref50_colSums <- colSums(eukaryotic_ec_uniref50 > 0)
boxplot(eukaryotic_ec_uniref50_colSums[eukaryotic_ec_uniref50_no_pro],
        eukaryotic_ec_uniref50_colSums[-which(names(eukaryotic_ec_uniref50_colSums) %in% eukaryotic_ec_uniref50_no_pro)])

eukaryotic_ec_uniref90_colSums <- colSums(eukaryotic_ec_uniref90 > 0)
boxplot(eukaryotic_ec_uniref90_colSums[eukaryotic_ec_uniref90_no_pro],
        eukaryotic_ec_uniref90_colSums[-which(names(eukaryotic_ec_uniref90_colSums) %in% eukaryotic_ec_uniref90_no_pro)])

# Test specifically prediction performance on these enriched ECs on blueberry dataset:
ec_blue_18S_nsti2_uniref90 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_18S_nsti2_uniref50 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline_uniref50_unstrat/EC_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                                         header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_18S_nsti0.05_uniref90 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv",
                                         header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_18S_nsti0.01_uniref90 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_0.01/pred_metagenome_unstrat.tsv",
                                            header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_18S_nsti0.05_uniref50 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline_uniref50_unstrat/EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv",
                                         header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)


ec_blue_mgs <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/humann2_output_final/humann2_level4ec_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")

colnames(ec_blue_18S_nsti2_uniref50) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti2_uniref50))
colnames(ec_blue_18S_nsti2_uniref90) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti2_uniref90))

rownames(ec_blue_mgs) <- gsub("^", "EC:", rownames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("_Abundance.RPKs", "", colnames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("BB", "Blue", colnames(ec_blue_mgs))

ec_blue_mgs <- ec_blue_mgs[-which(rownames(ec_blue_mgs) %in% c("EC:UNMAPPED", "EC:UNGROUPED")),]
ec_blue_mgs_nomiss <- add_missing_funcs(ec_blue_mgs, possible_mgs_ecs)

ec_blue_18S_nsti2_uniref50_nomiss <- add_missing_funcs(ec_blue_18S_nsti2_uniref50, colnames(eukaryotic_ec_uniref50))
ec_blue_18S_nsti2_uniref90_nomiss <- add_missing_funcs(ec_blue_18S_nsti2_uniref90, colnames(eukaryotic_ec_uniref90))

ec_blue_18S_nsti2_uniref50_nomiss <- add_missing_funcs(ec_blue_18S_nsti2_uniref50, colnames(eukaryotic_ec_uniref50))
ec_blue_18S_nsti2_uniref90_nomiss <- add_missing_funcs(ec_blue_18S_nsti2_uniref90, colnames(eukaryotic_ec_uniref90))

ec_blue_18S_nsti2_uniref50_newer_humann2 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline_uniref50_unstrat_new_humann2_test/ec_18S_uniref50_newer_humann2_test_metagenome_out/pred_metagenome_unstrat.tsv",
                                         header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

colnames(ec_blue_18S_nsti2_uniref50_newer_humann2) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti2_uniref50_newer_humann2))
ec_blue_18S_nsti2_uniref50_newer_humann2_nomiss <- add_missing_funcs(ec_blue_18S_nsti2_uniref50_newer_humann2, colnames(eukaryotic_ec_uniref50_newer_humann2))
SCC_ec_blue_18S_nsti2_uniref50_newer_humann2_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti2_uniref50_newer_humann2_nomiss,
                                                         tab2 = ec_blue_mgs_nomiss,
                                                         cat_string="uniref50",
                                                         metric="spearman")

ec_blue_18S_nsti0.05_uniref50_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.05_uniref50, colnames(eukaryotic_ec_uniref50))
ec_blue_18S_nsti0.05_uniref90_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.05_uniref90, colnames(eukaryotic_ec_uniref90))

SCC_ec_blue_18S_nsti0.05_uniref50_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti0.05_uniref50_nomiss,
                                                         tab2 = ec_blue_mgs_nomiss,
                                                         cat_string="uniref50",
                                                         metric="spearman")

colnames(ec_blue_18S_nsti0.05_uniref50) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti0.05_uniref50))
colnames(ec_blue_18S_nsti0.05_uniref90) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti0.05_uniref90))
colnames(ec_blue_18S_nsti0.01_uniref90) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti0.01_uniref90))

ec_blue_18S_nsti0.05_uniref50_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.05_uniref50, colnames(eukaryotic_ec_uniref50))
ec_blue_18S_nsti0.05_uniref90_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.05_uniref90, colnames(eukaryotic_ec_uniref90))
ec_blue_18S_nsti0.01_uniref90_nomiss <- add_missing_funcs(ec_blue_18S_nsti0.01_uniref90, colnames(eukaryotic_ec_uniref90))



SCC_ec_blue_18S_nsti0.05_uniref50_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti0.05_uniref50_nomiss,
                                                         tab2 = ec_blue_mgs_nomiss,
                                                         cat_string="uniref50",
                                                         metric="spearman")

SCC_ec_blue_18S_nsti0.05_uniref90_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti0.05_uniref90_nomiss,
                                                      tab2 = ec_blue_mgs_nomiss,
                                                      cat_string="uniref90",
                                                      metric="spearman")

SCC_ec_blue_18S_nsti0.01_uniref90_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti0.01_uniref90_nomiss,
                                                         tab2 = ec_blue_mgs_nomiss,
                                                         cat_string="uniref90",
                                                         metric="spearman")

acc_ec_blue_18S_nsti2_uniref90_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                               ec_blue_18S_nsti2_uniref90_nomiss,
                                                               category="uniref50")

acc_ec_blue_18S_null_uniref90_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                               blue_ec_uniref90_mgs_null_df,
                                                               category="uniref50")

acc_ec_blue_18S_null_uniref50_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                              blue_ec_uniref50_mgs_null_df,
                                                              category="uniref50")

acc_ec_blue_18S_nsti2_uniref90_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                                  ec_blue_18S_nsti2_uniref90_nomiss,
                                                                  category="uniref50")

acc_ec_blue_18S_nsti2_uniref50_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                               ec_blue_18S_nsti2_uniref50_nomiss,
                                                               category="test")

acc_ec_blue_18S_nsti0.05_uniref50_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                                  ec_blue_18S_nsti0.05_uniref50_nomiss,
                                                                  category="uniref50")

acc_ec_blue_18S_nsti0.05_uniref90_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                                  ec_blue_18S_nsti0.05_uniref90_nomiss,
                                                                             category="uniref50")

acc_ec_blue_18S_nsti0.01_uniref90_nomiss <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                                  ec_blue_18S_nsti0.01_uniref90_nomiss,
                                                                  category="uniref50")

par(mfrow=c(2,2))
boxplot(acc_ec_blue_18S_null_uniref90_nomiss$precision, acc_ec_blue_18S_nsti2_uniref90_nomiss$precision, acc_ec_blue_18S_nsti0.05_uniref90_nomiss$precision, acc_ec_blue_18S_nsti0.01_uniref90_nomiss$precision)
boxplot(acc_ec_blue_18S_null_uniref90_nomiss$recall, acc_ec_blue_18S_nsti2_uniref90_nomiss$recall, acc_ec_blue_18S_nsti0.05_uniref90_nomiss$recall, acc_ec_blue_18S_nsti0.01_uniref90_nomiss$recall)
boxplot(acc_ec_blue_18S_null_uniref50_nomiss$precision, acc_ec_blue_18S_nsti2_uniref50_nomiss$precision, acc_ec_blue_18S_nsti0.05_uniref50_nomiss$precision)
boxplot(acc_ec_blue_18S_null_uniref50_nomiss$recall, acc_ec_blue_18S_nsti2_uniref50_nomiss$recall, acc_ec_blue_18S_nsti0.05_uniref50_nomiss$recall)

boxplot(acc_ec_blue_18S_nsti2_uniref50_nomiss$precision, acc_ec_blue_18S_nsti0.05_uniref50_nomiss$precision)

#ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset <- ec_blue_18S_nsti2_uniref50_nomiss[eukaryotic_ec_uniref50_no_pro,]

#ec_blue_18S_nsti2_uniref90_no_pro_subset <- ec_blue_18S_nsti2_uniref90[eukaryotic_ec_uniref90_no_pro,]
blue_ec_uniref90_mgs_null_df <- generate_null_mean_db_funcs(db = eukaryotic_ec_uniref90, tab = ec_blue_mgs_nomiss)

SCC_ec_blue_18S_null_uniref90 <- cor_all_cols(tab1 = blue_ec_uniref90_mgs_null_df,
                                                      tab2 = ec_blue_mgs_nomiss,
                                                      cat_string="Null",
                                                      metric="spearman")

blue_ec_uniref50_mgs_null_df <- generate_null_mean_db_funcs(db = eukaryotic_ec_uniref50, tab = ec_blue_mgs_nomiss)

SCC_ec_blue_18S_null_uniref50 <- cor_all_cols(tab1 = blue_ec_uniref50_mgs_null_df,
                                              tab2 = ec_blue_mgs_nomiss,
                                              cat_string="Null",
                                              metric="spearman")

boxplot(SCC_ec_blue_18S_null_uniref90$metric, SCC_ec_blue_18S_null_uniref50$metric, SCC_ec_blue_18S_nsti2_uniref90_nomiss$metric, SCC_ec_blue_18S_nsti2_uniref50_nomiss$metric)

SCC_ec_blue_18S_nsti2_uniref50_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti2_uniref50_nomiss,
                                                             tab2 = ec_blue_mgs_nomiss,
                                                             cat_string="uniref50",
                                                             metric="spearman")

SCC_ec_blue_18S_nsti2_uniref90_nomiss <- cor_all_cols(tab1 = ec_blue_18S_nsti2_uniref90_nomiss,
                                                      tab2 = ec_blue_mgs_nomiss,
                                                      cat_string="uniref90",
                                                      metric="spearman")

acc_ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset <- calc_accuracy_metrics(ec_blue_mgs_nomiss,
                                                                             ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset,
                                                                             category="uniref50")

SCC_ec_blue_18S_nsti2_uniref90_no_pro_subset <- cor_all_cols(tab1 = ec_blue_18S_nsti2_uniref90,
                                                             tab2 = ec_blue_mgs,
                                                             cat_string="uniref90",
                                                             metric="spearman")


ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset <- ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset[eukaryotic_ec_uniref50_no_pro,]
tmp2 <- rownames(ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset)[which(rownames(ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset) %in% rownames(ec_blue_mgs_nomiss))]
tmp1 <- ec_blue_mgs_nomiss[tmp2, "Blue197"]
tmp3 <- ec_blue_18S_nsti2_uniref50_nomiss_no_pro_subset[tmp2, "Blue197"]
plot(tmp1, tmp3)
eukaryotic_ec_uniref90_no_pro_subset_colSums <- colSums(eukaryotic_ec_uniref90[,tmp2])

par(mfrow=c(2,1))
plot(tmp1, eukaryotic_ec_uniref90_no_pro_subset_colSums)
plot(tmp3, eukaryotic_ec_uniref90_no_pro_subset_colSums)


