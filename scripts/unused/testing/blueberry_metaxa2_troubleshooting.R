blueberry_metaxa2_summary <- as.data.frame(t(read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/blueberry_metaxa2_kingdom_summary.txt",
                                        header=T, sep="\t", row.names=1)))

blueberry_metaxa2_summary$Prokaryota <- blueberry_metaxa2_summary$Archaea + blueberry_metaxa2_summary$Bacteria
blueberry_metaxa2_summary <- blueberry_metaxa2_summary[, -c(which(colnames(blueberry_metaxa2_summary) %in% c("Bacteria", "Archaea")))]
blueberry_metaxa2_summary$Total <- rowSums(blueberry_metaxa2_summary)

blueberry_metaxa2_summary$Prokaryota_per <- (blueberry_metaxa2_summary$Prokaryota/blueberry_metaxa2_summary$Total)*100
blueberry_metaxa2_summary$Eukaryota_per <- (blueberry_metaxa2_summary$Eukaryota/blueberry_metaxa2_summary$Total)*100

blueberry_metaxa2_summary$Pro_to_euk <- log2(blueberry_metaxa2_summary$Prokaryota/blueberry_metaxa2_summary$Eukaryota)

blueberry_ko_spearman <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_ko_spearman_df.rds")
blueberry_ko_spearman <- blueberry_ko_spearman[which(blueberry_ko_spearman$cat=="NSTI=2"),]
rownames(blueberry_ko_spearman) <- as.character(blueberry_ko_spearman$sample_names)
rownames(blueberry_ko_spearman) <- gsub("Bact", "BB", rownames(blueberry_ko_spearman))
rownames(blueberry_ko_spearman) <- gsub("\\.", "_", rownames(blueberry_ko_spearman))

blueberry_ko_acc <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_acc_metrics.rds")
blueberry_ko_acc <- blueberry_ko_acc[which(blueberry_ko_acc == "NSTI=2"),]
rownames(blueberry_ko_acc) <- as.character(blueberry_ko_acc$sample)
rownames(blueberry_ko_acc) <- gsub("Bact", "BB", rownames(blueberry_ko_acc))
rownames(blueberry_ko_acc) <- gsub("\\.", "_", rownames(blueberry_ko_acc))

blue_18S_weighted_nsti <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline/ec_18S_metagenome_out/weighted_nsti.tsv",
                                     header=T, sep="\t", row.names=1)
rownames(blue_18S_weighted_nsti) <- gsub("^", "BB", rownames(blue_18S_weighted_nsti))
rownames(blue_18S_weighted_nsti) <- gsub("-", "_", rownames(blue_18S_weighted_nsti))

par(mfrow=c(2,2))
plot(blueberry_metaxa2_summary$Prokaryota_per,
     blueberry_ko_spearman[rownames(blueberry_metaxa2_summary), "metric"],
     ylab="Spearman", xlab="% Prokaryotic", pch=16)
plot(blueberry_metaxa2_summary$Prokaryota_per,
     blueberry_ko_acc[rownames(blueberry_metaxa2_summary), "precision"],
     ylab="Precision", xlab="% Prokaryotic", pch=16)
plot(blueberry_metaxa2_summary$Prokaryota_per,
     blueberry_ko_acc[rownames(blueberry_metaxa2_summary), "recall"],
     ylab="Recall", xlab="% Prokaryotic", pch=16)
plot(blueberry_metaxa2_summary$Prokaryota_per,
     blue_18S_weighted_nsti[rownames(blueberry_metaxa2_summary), "weighted_NSTI"],
     ylab="Weighted NSTI", xlab="% Prokaryotic", pch=16)