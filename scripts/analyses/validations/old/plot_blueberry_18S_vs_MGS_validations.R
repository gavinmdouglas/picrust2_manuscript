library(ggplot2)
library(reshape2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

blueberry_18S_ec_acc <- readRDS("blueberry_18S_ec_acc_metrics.rds")
blueberry_18S_ec_scc <- readRDS("blueberry_18S_ec_spearman_df.rds")
blueberry_18S_pathabun_acc <- readRDS("blueberry_18S_pathabun_acc_metrics.rds")
blueberry_18S_pathabun_scc <- readRDS("blueberry_18S_pathabun_spearman_df.rds")
blueberry_18S_pathcov_acc <- readRDS("blueberry_18S_pathcov_acc_metrics.rds")
blueberry_18S_pathcov_scc <- readRDS("blueberry_18S_pathcov_spearman_df.rds")

combined_16S_18S_scc_contrib <- readRDS("blueberry_ec_16S_18S_contributions.rds")

blueberry_18S_ec_scc_subset <- blueberry_18S_ec_scc[,c("sample_names", "metric", "cat")]
colnames(blueberry_18S_ec_scc_subset) <- c("Sample", "Spearman's Correlation", "Category")
blueberry_18S_ec_scc_subset <- blueberry_18S_ec_scc_subset[
  with(blueberry_18S_ec_scc_subset, order(Category, Sample)),]
blueberry_18S_ec_scc_subset_melt <- melt(blueberry_18S_ec_scc_subset)

blueberry_18S_ec_acc_subset <- blueberry_18S_ec_acc[,c("sample", "precision", "recall", "category")]
colnames(blueberry_18S_ec_acc_subset) <- c("Sample", "Precision", "Recall", "Category")
blueberry_18S_ec_acc_subset <- blueberry_18S_ec_acc_subset[
  with(blueberry_18S_ec_acc_subset, order(Category, Sample)),]
blueberry_18S_ec_acc_subset_melt <- melt(blueberry_18S_ec_acc_subset)

blueberry_18S_ec_scc_acc_subset_melt <- rbind(blueberry_18S_ec_scc_subset_melt, blueberry_18S_ec_acc_subset_melt)

blueberry_18S_ec_scc_acc_subset_melt$Category <- factor(blueberry_18S_ec_scc_acc_subset_melt$Category,
                                                        levels=c("Null", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

blueberry_18S_ec_scc_acc_boxplots <- ggplot(blueberry_18S_ec_scc_acc_subset_melt, aes(x=Category, y=value, fill="#F8766D")) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + guides(fill=FALSE)


combined_16S_18S_scc_contrib_mean <- aggregate(.~contrib + sample, data=combined_16S_18S_scc_contrib, mean)

combined_16S_18S_scc_contrib_mean_max <- data.frame(matrix(NA, nrow=length(unique(combined_16S_18S_scc_contrib_mean$sample)), ncol=4))
rownames(combined_16S_18S_scc_contrib_mean_max) <- unique(combined_16S_18S_scc_contrib_mean$sample)
colnames(combined_16S_18S_scc_contrib_mean_max) <- c("rho_contrib", "precision_contrib", "recall_contrib", "F1_contrib")
for(samp in unique(combined_16S_18S_scc_contrib_mean$sample)) {
  combined_16S_18S_scc_contrib_mean_subset <- combined_16S_18S_scc_contrib_mean[combined_16S_18S_scc_contrib_mean$sample==samp,]
  max_rho_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$rho == max(combined_16S_18S_scc_contrib_mean_subset$rho)), "contrib"][1]
  max_precision_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$precision == max(combined_16S_18S_scc_contrib_mean_subset$precision)), "contrib"][1]
  max_recall_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$recall == max(combined_16S_18S_scc_contrib_mean_subset$recall)), "contrib"][1]
  max_F1_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$F1 == max(combined_16S_18S_scc_contrib_mean_subset$F1)), "contrib"][1]
  combined_16S_18S_scc_contrib_mean_max[samp, ] <- c(max_rho_contrib, max_precision_contrib, max_recall_contrib, max_F1_contrib)
}

combined_16S_18S_scc_contrib_mean_max_percent <- combined_16S_18S_scc_contrib_mean_max / 100

kingdom_count <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/sortmerna_out/kingdom_counts.tsv",
                            header=F, sep="\t", row.names=1)
colnames(kingdom_count) <- c("archaea", "bacteria", "eukaryota", "mixed")
rownames(kingdom_count) <- gsub("_16S_sortmerna_out.fastq.blast", "", rownames(kingdom_count))
kingdom_percent <- sweep(kingdom_count, 1, rowSums(kingdom_count), '/') * 100
kingdom_percent_tmp <- kingdom_percent
rownames(kingdom_percent_tmp) <- gsub("BB", "Blue", rownames(kingdom_percent_tmp) )
kingdom_percent_tmp <- kingdom_percent_tmp[rownames(combined_16S_18S_scc_contrib_mean_max),]

par(mfrow=c(2,2))
plot(combined_16S_18S_scc_contrib_mean_max_percent$rho_contrib, kingdom_percent_tmp$eukaryota, ylab="% Eukaryota", xlab="% 16S predictions", xlim=c(0,100), main="SCC")
plot(combined_16S_18S_scc_contrib_mean_max_percent$recall_contrib, kingdom_percent_tmp$eukaryota, ylab="% Eukaryota", xlab="% 16S predictions", xlim=c(0,100), main="Recall")
plot(combined_16S_18S_scc_contrib_mean_max_percent$precision_contrib, kingdom_percent_tmp$eukaryota, ylab="% Eukaryota", xlab="% 16S predictions", xlim=c(0,100), main="Precision")
#plot(combined_16S_18S_scc_contrib_mean_max_percent$F1_contrib  , kingdom_percent_tmp$eukaryota)
