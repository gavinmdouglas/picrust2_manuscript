### Make figure contrasting accuracy metrics across HSP tools on each 16S validation dataset.
### Also test for statistical significance between these categories (and save wilcoxon output to RDS).

library(ggplot2)
library(reshape2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

hmp_ec_acc_metrics <- readRDS("hmp_ec_acc_metrics.rds")
hmp_pathabun_acc_metrics <- readRDS("hmp_pathabun_acc_metrics.rds")
hmp_pathcov_acc_metrics <- readRDS("hmp_pathcov_acc_metrics.rds")

mammal_ec_acc_metrics <- readRDS("mammal_ec_acc_metrics.rds")
mammal_pathabun_acc_metrics <- readRDS("mammal_pathabun_acc_metrics.rds")
mammal_pathcov_acc_metrics <- readRDS("mammal_pathcov_acc_metrics.rds")

ocean_ec_acc_metrics <- readRDS("ocean_ec_acc_metrics.rds")
ocean_pathabun_acc_metrics <- readRDS("ocean_pathabun_acc_metrics.rds")
ocean_pathcov_acc_metrics <- readRDS("ocean_pathcov_acc_metrics.rds")

blueberry_ec_acc_metrics <- readRDS("blueberry_ec_acc_metrics.rds")
blueberry_pathabun_acc_metrics <- readRDS("blueberry_pathabun_acc_metrics.rds")
blueberry_pathcov_acc_metrics <- readRDS("blueberry_pathcov_acc_metrics.rds")


hmp_ec_acc_metrics_subset <- hmp_ec_acc_metrics[,c("sample", "precision", "recall", "category")]
colnames(hmp_ec_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
hmp_ec_acc_metrics_subset <- hmp_ec_acc_metrics_subset[
  with(hmp_ec_acc_metrics_subset, order(Category, Sample)),]
hmp_ec_acc_metrics_subset$Database <- "Other"
hmp_ec_acc_metrics_subset$Database[grep("NSTI" , hmp_ec_acc_metrics_subset$Category)] <- "PICRUSt2"
hmp_ec_acc_metrics_subset$Database[which(hmp_ec_acc_metrics_subset$Category=="Null")] <- "Null"
hmp_ec_acc_metrics_subset$Category <- factor(hmp_ec_acc_metrics_subset$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_ec_acc_metrics_subset_melt <- melt(hmp_ec_acc_metrics_subset)
hmp_ec_acc_metrics_boxplots <- ggplot(hmp_ec_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
                                   facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
                                   ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("HMP")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

hmp_pathabun_acc_metrics_subset <- hmp_pathabun_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(hmp_pathabun_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
hmp_pathabun_acc_metrics_subset <- hmp_pathabun_acc_metrics_subset[
  with(hmp_pathabun_acc_metrics_subset, order(Category, Sample)),]
hmp_pathabun_acc_metrics_subset$Database <- "PICRUSt2"
hmp_pathabun_acc_metrics_subset$Database[which(hmp_pathabun_acc_metrics_subset$Category=="Null")] <- "Null"
hmp_pathabun_acc_metrics_subset$Category <- factor(hmp_pathabun_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_pathabun_acc_metrics_subset_melt <- melt(hmp_pathabun_acc_metrics_subset)
hmp_pathabun_acc_metrics_boxplots <- ggplot(hmp_pathabun_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("HMP")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))

hmp_pathcov_acc_metrics_subset <- hmp_pathcov_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(hmp_pathcov_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
hmp_pathcov_acc_metrics_subset <- hmp_pathcov_acc_metrics_subset[
  with(hmp_pathcov_acc_metrics_subset, order(Category, Sample)),]
hmp_pathcov_acc_metrics_subset$Database <- "PICRUSt2"
hmp_pathcov_acc_metrics_subset$Database[which(hmp_pathcov_acc_metrics_subset$Category=="Null")] <- "Null"
hmp_pathcov_acc_metrics_subset$Category <- factor(hmp_pathcov_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_pathcov_acc_metrics_subset_melt <- melt(hmp_pathcov_acc_metrics_subset)
hmp_pathcov_acc_metrics_boxplots <- ggplot(hmp_pathcov_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("HMP")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))




mammal_ec_acc_metrics_subset <- mammal_ec_acc_metrics[,c("sample", "precision", "recall", "category")]
colnames(mammal_ec_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
mammal_ec_acc_metrics_subset <- mammal_ec_acc_metrics_subset[
  with(mammal_ec_acc_metrics_subset, order(Category, Sample)),]
mammal_ec_acc_metrics_subset$Database <- "Other"
mammal_ec_acc_metrics_subset$Database[grep("NSTI" , mammal_ec_acc_metrics_subset$Category)] <- "PICRUSt2"
mammal_ec_acc_metrics_subset$Database[which(mammal_ec_acc_metrics_subset$Category=="Null")] <- "Null"
mammal_ec_acc_metrics_subset$Category <- factor(mammal_ec_acc_metrics_subset$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_ec_acc_metrics_subset_melt <- melt(mammal_ec_acc_metrics_subset)
mammal_ec_acc_metrics_boxplots <- ggplot(mammal_ec_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Mammal")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

mammal_pathabun_acc_metrics_subset <- mammal_pathabun_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(mammal_pathabun_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
mammal_pathabun_acc_metrics_subset <- mammal_pathabun_acc_metrics_subset[
  with(mammal_pathabun_acc_metrics_subset, order(Category, Sample)),]
mammal_pathabun_acc_metrics_subset$Database <- "PICRUSt2"
mammal_pathabun_acc_metrics_subset$Database[which(mammal_pathabun_acc_metrics_subset$Category=="Null")] <- "Null"
mammal_pathabun_acc_metrics_subset$Category <- factor(mammal_pathabun_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_pathabun_acc_metrics_subset_melt <- melt(mammal_pathabun_acc_metrics_subset)
mammal_pathabun_acc_metrics_boxplots <- ggplot(mammal_pathabun_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Mammal")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))

mammal_pathcov_acc_metrics_subset <- mammal_pathcov_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(mammal_pathcov_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
mammal_pathcov_acc_metrics_subset <- mammal_pathcov_acc_metrics_subset[
  with(mammal_pathcov_acc_metrics_subset, order(Category, Sample)),]
mammal_pathcov_acc_metrics_subset$Database <- "PICRUSt2"
mammal_pathcov_acc_metrics_subset$Database[which(mammal_pathcov_acc_metrics_subset$Category=="Null")] <- "Null"
mammal_pathcov_acc_metrics_subset$Category <- factor(mammal_pathcov_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_pathcov_acc_metrics_subset_melt <- melt(mammal_pathcov_acc_metrics_subset)
mammal_pathcov_acc_metrics_boxplots <- ggplot(mammal_pathcov_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Mammal")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))




ocean_ec_acc_metrics_subset <- ocean_ec_acc_metrics[,c("sample", "precision", "recall", "category")]
colnames(ocean_ec_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
ocean_ec_acc_metrics_subset <- ocean_ec_acc_metrics_subset[
  with(ocean_ec_acc_metrics_subset, order(Category, Sample)),]
ocean_ec_acc_metrics_subset$Database <- "Other"
ocean_ec_acc_metrics_subset$Database[grep("NSTI" , ocean_ec_acc_metrics_subset$Category)] <- "PICRUSt2"
ocean_ec_acc_metrics_subset$Database[which(ocean_ec_acc_metrics_subset$Category=="Null")] <- "Null"
ocean_ec_acc_metrics_subset$Category <- factor(ocean_ec_acc_metrics_subset$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_ec_acc_metrics_subset_melt <- melt(ocean_ec_acc_metrics_subset)
ocean_ec_acc_metrics_boxplots <- ggplot(ocean_ec_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Ocean")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

ocean_pathabun_acc_metrics_subset <- ocean_pathabun_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(ocean_pathabun_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
ocean_pathabun_acc_metrics_subset <- ocean_pathabun_acc_metrics_subset[
  with(ocean_pathabun_acc_metrics_subset, order(Category, Sample)),]
ocean_pathabun_acc_metrics_subset$Database <- "PICRUSt2"
ocean_pathabun_acc_metrics_subset$Database[which(ocean_pathabun_acc_metrics_subset$Category=="Null")] <- "Null"
ocean_pathabun_acc_metrics_subset$Category <- factor(ocean_pathabun_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_pathabun_acc_metrics_subset_melt <- melt(ocean_pathabun_acc_metrics_subset)
ocean_pathabun_acc_metrics_boxplots <- ggplot(ocean_pathabun_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Ocean")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))

ocean_pathcov_acc_metrics_subset <- ocean_pathcov_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(ocean_pathcov_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
ocean_pathcov_acc_metrics_subset <- ocean_pathcov_acc_metrics_subset[
  with(ocean_pathcov_acc_metrics_subset, order(Category, Sample)),]
ocean_pathcov_acc_metrics_subset$Database <- "PICRUSt2"
ocean_pathcov_acc_metrics_subset$Database[which(ocean_pathcov_acc_metrics_subset$Category=="Null")] <- "Null"
ocean_pathcov_acc_metrics_subset$Category <- factor(ocean_pathcov_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_pathcov_acc_metrics_subset_melt <- melt(ocean_pathcov_acc_metrics_subset)
ocean_pathcov_acc_metrics_boxplots <- ggplot(ocean_pathcov_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Ocean")  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))




blueberry_ec_acc_metrics_subset <- blueberry_ec_acc_metrics[,c("sample", "precision", "recall", "category")]
colnames(blueberry_ec_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
blueberry_ec_acc_metrics_subset <- blueberry_ec_acc_metrics_subset[
  with(blueberry_ec_acc_metrics_subset, order(Category, Sample)),]
blueberry_ec_acc_metrics_subset$Database <- "Other"
blueberry_ec_acc_metrics_subset$Database[grep("NSTI" , blueberry_ec_acc_metrics_subset$Category)] <- "PICRUSt2"
blueberry_ec_acc_metrics_subset$Database[which(blueberry_ec_acc_metrics_subset$Category=="Null")] <- "Null"
blueberry_ec_acc_metrics_subset$Category <- factor(blueberry_ec_acc_metrics_subset$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_ec_acc_metrics_subset_melt <- melt(blueberry_ec_acc_metrics_subset)
blueberry_ec_acc_metrics_boxplots <- ggplot(blueberry_ec_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Soil (Blueberry)") + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

blueberry_pathabun_acc_metrics_subset <- blueberry_pathabun_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(blueberry_pathabun_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
blueberry_pathabun_acc_metrics_subset <- blueberry_pathabun_acc_metrics_subset[
  with(blueberry_pathabun_acc_metrics_subset, order(Category, Sample)),]
blueberry_pathabun_acc_metrics_subset$Database <- "PICRUSt2"
blueberry_pathabun_acc_metrics_subset$Database[which(blueberry_pathabun_acc_metrics_subset$Category=="Null")] <- "Null"
blueberry_pathabun_acc_metrics_subset$Category <- factor(blueberry_pathabun_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_pathabun_acc_metrics_subset_melt <- melt(blueberry_pathabun_acc_metrics_subset)
blueberry_pathabun_acc_metrics_boxplots <- ggplot(blueberry_pathabun_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Soil (Blueberry)") + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))

blueberry_pathcov_acc_metrics_subset <- blueberry_pathcov_acc_metrics[,c("sample", "precision", "recall",  "category")]
colnames(blueberry_pathcov_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")
blueberry_pathcov_acc_metrics_subset <- blueberry_pathcov_acc_metrics_subset[
  with(blueberry_pathcov_acc_metrics_subset, order(Category, Sample)),]
blueberry_pathcov_acc_metrics_subset$Database <- "PICRUSt2"
blueberry_pathcov_acc_metrics_subset$Database[which(blueberry_pathcov_acc_metrics_subset$Category=="Null")] <- "Null"
blueberry_pathcov_acc_metrics_subset$Category <- factor(blueberry_pathcov_acc_metrics_subset$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_pathcov_acc_metrics_subset_melt <- melt(blueberry_pathcov_acc_metrics_subset)
blueberry_pathcov_acc_metrics_boxplots <- ggplot(blueberry_pathcov_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.2, 1.0)) + ylab(c("")) + ggtitle("Soil (Blueberry)") + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))


# Plot boxplots:
plot_grid(hmp_ec_acc_metrics_boxplots,
          mammal_ec_acc_metrics_boxplots,
          ocean_ec_acc_metrics_boxplots,
          blueberry_ec_acc_metrics_boxplots,
          labels=c("A", "B", "C", "D"))

plot_grid(hmp_pathabun_acc_metrics_boxplots,
          mammal_pathabun_acc_metrics_boxplots,
          ocean_pathabun_acc_metrics_boxplots,
          blueberry_pathabun_acc_metrics_boxplots,
          labels=c("A", "B", "C", "D"))

plot_grid(hmp_pathcov_acc_metrics_boxplots,
          mammal_pathcov_acc_metrics_boxplots,
          ocean_pathcov_acc_metrics_boxplots,
          blueberry_pathcov_acc_metrics_boxplots,
          labels=c("A", "B", "C", "D"))

# Get significance between groups.

hmp_ec_precision_nsti2_vs_null <- wilcox.test(hmp_ec_acc_metrics_subset$Precision[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_ec_acc_metrics_subset$Precision[which(hmp_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

hmp_ec_precision_nsti2_vs_paprica <- wilcox.test(hmp_ec_acc_metrics_subset$Precision[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_ec_acc_metrics_subset$Precision[which(hmp_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

hmp_ec_precision_nsti2_vs_nsti2gg <- wilcox.test(hmp_ec_acc_metrics_subset$Precision[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 hmp_ec_acc_metrics_subset$Precision[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


hmp_ec_recall_nsti2_vs_null <- wilcox.test(hmp_ec_acc_metrics_subset$Recall[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_ec_acc_metrics_subset$Recall[which(hmp_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

hmp_ec_recall_nsti2_vs_paprica <- wilcox.test(hmp_ec_acc_metrics_subset$Recall[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 hmp_ec_acc_metrics_subset$Recall[which(hmp_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

hmp_ec_recall_nsti2_vs_nsti2gg <- wilcox.test(hmp_ec_acc_metrics_subset$Recall[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 hmp_ec_acc_metrics_subset$Recall[which(hmp_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)

hmp_ec_precision_wilcox <- list(nsti2_vs_null=hmp_ec_precision_nsti2_vs_null,
                                nsti2_vs_paprica=hmp_ec_precision_nsti2_vs_paprica,
                                nsti2_vs_nsti2gg=hmp_ec_precision_nsti2_vs_nsti2gg)

hmp_ec_recall_wilcox <- list(nsti2_vs_null=hmp_ec_recall_nsti2_vs_null,
                                nsti2_vs_paprica=hmp_ec_recall_nsti2_vs_paprica,
                                nsti2_vs_nsti2gg=hmp_ec_recall_nsti2_vs_nsti2gg)


hmp_ec_metrics_wilcox <- list(precision=hmp_ec_precision_wilcox,
                              recall=hmp_ec_recall_wilcox)



hmp_pathabun_precision_nsti2_vs_null <- wilcox.test(hmp_pathabun_acc_metrics_subset$Precision[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_pathabun_acc_metrics_subset$Precision[which(hmp_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

hmp_pathabun_precision_nsti2_vs_paprica <- wilcox.test(hmp_pathabun_acc_metrics_subset$Precision[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                 hmp_pathabun_acc_metrics_subset$Precision[which(hmp_pathabun_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

hmp_pathabun_precision_nsti2_vs_nsti2gg <- wilcox.test(hmp_pathabun_acc_metrics_subset$Precision[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                 hmp_pathabun_acc_metrics_subset$Precision[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


hmp_pathabun_recall_nsti2_vs_null <- wilcox.test(hmp_pathabun_acc_metrics_subset$Recall[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                           hmp_pathabun_acc_metrics_subset$Recall[which(hmp_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

hmp_pathabun_recall_nsti2_vs_nsti2gg <- wilcox.test(hmp_pathabun_acc_metrics_subset$Recall[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_pathabun_acc_metrics_subset$Recall[which(hmp_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)

hmp_pathabun_precision_wilcox <- list(nsti2_vs_null=hmp_pathabun_precision_nsti2_vs_null,
                                nsti2_vs_nsti2gg=hmp_pathabun_precision_nsti2_vs_nsti2gg)

hmp_pathabun_recall_wilcox <- list(nsti2_vs_null=hmp_pathabun_recall_nsti2_vs_null,
                             nsti2_vs_nsti2gg=hmp_pathabun_recall_nsti2_vs_nsti2gg)

hmp_pathabun_metrics_wilcox <- list(precision=hmp_pathabun_precision_wilcox,
                              recall=hmp_pathabun_recall_wilcox)



hmp_pathcov_precision_nsti2_vs_null <- wilcox.test(hmp_pathcov_acc_metrics_subset$Precision[which(hmp_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_pathcov_acc_metrics_subset$Precision[which(hmp_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

hmp_pathcov_precision_nsti2_vs_nsti2gg <- wilcox.test(hmp_pathcov_acc_metrics_subset$Precision[which(hmp_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                 hmp_pathcov_acc_metrics_subset$Precision[which(hmp_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


hmp_pathcov_recall_nsti2_vs_null <- wilcox.test(hmp_pathcov_acc_metrics_subset$Recall[which(hmp_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                           hmp_pathcov_acc_metrics_subset$Recall[which(hmp_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

hmp_pathcov_recall_nsti2_vs_nsti2gg <- wilcox.test(hmp_pathcov_acc_metrics_subset$Recall[which(hmp_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                              hmp_pathcov_acc_metrics_subset$Recall[which(hmp_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)
hmp_pathcov_precision_wilcox <- list(nsti2_vs_null=hmp_pathcov_precision_nsti2_vs_null,
                                nsti2_vs_nsti2gg=hmp_pathcov_precision_nsti2_vs_nsti2gg)

hmp_pathcov_recall_wilcox <- list(nsti2_vs_null=hmp_pathcov_recall_nsti2_vs_null,
                             nsti2_vs_nsti2gg=hmp_pathcov_recall_nsti2_vs_nsti2gg)

hmp_pathcov_metrics_wilcox <- list(precision=hmp_pathcov_precision_wilcox,
                              recall=hmp_pathcov_recall_wilcox)


hmp_metrics_wilcox <- list(ec=hmp_ec_metrics_wilcox,
                           pathabun=hmp_pathabun_metrics_wilcox,
                           pathcov=hmp_pathcov_metrics_wilcox)



mammal_ec_precision_nsti2_vs_null <- wilcox.test(mammal_ec_acc_metrics_subset$Precision[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              mammal_ec_acc_metrics_subset$Precision[which(mammal_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

mammal_ec_precision_nsti2_vs_paprica <- wilcox.test(mammal_ec_acc_metrics_subset$Precision[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 mammal_ec_acc_metrics_subset$Precision[which(mammal_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

mammal_ec_precision_nsti2_vs_nsti2gg <- wilcox.test(mammal_ec_acc_metrics_subset$Precision[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 mammal_ec_acc_metrics_subset$Precision[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


mammal_ec_recall_nsti2_vs_null <- wilcox.test(mammal_ec_acc_metrics_subset$Recall[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2")],
                                           mammal_ec_acc_metrics_subset$Recall[which(mammal_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

mammal_ec_recall_nsti2_vs_paprica <- wilcox.test(mammal_ec_acc_metrics_subset$Recall[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              mammal_ec_acc_metrics_subset$Recall[which(mammal_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

mammal_ec_recall_nsti2_vs_nsti2gg <- wilcox.test(mammal_ec_acc_metrics_subset$Recall[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              mammal_ec_acc_metrics_subset$Recall[which(mammal_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)

mammal_ec_precision_wilcox <- list(nsti2_vs_null=mammal_ec_precision_nsti2_vs_null,
                                nsti2_vs_paprica=mammal_ec_precision_nsti2_vs_paprica,
                                nsti2_vs_nsti2gg=mammal_ec_precision_nsti2_vs_nsti2gg)

mammal_ec_recall_wilcox <- list(nsti2_vs_null=mammal_ec_recall_nsti2_vs_null,
                             nsti2_vs_paprica=mammal_ec_recall_nsti2_vs_paprica,
                             nsti2_vs_nsti2gg=mammal_ec_recall_nsti2_vs_nsti2gg)

mammal_ec_metrics_wilcox <- list(precision=mammal_ec_precision_wilcox,
                              recall=mammal_ec_recall_wilcox)



mammal_pathabun_precision_nsti2_vs_null <- wilcox.test(mammal_pathabun_acc_metrics_subset$Precision[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                    mammal_pathabun_acc_metrics_subset$Precision[which(mammal_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

mammal_pathabun_precision_nsti2_vs_paprica <- wilcox.test(mammal_pathabun_acc_metrics_subset$Precision[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                       mammal_pathabun_acc_metrics_subset$Precision[which(mammal_pathabun_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

mammal_pathabun_precision_nsti2_vs_nsti2gg <- wilcox.test(mammal_pathabun_acc_metrics_subset$Precision[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                       mammal_pathabun_acc_metrics_subset$Precision[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


mammal_pathabun_recall_nsti2_vs_null <- wilcox.test(mammal_pathabun_acc_metrics_subset$Recall[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                 mammal_pathabun_acc_metrics_subset$Recall[which(mammal_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

mammal_pathabun_recall_nsti2_vs_nsti2gg <- wilcox.test(mammal_pathabun_acc_metrics_subset$Recall[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                    mammal_pathabun_acc_metrics_subset$Recall[which(mammal_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


mammal_pathabun_precision_wilcox <- list(nsti2_vs_null=mammal_pathabun_precision_nsti2_vs_null,
                                      nsti2_vs_nsti2gg=mammal_pathabun_precision_nsti2_vs_nsti2gg)

mammal_pathabun_recall_wilcox <- list(nsti2_vs_null=mammal_pathabun_recall_nsti2_vs_null,
                                   nsti2_vs_nsti2gg=mammal_pathabun_recall_nsti2_vs_nsti2gg)

mammal_pathabun_metrics_wilcox <- list(precision=mammal_pathabun_precision_wilcox,
                                    recall=mammal_pathabun_recall_wilcox)



mammal_pathcov_precision_nsti2_vs_null <- wilcox.test(mammal_pathcov_acc_metrics_subset$Precision[which(mammal_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                   mammal_pathcov_acc_metrics_subset$Precision[which(mammal_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

mammal_pathcov_precision_nsti2_vs_nsti2gg <- wilcox.test(mammal_pathcov_acc_metrics_subset$Precision[which(mammal_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                      mammal_pathcov_acc_metrics_subset$Precision[which(mammal_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


mammal_pathcov_recall_nsti2_vs_null <- wilcox.test(mammal_pathcov_acc_metrics_subset$Recall[which(mammal_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                mammal_pathcov_acc_metrics_subset$Recall[which(mammal_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

mammal_pathcov_recall_nsti2_vs_nsti2gg <- wilcox.test(mammal_pathcov_acc_metrics_subset$Recall[which(mammal_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                   mammal_pathcov_acc_metrics_subset$Recall[which(mammal_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


mammal_pathcov_precision_wilcox <- list(nsti2_vs_null=mammal_pathcov_precision_nsti2_vs_null,
                                     nsti2_vs_nsti2gg=mammal_pathcov_precision_nsti2_vs_nsti2gg)

mammal_pathcov_recall_wilcox <- list(nsti2_vs_null=mammal_pathcov_recall_nsti2_vs_null,
                                  nsti2_vs_nsti2gg=mammal_pathcov_recall_nsti2_vs_nsti2gg)

mammal_pathcov_metrics_wilcox <- list(precision=mammal_pathcov_precision_wilcox,
                                   recall=mammal_pathcov_recall_wilcox)


mammal_metrics_wilcox <- list(ec=mammal_ec_metrics_wilcox,
                           pathabun=mammal_pathabun_metrics_wilcox,
                           pathcov=mammal_pathcov_metrics_wilcox)



ocean_ec_precision_nsti2_vs_null <- wilcox.test(ocean_ec_acc_metrics_subset$Precision[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              ocean_ec_acc_metrics_subset$Precision[which(ocean_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

ocean_ec_precision_nsti2_vs_paprica <- wilcox.test(ocean_ec_acc_metrics_subset$Precision[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 ocean_ec_acc_metrics_subset$Precision[which(ocean_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

ocean_ec_precision_nsti2_vs_nsti2gg <- wilcox.test(ocean_ec_acc_metrics_subset$Precision[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                 ocean_ec_acc_metrics_subset$Precision[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


ocean_ec_recall_nsti2_vs_null <- wilcox.test(ocean_ec_acc_metrics_subset$Recall[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2")],
                                           ocean_ec_acc_metrics_subset$Recall[which(ocean_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

ocean_ec_recall_nsti2_vs_paprica <- wilcox.test(ocean_ec_acc_metrics_subset$Recall[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              ocean_ec_acc_metrics_subset$Recall[which(ocean_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

ocean_ec_recall_nsti2_vs_nsti2gg <- wilcox.test(ocean_ec_acc_metrics_subset$Recall[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2")],
                                              ocean_ec_acc_metrics_subset$Recall[which(ocean_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


ocean_ec_precision_wilcox <- list(nsti2_vs_null=ocean_ec_precision_nsti2_vs_null,
                                nsti2_vs_paprica=ocean_ec_precision_nsti2_vs_paprica,
                                nsti2_vs_nsti2gg=ocean_ec_precision_nsti2_vs_nsti2gg)

ocean_ec_recall_wilcox <- list(nsti2_vs_null=ocean_ec_recall_nsti2_vs_null,
                             nsti2_vs_paprica=ocean_ec_recall_nsti2_vs_paprica,
                             nsti2_vs_nsti2gg=ocean_ec_recall_nsti2_vs_nsti2gg)


ocean_ec_metrics_wilcox <- list(precision=ocean_ec_precision_wilcox,
                              recall=ocean_ec_recall_wilcox)

ocean_pathabun_precision_nsti2_vs_null <- wilcox.test(ocean_pathabun_acc_metrics_subset$Precision[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                     ocean_pathabun_acc_metrics_subset$Precision[which(ocean_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

ocean_pathabun_precision_nsti2_vs_paprica <- wilcox.test(ocean_pathabun_acc_metrics_subset$Precision[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                        ocean_pathabun_acc_metrics_subset$Precision[which(ocean_pathabun_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

ocean_pathabun_precision_nsti2_vs_nsti2gg <- wilcox.test(ocean_pathabun_acc_metrics_subset$Precision[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                        ocean_pathabun_acc_metrics_subset$Precision[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


ocean_pathabun_recall_nsti2_vs_null <- wilcox.test(ocean_pathabun_acc_metrics_subset$Recall[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                  ocean_pathabun_acc_metrics_subset$Recall[which(ocean_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

ocean_pathabun_recall_nsti2_vs_nsti2gg <- wilcox.test(ocean_pathabun_acc_metrics_subset$Recall[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                     ocean_pathabun_acc_metrics_subset$Recall[which(ocean_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


ocean_pathabun_precision_wilcox <- list(nsti2_vs_null=ocean_pathabun_precision_nsti2_vs_null,
                                       nsti2_vs_nsti2gg=ocean_pathabun_precision_nsti2_vs_nsti2gg)

ocean_pathabun_recall_wilcox <- list(nsti2_vs_null=ocean_pathabun_recall_nsti2_vs_null,
                                    nsti2_vs_nsti2gg=ocean_pathabun_recall_nsti2_vs_nsti2gg)


ocean_pathabun_metrics_wilcox <- list(precision=ocean_pathabun_precision_wilcox,
                                     recall=ocean_pathabun_recall_wilcox)



ocean_pathcov_precision_nsti2_vs_null <- wilcox.test(ocean_pathcov_acc_metrics_subset$Precision[which(ocean_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                    ocean_pathcov_acc_metrics_subset$Precision[which(ocean_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

ocean_pathcov_precision_nsti2_vs_nsti2gg <- wilcox.test(ocean_pathcov_acc_metrics_subset$Precision[which(ocean_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                       ocean_pathcov_acc_metrics_subset$Precision[which(ocean_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


ocean_pathcov_recall_nsti2_vs_null <- wilcox.test(ocean_pathcov_acc_metrics_subset$Recall[which(ocean_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                 ocean_pathcov_acc_metrics_subset$Recall[which(ocean_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

ocean_pathcov_recall_nsti2_vs_nsti2gg <- wilcox.test(ocean_pathcov_acc_metrics_subset$Recall[which(ocean_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                    ocean_pathcov_acc_metrics_subset$Recall[which(ocean_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


ocean_pathcov_precision_wilcox <- list(nsti2_vs_null=ocean_pathcov_precision_nsti2_vs_null,
                                      nsti2_vs_nsti2gg=ocean_pathcov_precision_nsti2_vs_nsti2gg)

ocean_pathcov_recall_wilcox <- list(nsti2_vs_null=ocean_pathcov_recall_nsti2_vs_null,
                                   nsti2_vs_nsti2gg=ocean_pathcov_recall_nsti2_vs_nsti2gg)

ocean_pathcov_metrics_wilcox <- list(precision=ocean_pathcov_precision_wilcox,
                                    recall=ocean_pathcov_recall_wilcox)

ocean_metrics_wilcox <- list(ec=ocean_ec_metrics_wilcox,
                            pathabun=ocean_pathabun_metrics_wilcox,
                            pathcov=ocean_pathcov_metrics_wilcox)


blueberry_ec_precision_nsti2_vs_null <- wilcox.test(blueberry_ec_acc_metrics_subset$Precision[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                blueberry_ec_acc_metrics_subset$Precision[which(blueberry_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

blueberry_ec_precision_nsti2_vs_paprica <- wilcox.test(blueberry_ec_acc_metrics_subset$Precision[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                   blueberry_ec_acc_metrics_subset$Precision[which(blueberry_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

blueberry_ec_precision_nsti2_vs_nsti2gg <- wilcox.test(blueberry_ec_acc_metrics_subset$Precision[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                   blueberry_ec_acc_metrics_subset$Precision[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


blueberry_ec_recall_nsti2_vs_null <- wilcox.test(blueberry_ec_acc_metrics_subset$Recall[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2")],
                                             blueberry_ec_acc_metrics_subset$Recall[which(blueberry_ec_acc_metrics_subset$Category=="Null")], paired=TRUE)

blueberry_ec_recall_nsti2_vs_paprica <- wilcox.test(blueberry_ec_acc_metrics_subset$Recall[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                blueberry_ec_acc_metrics_subset$Recall[which(blueberry_ec_acc_metrics_subset$Category=="PAPRICA")], paired=TRUE)

blueberry_ec_recall_nsti2_vs_nsti2gg <- wilcox.test(blueberry_ec_acc_metrics_subset$Recall[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2")],
                                                blueberry_ec_acc_metrics_subset$Recall[which(blueberry_ec_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


blueberry_ec_precision_wilcox <- list(nsti2_vs_null=blueberry_ec_precision_nsti2_vs_null,
                                  nsti2_vs_paprica=blueberry_ec_precision_nsti2_vs_paprica,
                                  nsti2_vs_nsti2gg=blueberry_ec_precision_nsti2_vs_nsti2gg)

blueberry_ec_recall_wilcox <- list(nsti2_vs_null=blueberry_ec_recall_nsti2_vs_null,
                               nsti2_vs_paprica=blueberry_ec_recall_nsti2_vs_paprica,
                               nsti2_vs_nsti2gg=blueberry_ec_recall_nsti2_vs_nsti2gg)

blueberry_ec_metrics_wilcox <- list(precision=blueberry_ec_precision_wilcox,
                                recall=blueberry_ec_recall_wilcox)


blueberry_pathabun_precision_nsti2_vs_null <- wilcox.test(blueberry_pathabun_acc_metrics_subset$Precision[which(blueberry_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                    blueberry_pathabun_acc_metrics_subset$Precision[which(blueberry_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

blueberry_pathabun_precision_nsti2_vs_nsti2gg <- wilcox.test(blueberry_pathabun_acc_metrics_subset$Precision[which(blueberry_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                       blueberry_pathabun_acc_metrics_subset$Precision[which(blueberry_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


blueberry_pathabun_recall_nsti2_vs_null <- wilcox.test(blueberry_pathabun_acc_metrics_subset$Recall[which(blueberry_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                 blueberry_pathabun_acc_metrics_subset$Recall[which(blueberry_pathabun_acc_metrics_subset$Category=="Null")], paired=TRUE)

blueberry_pathabun_recall_nsti2_vs_nsti2gg <- wilcox.test(blueberry_pathabun_acc_metrics_subset$Recall[which(blueberry_pathabun_acc_metrics_subset$Category=="NSTI=2")],
                                                    blueberry_pathabun_acc_metrics_subset$Recall[which(blueberry_pathabun_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_pathabun_precision_wilcox <- list(nsti2_vs_null=blueberry_pathabun_precision_nsti2_vs_null,
                                      nsti2_vs_nsti2gg=blueberry_pathabun_precision_nsti2_vs_nsti2gg)

blueberry_pathabun_recall_wilcox <- list(nsti2_vs_null=blueberry_pathabun_recall_nsti2_vs_null,
                                   nsti2_vs_nsti2gg=blueberry_pathabun_recall_nsti2_vs_nsti2gg)

blueberry_pathabun_metrics_wilcox <- list(precision=blueberry_pathabun_precision_wilcox,
                                    recall=blueberry_pathabun_recall_wilcox)



blueberry_pathcov_precision_nsti2_vs_null <- wilcox.test(blueberry_pathcov_acc_metrics_subset$Precision[which(blueberry_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                   blueberry_pathcov_acc_metrics_subset$Precision[which(blueberry_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

blueberry_pathcov_precision_nsti2_vs_nsti2gg <- wilcox.test(blueberry_pathcov_acc_metrics_subset$Precision[which(blueberry_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                      blueberry_pathcov_acc_metrics_subset$Precision[which(blueberry_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)


blueberry_pathcov_recall_nsti2_vs_null <- wilcox.test(blueberry_pathcov_acc_metrics_subset$Recall[which(blueberry_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                blueberry_pathcov_acc_metrics_subset$Recall[which(blueberry_pathcov_acc_metrics_subset$Category=="Null")], paired=TRUE)

blueberry_pathcov_recall_nsti2_vs_nsti2gg <- wilcox.test(blueberry_pathcov_acc_metrics_subset$Recall[which(blueberry_pathcov_acc_metrics_subset$Category=="NSTI=2")],
                                                   blueberry_pathcov_acc_metrics_subset$Recall[which(blueberry_pathcov_acc_metrics_subset$Category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_pathcov_precision_wilcox <- list(nsti2_vs_null=blueberry_pathcov_precision_nsti2_vs_null,
                                     nsti2_vs_nsti2gg=blueberry_pathcov_precision_nsti2_vs_nsti2gg)

blueberry_pathcov_recall_wilcox <- list(nsti2_vs_null=blueberry_pathcov_recall_nsti2_vs_null,
                                  nsti2_vs_nsti2gg=blueberry_pathcov_recall_nsti2_vs_nsti2gg)

blueberry_pathcov_metrics_wilcox <- list(precision=blueberry_pathcov_precision_wilcox,
                                   recall=blueberry_pathcov_recall_wilcox)

blueberry_metrics_wilcox <- list(ec=blueberry_ec_metrics_wilcox,
                           pathabun=blueberry_pathabun_metrics_wilcox,
                           pathcov=blueberry_pathcov_metrics_wilcox)

combined_metrics_wilcox <- list(hmp=hmp_metrics_wilcox,
                                mammal=mammal_metrics_wilcox,
                                ocean=ocean_metrics_wilcox,
                                blueberry=blueberry_metrics_wilcox)

saveRDS(object = combined_metrics_wilcox, file = "combined_metrics_wilcoxon_ec_path.rds")

### Read in when needed: ###
combined_metrics_wilcox <- readRDS("combined_metrics_wilcoxon_ec_path.rds")
sapply(combined_metrics_wilcox$hmp$ec$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$hmp$ec$recall, function(x){x$p.value})

sapply(combined_metrics_wilcox$mammal$ec$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$mammal$ec$recall, function(x){x$p.value})


sapply(combined_metrics_wilcox$ocean$ec$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$ocean$ec$recall, function(x){x$p.value})


sapply(combined_metrics_wilcox$blueberry$ec$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$blueberry$ec$recall, function(x){x$p.value})



sapply(combined_metrics_wilcox$hmp$pathabun$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$hmp$pathabun$recall, function(x){x$p.value})

sapply(combined_metrics_wilcox$mammal$pathabun$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$mammal$pathabun$recall, function(x){x$p.value})


sapply(combined_metrics_wilcox$ocean$pathabun$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$ocean$pathabun$recall, function(x){x$p.value})


sapply(combined_metrics_wilcox$blueberry$pathabun$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$blueberry$pathabun$recall, function(x){x$p.value})





sapply(combined_metrics_wilcox$hmp$pathcov$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$hmp$pathcov$recall, function(x){x$p.value})

sapply(combined_metrics_wilcox$mammal$pathcov$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$mammal$pathcov$recall, function(x){x$p.value})


sapply(combined_metrics_wilcox$ocean$pathcov$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$ocean$pathcov$recall, function(x){x$p.value})


sapply(combined_metrics_wilcox$blueberry$pathcov$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$blueberry$pathcov$recall, function(x){x$p.value})



# Compare PAPRICA and null:
wilcox.test(hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "PAPRICA"), "precision"],
            hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "precision"])

wilcox.test(mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "PAPRICA"), "precision"],
            mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "precision"])

wilcox.test(ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "PAPRICA"), "precision"],
            ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "precision"])

wilcox.test(blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "PAPRICA"), "precision"],
            blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "precision"])


# What % increase in recall over each dataset?
mean(((hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "NSTI=2"), "recall"] - hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "recall"])/hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "recall"])*100)
mean(((mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "NSTI=2"), "recall"] - mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "recall"])/mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "recall"])*100)
mean(((ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "NSTI=2"), "recall"] - ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "recall"])/ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "recall"])*100)
mean(((blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "NSTI=2"), "recall"] - blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "recall"])/blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "recall"])*100)


sd(((hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "NSTI=2"), "recall"] - hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "recall"])/hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "recall"])*100)
sd(((mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "NSTI=2"), "recall"] - mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "recall"])/mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "recall"])*100)
sd(((ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "NSTI=2"), "recall"] - ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "recall"])/ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "recall"])*100)
sd(((blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "NSTI=2"), "recall"] - blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "recall"])/blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "recall"])*100)


mean(((hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "NSTI=2"), "precision"] - hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "precision"])/hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "precision"])*100)
mean(((mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "NSTI=2"), "precision"] - mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "precision"])/mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "precision"])*100)
mean(((ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "NSTI=2"), "precision"] - ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "precision"])/ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "precision"])*100)
mean(((blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "NSTI=2"), "precision"] - blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "precision"])/blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "precision"])*100)


sd(((hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "NSTI=2"), "precision"] - hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "precision"])/hmp_ec_acc_metrics[which(hmp_ec_acc_metrics$category == "Null"), "precision"])*100)
sd(((mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "NSTI=2"), "precision"] - mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "precision"])/mammal_ec_acc_metrics[which(mammal_ec_acc_metrics$category == "Null"), "precision"])*100)
sd(((ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "NSTI=2"), "precision"] - ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "precision"])/ocean_ec_acc_metrics[which(ocean_ec_acc_metrics$category == "Null"), "precision"])*100)
sd(((blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "NSTI=2"), "precision"] - blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "precision"])/blueberry_ec_acc_metrics[which(blueberry_ec_acc_metrics$category == "Null"), "precision"])*100)


