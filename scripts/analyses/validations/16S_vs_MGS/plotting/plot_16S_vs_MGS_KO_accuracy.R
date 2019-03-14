### Make figure contrasting accuracy metrics across HSP tools on each 16S validation dataset.
### Also test for statistical significance between these categories (and save wilcoxon output to RDS).

library(ggplot2)
library(reshape2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

hmp_acc_metrics <- readRDS("hmp_acc_metrics.rds")

hmp_acc_metrics_subset <- hmp_acc_metrics[,c("sample", "precision", "recall", "category")]

colnames(hmp_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")

hmp_acc_metrics_subset$Database <- "Other"

hmp_acc_metrics_subset$Database[grep("NSTI" , hmp_acc_metrics_subset$Category)] <- "PICRUSt2"

hmp_acc_metrics_subset$Database[which(hmp_acc_metrics_subset$Category=="Null")] <- "Null"

hmp_acc_metrics_subset$Category <- factor(hmp_acc_metrics_subset$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

hmp_acc_metrics_subset_melt <- melt(hmp_acc_metrics_subset)

hmp_acc_metrics_boxplots <- ggplot(hmp_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
                                   facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
                                   ylim(c(0.25, 1.0)) + ylab(c("")) + ggtitle("HMP")  + guides(fill=FALSE) +
                                   scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


mammal_acc_metrics <- readRDS("mammal_acc_metrics.rds")

mammal_acc_metrics_subset <- mammal_acc_metrics[,c("sample", "precision", "recall", "category")]

colnames(mammal_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")

mammal_acc_metrics_subset$Database <- "Other"

mammal_acc_metrics_subset$Database[grep("NSTI" , mammal_acc_metrics_subset$Category)] <- "PICRUSt2"

mammal_acc_metrics_subset$Database[which(mammal_acc_metrics_subset$Category=="Null")] <- "Null"

mammal_acc_metrics_subset$Category <- factor(mammal_acc_metrics_subset$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

mammal_acc_metrics_subset_melt <- melt(mammal_acc_metrics_subset)

mammal_acc_metrics_boxplots <- ggplot(mammal_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.25, 1.0)) + ylab(c("")) + ggtitle("Mammal") + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


ocean_acc_metrics <- readRDS("ocean_acc_metrics.rds")

ocean_acc_metrics_subset <- ocean_acc_metrics[,c("sample", "precision", "recall", "category")]

colnames(ocean_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")

ocean_acc_metrics_subset$Database <- "Other"

ocean_acc_metrics_subset$Database[grep("NSTI" , ocean_acc_metrics_subset$Category)] <- "PICRUSt2"

ocean_acc_metrics_subset$Database[which(ocean_acc_metrics_subset$Category=="Null")] <- "Null"

ocean_acc_metrics_subset$Category <- factor(ocean_acc_metrics_subset$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

ocean_acc_metrics_subset_melt <- melt(ocean_acc_metrics_subset)

ocean_acc_metrics_boxplots <- ggplot(ocean_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.25, 1.0)) + ylab(c("")) + ggtitle("Ocean") + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


blueberry_acc_metrics <- readRDS("blueberry_acc_metrics.rds")

blueberry_acc_metrics_subset <- blueberry_acc_metrics[,c("sample", "precision", "recall", "category")]

colnames(blueberry_acc_metrics_subset) <- c("Sample", "Precision", "Recall", "Category")

blueberry_acc_metrics_subset$Database <- "Other"

blueberry_acc_metrics_subset$Database[grep("NSTI" , blueberry_acc_metrics_subset$Category)] <- "PICRUSt2"

blueberry_acc_metrics_subset$Database[which(blueberry_acc_metrics_subset$Category=="Null")] <- "Null"

blueberry_acc_metrics_subset$Category <- factor(blueberry_acc_metrics_subset$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))

blueberry_acc_metrics_subset_melt <- melt(blueberry_acc_metrics_subset)

blueberry_acc_metrics_boxplots <- ggplot(blueberry_acc_metrics_subset_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() +
  ylim(c(0.25, 1.0)) + ylab(c("")) + ggtitle("Soil (Blueberry)") + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


# Plot boxplots:
plot_grid(hmp_acc_metrics_boxplots,
          mammal_acc_metrics_boxplots,
          ocean_acc_metrics_boxplots,
          blueberry_acc_metrics_boxplots,
          labels=c("A", "B", "C", "D"))

# Get significance between groups.
hmp_acc_metrics_subset_sorted <- hmp_acc_metrics_subset[
  with(hmp_acc_metrics_subset, order(Category, Sample)),]

hmp_precision_nsti2_vs_null <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2")],
                                           hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="Null")], paired=TRUE)

hmp_precision_nsti2_vs_tax4fun <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2")],
                                           hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="Tax4Fun")], paired=TRUE)

hmp_precision_nsti2_vs_panfp <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2")],
                                           hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="PanFP")], paired=TRUE)

hmp_precision_nsti2_vs_piphillin <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2")],
            hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="Piphillin")], paired=TRUE)

hmp_precision_nsti2_vs_picrust1 <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2")],
            hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="PICRUSt1")], paired=TRUE)

hmp_precision_nsti2_vs_nsti2gg <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2")],
            hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)


hmp_recall_nsti2_vs_null <- wilcox.test(hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2")],
                                           hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="Null")], paired=TRUE)

hmp_recall_nsti2_vs_tax4fun <- wilcox.test(hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2")],
                                              hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="Tax4Fun")], paired=TRUE)

hmp_recall_nsti2_vs_panfp <- wilcox.test(hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2")],
                                            hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="PanFP")], paired=TRUE)

hmp_recall_nsti2_vs_piphillin <- wilcox.test(hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2")],
                                                hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="Piphillin")], paired=TRUE)

hmp_recall_nsti2_vs_picrust1 <- wilcox.test(hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2")],
                                               hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="PICRUSt1")], paired=TRUE)

hmp_recall_nsti2_vs_nsti2gg <- wilcox.test(hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2")],
                                              hmp_acc_metrics$recall[which(hmp_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

hmp_fpr_nsti2_vs_null <- wilcox.test(hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2")],
                                           hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="Null")], paired=TRUE)

hmp_fpr_nsti2_vs_tax4fun <- wilcox.test(hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2")],
                                              hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="Tax4Fun")], paired=TRUE)

hmp_fpr_nsti2_vs_panfp <- wilcox.test(hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2")],
                                            hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="PanFP")], paired=TRUE)

hmp_fpr_nsti2_vs_piphillin <- wilcox.test(hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2")],
                                                hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="Piphillin")], paired=TRUE)

hmp_fpr_nsti2_vs_picrust1 <- wilcox.test(hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2")],
                                               hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="PICRUSt1")], paired=TRUE)

hmp_fpr_nsti2_vs_nsti2gg <- wilcox.test(hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2")],
                                              hmp_acc_metrics$fpr[which(hmp_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

hmp_precision_wilcox <- list(nsti2_vs_null=hmp_precision_nsti2_vs_null,
                             nsti2_vs_tax4fun=hmp_precision_nsti2_vs_tax4fun,
                             nsti2_vs_panfp=hmp_precision_nsti2_vs_panfp,
                             nsti2_vs_piphillin=hmp_precision_nsti2_vs_piphillin,
                             nsti2_vs_picrust1=hmp_precision_nsti2_vs_picrust1,
                             nsti2_vs_nsti2gg=hmp_precision_nsti2_vs_nsti2gg)

hmp_recall_wilcox <- list(nsti2_vs_null=hmp_recall_nsti2_vs_null,
                             nsti2_vs_tax4fun=hmp_recall_nsti2_vs_tax4fun,
                             nsti2_vs_panfp=hmp_recall_nsti2_vs_panfp,
                             nsti2_vs_piphillin=hmp_recall_nsti2_vs_piphillin,
                             nsti2_vs_picrust1=hmp_recall_nsti2_vs_picrust1,
                             nsti2_vs_nsti2gg=hmp_recall_nsti2_vs_nsti2gg)

hmp_fpr_wilcox <- list(nsti2_vs_null=hmp_fpr_nsti2_vs_null,
                             nsti2_vs_tax4fun=hmp_fpr_nsti2_vs_tax4fun,
                             nsti2_vs_panfp=hmp_fpr_nsti2_vs_panfp,
                             nsti2_vs_piphillin=hmp_fpr_nsti2_vs_piphillin,
                             nsti2_vs_picrust1=hmp_fpr_nsti2_vs_picrust1,
                             nsti2_vs_nsti2gg=hmp_fpr_nsti2_vs_nsti2gg)

hmp_metrics_wilcox <- list(precision=hmp_precision_wilcox, recall=hmp_recall_wilcox, fpr=hmp_fpr_wilcox)


mammal_acc_metrics_subset_sorted <- mammal_acc_metrics_subset[
  with(mammal_acc_metrics_subset, order(Category, Sample)),]

mammal_precision_nsti2_vs_null <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2")],
                                           mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="Null")], paired=TRUE)

mammal_precision_nsti2_vs_tax4fun <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2")],
                                              mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="Tax4Fun")], paired=TRUE)

mammal_precision_nsti2_vs_panfp <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2")],
                                            mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="PanFP")], paired=TRUE)

mammal_precision_nsti2_vs_piphillin <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2")],
                                                mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="Piphillin")], paired=TRUE)

mammal_precision_nsti2_vs_picrust1 <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2")],
                                               mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="PICRUSt1")], paired=TRUE)

mammal_precision_nsti2_vs_nsti2gg <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2")],
                                              mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)


mammal_recall_nsti2_vs_null <- wilcox.test(mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2")],
                                        mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="Null")], paired=TRUE)

mammal_recall_nsti2_vs_tax4fun <- wilcox.test(mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2")],
                                           mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="Tax4Fun")], paired=TRUE)

mammal_recall_nsti2_vs_panfp <- wilcox.test(mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2")],
                                         mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="PanFP")], paired=TRUE)

mammal_recall_nsti2_vs_piphillin <- wilcox.test(mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2")],
                                             mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="Piphillin")], paired=TRUE)

mammal_recall_nsti2_vs_picrust1 <- wilcox.test(mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2")],
                                            mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="PICRUSt1")], paired=TRUE)

mammal_recall_nsti2_vs_nsti2gg <- wilcox.test(mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2")],
                                           mammal_acc_metrics$recall[which(mammal_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

mammal_fpr_nsti2_vs_null <- wilcox.test(mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2")],
                                     mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="Null")], paired=TRUE)

mammal_fpr_nsti2_vs_tax4fun <- wilcox.test(mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2")],
                                        mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="Tax4Fun")], paired=TRUE)

mammal_fpr_nsti2_vs_panfp <- wilcox.test(mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2")],
                                      mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="PanFP")], paired=TRUE)

mammal_fpr_nsti2_vs_piphillin <- wilcox.test(mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2")],
                                          mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="Piphillin")], paired=TRUE)

mammal_fpr_nsti2_vs_picrust1 <- wilcox.test(mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2")],
                                         mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="PICRUSt1")], paired=TRUE)

mammal_fpr_nsti2_vs_nsti2gg <- wilcox.test(mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2")],
                                        mammal_acc_metrics$fpr[which(mammal_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

mammal_precision_wilcox <- list(nsti2_vs_null=mammal_precision_nsti2_vs_null,
                             nsti2_vs_tax4fun=mammal_precision_nsti2_vs_tax4fun,
                             nsti2_vs_panfp=mammal_precision_nsti2_vs_panfp,
                             nsti2_vs_piphillin=mammal_precision_nsti2_vs_piphillin,
                             nsti2_vs_picrust1=mammal_precision_nsti2_vs_picrust1,
                             nsti2_vs_nsti2gg=mammal_precision_nsti2_vs_nsti2gg)

mammal_recall_wilcox <- list(nsti2_vs_null=mammal_recall_nsti2_vs_null,
                          nsti2_vs_tax4fun=mammal_recall_nsti2_vs_tax4fun,
                          nsti2_vs_panfp=mammal_recall_nsti2_vs_panfp,
                          nsti2_vs_piphillin=mammal_recall_nsti2_vs_piphillin,
                          nsti2_vs_picrust1=mammal_recall_nsti2_vs_picrust1,
                          nsti2_vs_nsti2gg=mammal_recall_nsti2_vs_nsti2gg)

mammal_fpr_wilcox <- list(nsti2_vs_null=mammal_fpr_nsti2_vs_null,
                       nsti2_vs_tax4fun=mammal_fpr_nsti2_vs_tax4fun,
                       nsti2_vs_panfp=mammal_fpr_nsti2_vs_panfp,
                       nsti2_vs_piphillin=mammal_fpr_nsti2_vs_piphillin,
                       nsti2_vs_picrust1=mammal_fpr_nsti2_vs_picrust1,
                       nsti2_vs_nsti2gg=mammal_fpr_nsti2_vs_nsti2gg)

mammal_metrics_wilcox <- list(precision=mammal_precision_wilcox, recall=mammal_recall_wilcox, fpr=mammal_fpr_wilcox)


ocean_acc_metrics_subset_sorted <- ocean_acc_metrics_subset[
  with(ocean_acc_metrics_subset, order(Category, Sample)),]

ocean_precision_nsti2_vs_null <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2")],
                                           ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="Null")], paired=TRUE)

ocean_precision_nsti2_vs_tax4fun <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2")],
                                              ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="Tax4Fun")], paired=TRUE)

ocean_precision_nsti2_vs_panfp <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2")],
                                            ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="PanFP")], paired=TRUE)

ocean_precision_nsti2_vs_piphillin <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2")],
                                                ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="Piphillin")], paired=TRUE)

ocean_precision_nsti2_vs_picrust1 <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2")],
                                               ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="PICRUSt1")], paired=TRUE)

ocean_precision_nsti2_vs_nsti2gg <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2")],
                                              ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)


ocean_recall_nsti2_vs_null <- wilcox.test(ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2")],
                                        ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="Null")], paired=TRUE)

ocean_recall_nsti2_vs_tax4fun <- wilcox.test(ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2")],
                                           ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="Tax4Fun")], paired=TRUE)

ocean_recall_nsti2_vs_panfp <- wilcox.test(ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2")],
                                         ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="PanFP")], paired=TRUE)

ocean_recall_nsti2_vs_piphillin <- wilcox.test(ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2")],
                                             ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="Piphillin")], paired=TRUE)

ocean_recall_nsti2_vs_picrust1 <- wilcox.test(ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2")],
                                            ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="PICRUSt1")], paired=TRUE)

ocean_recall_nsti2_vs_nsti2gg <- wilcox.test(ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2")],
                                           ocean_acc_metrics$recall[which(ocean_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

ocean_fpr_nsti2_vs_null <- wilcox.test(ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2")],
                                     ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="Null")], paired=TRUE)

ocean_fpr_nsti2_vs_tax4fun <- wilcox.test(ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2")],
                                        ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="Tax4Fun")], paired=TRUE)

ocean_fpr_nsti2_vs_panfp <- wilcox.test(ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2")],
                                      ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="PanFP")], paired=TRUE)

ocean_fpr_nsti2_vs_piphillin <- wilcox.test(ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2")],
                                          ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="Piphillin")], paired=TRUE)

ocean_fpr_nsti2_vs_picrust1 <- wilcox.test(ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2")],
                                         ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="PICRUSt1")], paired=TRUE)

ocean_fpr_nsti2_vs_nsti2gg <- wilcox.test(ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2")],
                                        ocean_acc_metrics$fpr[which(ocean_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

ocean_precision_wilcox <- list(nsti2_vs_null=ocean_precision_nsti2_vs_null,
                             nsti2_vs_tax4fun=ocean_precision_nsti2_vs_tax4fun,
                             nsti2_vs_panfp=ocean_precision_nsti2_vs_panfp,
                             nsti2_vs_piphillin=ocean_precision_nsti2_vs_piphillin,
                             nsti2_vs_picrust1=ocean_precision_nsti2_vs_picrust1,
                             nsti2_vs_nsti2gg=ocean_precision_nsti2_vs_nsti2gg)

ocean_recall_wilcox <- list(nsti2_vs_null=ocean_recall_nsti2_vs_null,
                          nsti2_vs_tax4fun=ocean_recall_nsti2_vs_tax4fun,
                          nsti2_vs_panfp=ocean_recall_nsti2_vs_panfp,
                          nsti2_vs_piphillin=ocean_recall_nsti2_vs_piphillin,
                          nsti2_vs_picrust1=ocean_recall_nsti2_vs_picrust1,
                          nsti2_vs_nsti2gg=ocean_recall_nsti2_vs_nsti2gg)

ocean_fpr_wilcox <- list(nsti2_vs_null=ocean_fpr_nsti2_vs_null,
                       nsti2_vs_tax4fun=ocean_fpr_nsti2_vs_tax4fun,
                       nsti2_vs_panfp=ocean_fpr_nsti2_vs_panfp,
                       nsti2_vs_piphillin=ocean_fpr_nsti2_vs_piphillin,
                       nsti2_vs_picrust1=ocean_fpr_nsti2_vs_picrust1,
                       nsti2_vs_nsti2gg=ocean_fpr_nsti2_vs_nsti2gg)

ocean_metrics_wilcox <- list(precision=ocean_precision_wilcox, recall=ocean_recall_wilcox, fpr=ocean_fpr_wilcox)


blueberry_acc_metrics_subset_sorted <- blueberry_acc_metrics_subset[
  with(blueberry_acc_metrics_subset, order(Category, Sample)),]

blueberry_precision_nsti2_vs_null <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2")],
                                           blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="Null")], paired=TRUE)

blueberry_precision_nsti2_vs_tax4fun <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2")],
                                              blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="Tax4Fun")], paired=TRUE)

blueberry_precision_nsti2_vs_panfp <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2")],
                                            blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="PanFP")], paired=TRUE)

blueberry_precision_nsti2_vs_piphillin <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2")],
                                                blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="Piphillin")], paired=TRUE)

blueberry_precision_nsti2_vs_picrust1 <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2")],
                                               blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="PICRUSt1")], paired=TRUE)

blueberry_precision_nsti2_vs_nsti2gg <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2")],
                                              blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)


blueberry_recall_nsti2_vs_null <- wilcox.test(blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2")],
                                        blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="Null")], paired=TRUE)

blueberry_recall_nsti2_vs_tax4fun <- wilcox.test(blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2")],
                                           blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="Tax4Fun")], paired=TRUE)

blueberry_recall_nsti2_vs_panfp <- wilcox.test(blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2")],
                                         blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="PanFP")], paired=TRUE)

blueberry_recall_nsti2_vs_piphillin <- wilcox.test(blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2")],
                                             blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="Piphillin")], paired=TRUE)

blueberry_recall_nsti2_vs_picrust1 <- wilcox.test(blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2")],
                                            blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="PICRUSt1")], paired=TRUE)

blueberry_recall_nsti2_vs_nsti2gg <- wilcox.test(blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2")],
                                           blueberry_acc_metrics$recall[which(blueberry_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_fpr_nsti2_vs_null <- wilcox.test(blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2")],
                                     blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="Null")], paired=TRUE)

blueberry_fpr_nsti2_vs_tax4fun <- wilcox.test(blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2")],
                                        blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="Tax4Fun")], paired=TRUE)

blueberry_fpr_nsti2_vs_panfp <- wilcox.test(blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2")],
                                      blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="PanFP")], paired=TRUE)

blueberry_fpr_nsti2_vs_piphillin <- wilcox.test(blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2")],
                                          blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="Piphillin")], paired=TRUE)

blueberry_fpr_nsti2_vs_picrust1 <- wilcox.test(blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2")],
                                         blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="PICRUSt1")], paired=TRUE)

blueberry_fpr_nsti2_vs_nsti2gg <- wilcox.test(blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2")],
                                        blueberry_acc_metrics$fpr[which(blueberry_acc_metrics$category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_precision_wilcox <- list(nsti2_vs_null=blueberry_precision_nsti2_vs_null,
                             nsti2_vs_tax4fun=blueberry_precision_nsti2_vs_tax4fun,
                             nsti2_vs_panfp=blueberry_precision_nsti2_vs_panfp,
                             nsti2_vs_piphillin=blueberry_precision_nsti2_vs_piphillin,
                             nsti2_vs_picrust1=blueberry_precision_nsti2_vs_picrust1,
                             nsti2_vs_nsti2gg=blueberry_precision_nsti2_vs_nsti2gg)

blueberry_recall_wilcox <- list(nsti2_vs_null=blueberry_recall_nsti2_vs_null,
                          nsti2_vs_tax4fun=blueberry_recall_nsti2_vs_tax4fun,
                          nsti2_vs_panfp=blueberry_recall_nsti2_vs_panfp,
                          nsti2_vs_piphillin=blueberry_recall_nsti2_vs_piphillin,
                          nsti2_vs_picrust1=blueberry_recall_nsti2_vs_picrust1,
                          nsti2_vs_nsti2gg=blueberry_recall_nsti2_vs_nsti2gg)

blueberry_fpr_wilcox <- list(nsti2_vs_null=blueberry_fpr_nsti2_vs_null,
                       nsti2_vs_tax4fun=blueberry_fpr_nsti2_vs_tax4fun,
                       nsti2_vs_panfp=blueberry_fpr_nsti2_vs_panfp,
                       nsti2_vs_piphillin=blueberry_fpr_nsti2_vs_piphillin,
                       nsti2_vs_picrust1=blueberry_fpr_nsti2_vs_picrust1,
                       nsti2_vs_nsti2gg=blueberry_fpr_nsti2_vs_nsti2gg)

blueberry_metrics_wilcox <- list(precision=blueberry_precision_wilcox, recall=blueberry_recall_wilcox, fpr=blueberry_fpr_wilcox)

combined_metrics_wilcox <- list(hmp=hmp_metrics_wilcox,
                               mammal=mammal_metrics_wilcox,
                               ocean=ocean_metrics_wilcox,
                               blueberry=blueberry_metrics_wilcox)

saveRDS(object = combined_metrics_wilcox, file = "ko_combined_metrics_wilcoxon.rds")

 setwd("/home/gavin/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")
 combined_metrics_wilcox <- readRDS("ko_combined_metrics_wilcoxon.rds")

sapply(combined_metrics_wilcox$hmp$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$hmp$recall, function(x){x$p.value})

sapply(combined_metrics_wilcox$mammal$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$mammal$recall, function(x){x$p.value})

sapply(combined_metrics_wilcox$ocean$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$ocean$recall, function(x){x$p.value})

sapply(combined_metrics_wilcox$blueberry$precision, function(x){x$p.value})
sapply(combined_metrics_wilcox$blueberry$recall, function(x){x$p.value})

# Test for NSTI=0.05 vs PICRUSt1
hmp_nsti0.05_vs_picrust1_wilcox <- wilcox.test(hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="NSTI=0.05")],
                                                  hmp_acc_metrics$precision[which(hmp_acc_metrics$category=="PICRUSt1")], paired=TRUE)


mammal_nsti0.05_vs_picrust1_wilcox <- wilcox.test(mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="NSTI=0.05")],
                                                  mammal_acc_metrics$precision[which(mammal_acc_metrics$category=="PICRUSt1")], paired=TRUE)

ocean_nsti0.05_vs_picrust1_wilcox <- wilcox.test(ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="NSTI=0.05")],
                                                  ocean_acc_metrics$precision[which(ocean_acc_metrics$category=="PICRUSt1")], paired=TRUE)

blueberry_nsti0.05_vs_picrust1_wilcox <- wilcox.test(blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="NSTI=0.05")],
                                                  blueberry_acc_metrics$precision[which(blueberry_acc_metrics$category=="PICRUSt1")], paired=TRUE)


# Mean and sd precision per dataset.
mean(hmp_acc_metrics[which(hmp_acc_metrics$category == "NSTI=2"), "precision"])
sd(hmp_acc_metrics[which(hmp_acc_metrics$category == "NSTI=2"), "precision"])

mean(mammal_acc_metrics[which(mammal_acc_metrics$category == "NSTI=2"), "precision"])
sd(mammal_acc_metrics[which(mammal_acc_metrics$category == "NSTI=2"), "precision"])

mean(ocean_acc_metrics[which(ocean_acc_metrics$category == "NSTI=2"), "precision"])
sd(ocean_acc_metrics[which(ocean_acc_metrics$category == "NSTI=2"), "precision"])

mean(blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "precision"])
sd(blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "precision"])

# What % increase in recall over each dataset?
mean(((hmp_acc_metrics[which(hmp_acc_metrics$category == "NSTI=2"), "recall"] - hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "recall"])/hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "recall"])*100)
mean(((mammal_acc_metrics[which(mammal_acc_metrics$category == "NSTI=2"), "recall"] - mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "recall"])/mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "recall"])*100)
mean(((ocean_acc_metrics[which(ocean_acc_metrics$category == "NSTI=2"), "recall"] - ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "recall"])/ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "recall"])*100)
mean(((blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "recall"] - blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "recall"])/blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "recall"])*100)


sd(((hmp_acc_metrics[which(hmp_acc_metrics$category == "NSTI=2"), "recall"] - hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "recall"])/hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "recall"])*100)
sd(((mammal_acc_metrics[which(mammal_acc_metrics$category == "NSTI=2"), "recall"] - mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "recall"])/mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "recall"])*100)
sd(((ocean_acc_metrics[which(ocean_acc_metrics$category == "NSTI=2"), "recall"] - ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "recall"])/ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "recall"])*100)
sd(((blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "recall"] - blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "recall"])/blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "recall"])*100)


mean(((hmp_acc_metrics[which(hmp_acc_metrics$category == "NSTI=2"), "precision"] - hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "precision"])/hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "precision"])*100)
mean(((mammal_acc_metrics[which(mammal_acc_metrics$category == "NSTI=2"), "precision"] - mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "precision"])/mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "precision"])*100)
mean(((ocean_acc_metrics[which(ocean_acc_metrics$category == "NSTI=2"), "precision"] - ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "precision"])/ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "precision"])*100)
mean(((blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "precision"] - blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "precision"])/blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "precision"])*100)


sd(((hmp_acc_metrics[which(hmp_acc_metrics$category == "NSTI=2"), "precision"] - hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "precision"])/hmp_acc_metrics[which(hmp_acc_metrics$category == "Null"), "precision"])*100)
sd(((mammal_acc_metrics[which(mammal_acc_metrics$category == "NSTI=2"), "precision"] - mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "precision"])/mammal_acc_metrics[which(mammal_acc_metrics$category == "Null"), "precision"])*100)
sd(((ocean_acc_metrics[which(ocean_acc_metrics$category == "NSTI=2"), "precision"] - ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "precision"])/ocean_acc_metrics[which(ocean_acc_metrics$category == "Null"), "precision"])*100)
sd(((blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "precision"] - blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "precision"])/blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "precision"])*100)


mean(((blueberry_acc_metrics[which(blueberry_acc_metrics$category == "NSTI=2"), "precision"] - blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "precision"])/blueberry_acc_metrics[which(blueberry_acc_metrics$category == "Null"), "precision"])*100)
