### Plot 16S vs MGS validations from previously saved RDS files for EC and metacyc pathways.

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

library(ggplot2)
library(reshape2)
library(cowplot)

# Plot KO validations 16S vs MGS.
hmp_ec_spearman <- readRDS("hmp_ec_spearman_df.rds")
mammal_ec_spearman <- readRDS("mammal_ec_spearman_df.rds")
ocean_ec_spearman <- readRDS("ocean_ec_spearman_df.rds")
blueberry_ec_spearman <- readRDS("blueberry_ec_spearman_df.rds")

hmp_pathabun_spearman <- readRDS("hmp_pathabun_spearman_df.rds")
mammal_pathabun_spearman <- readRDS("mammal_pathabun_spearman_df.rds")
ocean_pathabun_spearman <- readRDS("ocean_pathabun_spearman_df.rds")
blueberry_pathabun_spearman <- readRDS("blueberry_pathabun_spearman_df.rds")

hmp_pathcov_spearman <- readRDS("hmp_pathcov_spearman_df.rds")
mammal_pathcov_spearman <- readRDS("mammal_pathcov_spearman_df.rds")
ocean_pathcov_spearman <- readRDS("ocean_pathcov_spearman_df.rds")
blueberry_pathcov_spearman <- readRDS("blueberry_pathcov_spearman_df.rds")

colnames(hmp_ec_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
hmp_ec_spearman <- hmp_ec_spearman[with(hmp_ec_spearman, order(Category, Sample)),]
hmp_ec_spearman$Database <- "PAPRICA"
hmp_ec_spearman$Database[grep("NSTI" , hmp_ec_spearman$Category)] <- "PICRUSt2"
hmp_ec_spearman$Database[which(hmp_ec_spearman$Category=="Null")] <- "Null"
hmp_ec_spearman$Category <- factor(hmp_ec_spearman$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_ec_spearman_melt <- melt(hmp_ec_spearman)

hmp_ec_spearman_boxplots <- ggplot(hmp_ec_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.3, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("HMP")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


colnames(hmp_pathabun_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
hmp_pathabun_spearman <- hmp_pathabun_spearman[with(hmp_pathabun_spearman, order(Category, Sample)),]
hmp_pathabun_spearman$Category <- factor(hmp_pathabun_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_pathabun_spearman$Database <- "PICRUSt2"
hmp_pathabun_spearman$Database[which(hmp_pathabun_spearman$Category=="Null")] <- "Null"
hmp_pathabun_spearman_melt <- melt(hmp_pathabun_spearman)
hmp_pathabun_spearman_boxplots <- ggplot(hmp_pathabun_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("HMP")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))

colnames(hmp_pathcov_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
hmp_pathcov_spearman <- hmp_pathcov_spearman[with(hmp_pathcov_spearman, order(Category, Sample)),]
hmp_pathcov_spearman$Category <- factor(hmp_pathcov_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_pathcov_spearman$Database <- "PICRUSt2"
hmp_pathcov_spearman$Database[which(hmp_pathcov_spearman$Category=="Null")] <- "Null"
hmp_pathcov_spearman_melt <- melt(hmp_pathcov_spearman)
hmp_pathcov_spearman_boxplots <- ggplot(hmp_pathcov_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("HMP")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))



colnames(mammal_ec_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
mammal_ec_spearman <- mammal_ec_spearman[with(mammal_ec_spearman, order(Category, Sample)),]
mammal_ec_spearman$Database <- "PAPRICA"
mammal_ec_spearman$Database[grep("NSTI" , mammal_ec_spearman$Category)] <- "PICRUSt2"
mammal_ec_spearman$Database[which(mammal_ec_spearman$Category=="Null")] <- "Null"
mammal_ec_spearman$Category <- factor(mammal_ec_spearman$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_ec_spearman_melt <- melt(mammal_ec_spearman)

mammal_ec_spearman_boxplots <- ggplot(mammal_ec_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.3, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Mammal")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


colnames(mammal_pathabun_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
mammal_pathabun_spearman <- mammal_pathabun_spearman[with(mammal_pathabun_spearman, order(Category, Sample)),]
mammal_pathabun_spearman$Category <- factor(mammal_pathabun_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_pathabun_spearman$Database <- "PICRUSt2"
mammal_pathabun_spearman$Database[which(mammal_pathabun_spearman$Category=="Null")] <- "Null"
mammal_pathabun_spearman_melt <- melt(mammal_pathabun_spearman)
mammal_pathabun_spearman_boxplots <- ggplot(mammal_pathabun_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Mammal")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))

colnames(mammal_pathcov_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
mammal_pathcov_spearman <- mammal_pathcov_spearman[with(mammal_pathcov_spearman, order(Category, Sample)),]
mammal_pathcov_spearman$Category <- factor(mammal_pathcov_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_pathcov_spearman$Database <- "PICRUSt2"
mammal_pathcov_spearman$Database[which(mammal_pathcov_spearman$Category=="Null")] <- "Null"
mammal_pathcov_spearman_melt <- melt(mammal_pathcov_spearman)
mammal_pathcov_spearman_boxplots <- ggplot(mammal_pathcov_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Mammal")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))


colnames(ocean_ec_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
ocean_ec_spearman <- ocean_ec_spearman[with(ocean_ec_spearman, order(Category, Sample)),]
ocean_ec_spearman$Database <- "PAPRICA"
ocean_ec_spearman$Database[grep("NSTI" , ocean_ec_spearman$Category)] <- "PICRUSt2"
ocean_ec_spearman$Database[which(ocean_ec_spearman$Category=="Null")] <- "Null"
ocean_ec_spearman$Category <- factor(ocean_ec_spearman$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_ec_spearman_melt <- melt(ocean_ec_spearman)

ocean_ec_spearman_boxplots <- ggplot(ocean_ec_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.3, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Ocean")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


colnames(ocean_pathabun_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
ocean_pathabun_spearman <- ocean_pathabun_spearman[with(ocean_pathabun_spearman, order(Category, Sample)),]
ocean_pathabun_spearman$Category <- factor(ocean_pathabun_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_pathabun_spearman$Database <- "PICRUSt2"
ocean_pathabun_spearman$Database[which(ocean_pathabun_spearman$Category=="Null")] <- "Null"
ocean_pathabun_spearman_melt <- melt(ocean_pathabun_spearman)
ocean_pathabun_spearman_boxplots <- ggplot(ocean_pathabun_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Ocean")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))

colnames(ocean_pathcov_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
ocean_pathcov_spearman <- ocean_pathcov_spearman[with(ocean_pathcov_spearman, order(Category, Sample)),]
ocean_pathcov_spearman$Category <- factor(ocean_pathcov_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_pathcov_spearman$Database <- "PICRUSt2"
ocean_pathcov_spearman$Database[which(ocean_pathcov_spearman$Category=="Null")] <- "Null"
ocean_pathcov_spearman_melt <- melt(ocean_pathcov_spearman)
ocean_pathcov_spearman_boxplots <- ggplot(ocean_pathcov_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Ocean")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))


colnames(blueberry_ec_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
blueberry_ec_spearman <- blueberry_ec_spearman[with(blueberry_ec_spearman, order(Category, Sample)),]
blueberry_ec_spearman$Database <- "PAPRICA"
blueberry_ec_spearman$Database[grep("NSTI" , blueberry_ec_spearman$Category)] <- "PICRUSt2"
blueberry_ec_spearman$Database[which(blueberry_ec_spearman$Category=="Null")] <- "Null"
blueberry_ec_spearman$Category <- factor(blueberry_ec_spearman$Category, levels=c("Null", "PAPRICA", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_ec_spearman_melt <- melt(blueberry_ec_spearman)

blueberry_ec_spearman_boxplots <- ggplot(blueberry_ec_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.3, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Soil (Blueberry)")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))


colnames(blueberry_pathabun_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
blueberry_pathabun_spearman <- blueberry_pathabun_spearman[with(blueberry_pathabun_spearman, order(Category, Sample)),]
blueberry_pathabun_spearman$Category <- factor(blueberry_pathabun_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_pathabun_spearman$Database <- "PICRUSt2"
blueberry_pathabun_spearman$Database[which(blueberry_pathabun_spearman$Category=="Null")] <- "Null"
blueberry_pathabun_spearman_melt <- melt(blueberry_pathabun_spearman)
blueberry_pathabun_spearman_boxplots <- ggplot(blueberry_pathabun_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Soil (Blueberry)") + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))

colnames(blueberry_pathcov_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
blueberry_pathcov_spearman <- blueberry_pathcov_spearman[with(blueberry_pathcov_spearman, order(Category, Sample)),]
blueberry_pathcov_spearman$Category <- factor(blueberry_pathcov_spearman$Category, levels=c("Null", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_pathcov_spearman$Database <- "PICRUSt2"
blueberry_pathcov_spearman$Database[which(blueberry_pathcov_spearman$Category=="Null")] <- "Null"
blueberry_pathcov_spearman_melt <- melt(blueberry_pathcov_spearman)
blueberry_pathcov_spearman_boxplots <- ggplot(blueberry_pathcov_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.0, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Soil (Blueberry)") + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#00BFC4"))

plot_grid(hmp_ec_spearman_boxplots, mammal_ec_spearman_boxplots,
          ocean_ec_spearman_boxplots, blueberry_ec_spearman_boxplots,
          labels=c("A", "B", "C", "D"))

plot_grid(hmp_pathabun_spearman_boxplots, mammal_pathabun_spearman_boxplots,
          ocean_pathabun_spearman_boxplots, blueberry_pathabun_spearman_boxplots,
          labels=c("A", "B", "C", "D"))

plot_grid(hmp_pathcov_spearman_boxplots, mammal_pathcov_spearman_boxplots,
          ocean_pathcov_spearman_boxplots, blueberry_pathcov_spearman_boxplots,
          labels=c("A", "B", "C", "D"))


hmp_ec_spearman_nsti2_vs_null <- wilcox.test(hmp_ec_spearman$`Spearman correlation coefficient`[which(hmp_ec_spearman$Category=="NSTI=2")],
                                                    hmp_ec_spearman$`Spearman correlation coefficient`[which(hmp_ec_spearman$Category=="Null")], paired=TRUE)

hmp_ec_spearman_nsti2_vs_paprica <- wilcox.test(hmp_ec_spearman$`Spearman correlation coefficient`[which(hmp_ec_spearman$Category=="NSTI=2")],
                                                    hmp_ec_spearman$`Spearman correlation coefficient`[which(hmp_ec_spearman$Category=="PAPRICA")], paired=TRUE)

hmp_ec_spearman_nsti2_vs_nsti2gg <- wilcox.test(hmp_ec_spearman$`Spearman correlation coefficient`[which(hmp_ec_spearman$Category=="NSTI=2")],
                                                    hmp_ec_spearman$`Spearman correlation coefficient`[which(hmp_ec_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

hmp_ec_spearman_wilcox <- list(nsti2_vs_null=hmp_ec_spearman_nsti2_vs_null,
                            nsti2_vs_paprica=hmp_ec_spearman_nsti2_vs_paprica,
                            nsti2_vs_nsti2gg=hmp_ec_spearman_nsti2_vs_nsti2gg)

hmp_pathabun_spearman_nsti2_vs_null <- wilcox.test(hmp_pathabun_spearman$`Spearman correlation coefficient`[which(hmp_pathabun_spearman$Category=="NSTI=2")],
                                             hmp_pathabun_spearman$`Spearman correlation coefficient`[which(hmp_pathabun_spearman$Category=="Null")], paired=TRUE)

hmp_pathabun_spearman_nsti2_vs_nsti2gg <- wilcox.test(hmp_pathabun_spearman$`Spearman correlation coefficient`[which(hmp_pathabun_spearman$Category=="NSTI=2")],
                                                hmp_pathabun_spearman$`Spearman correlation coefficient`[which(hmp_pathabun_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

hmp_pathabun_spearman_wilcox <- list(nsti2_vs_null=hmp_pathabun_spearman_nsti2_vs_null,
                               nsti2_vs_nsti2gg=hmp_pathabun_spearman_nsti2_vs_nsti2gg)


hmp_pathcov_spearman_nsti2_vs_null <- wilcox.test(hmp_pathcov_spearman$`Spearman correlation coefficient`[which(hmp_pathcov_spearman$Category=="NSTI=2")],
                                             hmp_pathcov_spearman$`Spearman correlation coefficient`[which(hmp_pathcov_spearman$Category=="Null")], paired=TRUE)

hmp_pathcov_spearman_nsti2_vs_nsti2gg <- wilcox.test(hmp_pathcov_spearman$`Spearman correlation coefficient`[which(hmp_pathcov_spearman$Category=="NSTI=2")],
                                                hmp_pathcov_spearman$`Spearman correlation coefficient`[which(hmp_pathcov_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

hmp_pathcov_spearman_wilcox <- list(nsti2_vs_null=hmp_pathcov_spearman_nsti2_vs_null,
                               nsti2_vs_nsti2gg=hmp_pathcov_spearman_nsti2_vs_nsti2gg)

hmp_spearman_wilcox <- list(ec=hmp_ec_spearman_wilcox, pathabun=hmp_pathabun_spearman_wilcox, pathcov=hmp_pathcov_spearman_wilcox)



mammal_ec_spearman_nsti2_vs_null <- wilcox.test(mammal_ec_spearman$`Spearman correlation coefficient`[which(mammal_ec_spearman$Category=="NSTI=2")],
                                             mammal_ec_spearman$`Spearman correlation coefficient`[which(mammal_ec_spearman$Category=="Null")], paired=TRUE)

mammal_ec_spearman_nsti2_vs_paprica <- wilcox.test(mammal_ec_spearman$`Spearman correlation coefficient`[which(mammal_ec_spearman$Category=="NSTI=2")],
                                                mammal_ec_spearman$`Spearman correlation coefficient`[which(mammal_ec_spearman$Category=="PAPRICA")], paired=TRUE)

mammal_ec_spearman_nsti2_vs_nsti2gg <- wilcox.test(mammal_ec_spearman$`Spearman correlation coefficient`[which(mammal_ec_spearman$Category=="NSTI=2")],
                                                mammal_ec_spearman$`Spearman correlation coefficient`[which(mammal_ec_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

mammal_ec_spearman_wilcox <- list(nsti2_vs_null=mammal_ec_spearman_nsti2_vs_null,
                               nsti2_vs_paprica=mammal_ec_spearman_nsti2_vs_paprica,
                               nsti2_vs_nsti2gg=mammal_ec_spearman_nsti2_vs_nsti2gg)

mammal_pathabun_spearman_nsti2_vs_null <- wilcox.test(mammal_pathabun_spearman$`Spearman correlation coefficient`[which(mammal_pathabun_spearman$Category=="NSTI=2")],
                                                   mammal_pathabun_spearman$`Spearman correlation coefficient`[which(mammal_pathabun_spearman$Category=="Null")], paired=TRUE)

mammal_pathabun_spearman_nsti2_vs_nsti2gg <- wilcox.test(mammal_pathabun_spearman$`Spearman correlation coefficient`[which(mammal_pathabun_spearman$Category=="NSTI=2")],
                                                      mammal_pathabun_spearman$`Spearman correlation coefficient`[which(mammal_pathabun_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

mammal_pathabun_spearman_wilcox <- list(nsti2_vs_null=mammal_pathabun_spearman_nsti2_vs_null,
                                     nsti2_vs_nsti2gg=mammal_pathabun_spearman_nsti2_vs_nsti2gg)


mammal_pathcov_spearman_nsti2_vs_null <- wilcox.test(mammal_pathcov_spearman$`Spearman correlation coefficient`[which(mammal_pathcov_spearman$Category=="NSTI=2")],
                                                  mammal_pathcov_spearman$`Spearman correlation coefficient`[which(mammal_pathcov_spearman$Category=="Null")], paired=TRUE)

mammal_pathcov_spearman_nsti2_vs_nsti2gg <- wilcox.test(mammal_pathcov_spearman$`Spearman correlation coefficient`[which(mammal_pathcov_spearman$Category=="NSTI=2")],
                                                     mammal_pathcov_spearman$`Spearman correlation coefficient`[which(mammal_pathcov_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

mammal_pathcov_spearman_wilcox <- list(nsti2_vs_null=mammal_pathcov_spearman_nsti2_vs_null,
                                    nsti2_vs_nsti2gg=mammal_pathcov_spearman_nsti2_vs_nsti2gg)

mammal_spearman_wilcox <- list(ec=mammal_ec_spearman_wilcox, pathabun=mammal_pathabun_spearman_wilcox, pathcov=mammal_pathcov_spearman_wilcox)



ocean_ec_spearman_nsti2_vs_null <- wilcox.test(ocean_ec_spearman$`Spearman correlation coefficient`[which(ocean_ec_spearman$Category=="NSTI=2")],
                                             ocean_ec_spearman$`Spearman correlation coefficient`[which(ocean_ec_spearman$Category=="Null")], paired=TRUE)

ocean_ec_spearman_nsti2_vs_paprica <- wilcox.test(ocean_ec_spearman$`Spearman correlation coefficient`[which(ocean_ec_spearman$Category=="NSTI=2")],
                                                ocean_ec_spearman$`Spearman correlation coefficient`[which(ocean_ec_spearman$Category=="PAPRICA")], paired=TRUE)

ocean_ec_spearman_nsti2_vs_nsti2gg <- wilcox.test(ocean_ec_spearman$`Spearman correlation coefficient`[which(ocean_ec_spearman$Category=="NSTI=2")],
                                                ocean_ec_spearman$`Spearman correlation coefficient`[which(ocean_ec_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

ocean_ec_spearman_wilcox <- list(nsti2_vs_null=ocean_ec_spearman_nsti2_vs_null,
                               nsti2_vs_paprica=ocean_ec_spearman_nsti2_vs_paprica,
                               nsti2_vs_nsti2gg=ocean_ec_spearman_nsti2_vs_nsti2gg)

ocean_pathabun_spearman_nsti2_vs_null <- wilcox.test(ocean_pathabun_spearman$`Spearman correlation coefficient`[which(ocean_pathabun_spearman$Category=="NSTI=2")],
                                                   ocean_pathabun_spearman$`Spearman correlation coefficient`[which(ocean_pathabun_spearman$Category=="Null")], paired=TRUE)

ocean_pathabun_spearman_nsti2_vs_nsti2gg <- wilcox.test(ocean_pathabun_spearman$`Spearman correlation coefficient`[which(ocean_pathabun_spearman$Category=="NSTI=2")],
                                                      ocean_pathabun_spearman$`Spearman correlation coefficient`[which(ocean_pathabun_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

ocean_pathabun_spearman_wilcox <- list(nsti2_vs_null=ocean_pathabun_spearman_nsti2_vs_null,
                                     nsti2_vs_nsti2gg=ocean_pathabun_spearman_nsti2_vs_nsti2gg)


ocean_pathcov_spearman_nsti2_vs_null <- wilcox.test(ocean_pathcov_spearman$`Spearman correlation coefficient`[which(ocean_pathcov_spearman$Category=="NSTI=2")],
                                                  ocean_pathcov_spearman$`Spearman correlation coefficient`[which(ocean_pathcov_spearman$Category=="Null")], paired=TRUE)

ocean_pathcov_spearman_nsti2_vs_nsti2gg <- wilcox.test(ocean_pathcov_spearman$`Spearman correlation coefficient`[which(ocean_pathcov_spearman$Category=="NSTI=2")],
                                                     ocean_pathcov_spearman$`Spearman correlation coefficient`[which(ocean_pathcov_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

ocean_pathcov_spearman_wilcox <- list(nsti2_vs_null=ocean_pathcov_spearman_nsti2_vs_null,
                                    nsti2_vs_nsti2gg=ocean_pathcov_spearman_nsti2_vs_nsti2gg)

ocean_spearman_wilcox <- list(ec=ocean_ec_spearman_wilcox, pathabun=ocean_pathabun_spearman_wilcox, pathcov=ocean_pathcov_spearman_wilcox)




blueberry_ec_spearman_nsti2_vs_null <- wilcox.test(blueberry_ec_spearman$`Spearman correlation coefficient`[which(blueberry_ec_spearman$Category=="NSTI=2")],
                                             blueberry_ec_spearman$`Spearman correlation coefficient`[which(blueberry_ec_spearman$Category=="Null")], paired=TRUE)

blueberry_ec_spearman_nsti2_vs_paprica <- wilcox.test(blueberry_ec_spearman$`Spearman correlation coefficient`[which(blueberry_ec_spearman$Category=="NSTI=2")],
                                                blueberry_ec_spearman$`Spearman correlation coefficient`[which(blueberry_ec_spearman$Category=="PAPRICA")], paired=TRUE)

blueberry_ec_spearman_nsti2_vs_nsti2gg <- wilcox.test(blueberry_ec_spearman$`Spearman correlation coefficient`[which(blueberry_ec_spearman$Category=="NSTI=2")],
                                                blueberry_ec_spearman$`Spearman correlation coefficient`[which(blueberry_ec_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_ec_spearman_wilcox <- list(nsti2_vs_null=blueberry_ec_spearman_nsti2_vs_null,
                               nsti2_vs_paprica=blueberry_ec_spearman_nsti2_vs_paprica,
                               nsti2_vs_nsti2gg=blueberry_ec_spearman_nsti2_vs_nsti2gg)

blueberry_pathabun_spearman_nsti2_vs_null <- wilcox.test(blueberry_pathabun_spearman$`Spearman correlation coefficient`[which(blueberry_pathabun_spearman$Category=="NSTI=2")],
                                                   blueberry_pathabun_spearman$`Spearman correlation coefficient`[which(blueberry_pathabun_spearman$Category=="Null")], paired=TRUE)

blueberry_pathabun_spearman_nsti2_vs_nsti2gg <- wilcox.test(blueberry_pathabun_spearman$`Spearman correlation coefficient`[which(blueberry_pathabun_spearman$Category=="NSTI=2")],
                                                      blueberry_pathabun_spearman$`Spearman correlation coefficient`[which(blueberry_pathabun_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_pathabun_spearman_wilcox <- list(nsti2_vs_null=blueberry_pathabun_spearman_nsti2_vs_null,
                                     nsti2_vs_nsti2gg=blueberry_pathabun_spearman_nsti2_vs_nsti2gg)


blueberry_pathcov_spearman_nsti2_vs_null <- wilcox.test(blueberry_pathcov_spearman$`Spearman correlation coefficient`[which(blueberry_pathcov_spearman$Category=="NSTI=2")],
                                                  blueberry_pathcov_spearman$`Spearman correlation coefficient`[which(blueberry_pathcov_spearman$Category=="Null")], paired=TRUE)

blueberry_pathcov_spearman_nsti2_vs_nsti2gg <- wilcox.test(blueberry_pathcov_spearman$`Spearman correlation coefficient`[which(blueberry_pathcov_spearman$Category=="NSTI=2")],
                                                     blueberry_pathcov_spearman$`Spearman correlation coefficient`[which(blueberry_pathcov_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_pathcov_spearman_wilcox <- list(nsti2_vs_null=blueberry_pathcov_spearman_nsti2_vs_null,
                                    nsti2_vs_nsti2gg=blueberry_pathcov_spearman_nsti2_vs_nsti2gg)

blueberry_spearman_wilcox <- list(ec=blueberry_ec_spearman_wilcox, pathabun=blueberry_pathabun_spearman_wilcox, pathcov=blueberry_pathcov_spearman_wilcox)


combined_spearman_wilcox_ec_path <- list(hmp=hmp_spearman_wilcox,
                                 mammal=mammal_spearman_wilcox,
                                 ocean=ocean_spearman_wilcox,
                                 blueberry=blueberry_spearman_wilcox)

saveRDS(object = combined_spearman_wilcox_ec_path, file="combined_spearman_wilcox_ec_path.rds")

### Read in when needed: ###
combined_spearman_wilcox_ec_path <- readRDS("combined_spearman_wilcox_ec_path.rds")

sapply(combined_spearman_wilcox_ec_path$hmp$ec, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$mammal$ec, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$ocean$ec, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$blueberry$ec, function(x){x$p.value})


sapply(combined_spearman_wilcox_ec_path$hmp$pathabun, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$mammal$pathabun, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$ocean$pathabun, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$blueberry$pathabun, function(x){x$p.value})



sapply(combined_spearman_wilcox_ec_path$hmp$pathcov, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$mammal$pathcov, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$ocean$pathcov, function(x){x$p.value})
sapply(combined_spearman_wilcox_ec_path$blueberry$pathcov, function(x){x$p.value})


wilcox.test(hmp_ec_spearman[which(hmp_ec_spearman$cat == "PAPRICA"), "metric"],
            hmp_ec_spearman[which(hmp_ec_spearman$cat == "Null"), "metric"])

wilcox.test(mammal_ec_spearman[which(mammal_ec_spearman$cat == "PAPRICA"), "metric"],
            mammal_ec_spearman[which(mammal_ec_spearman$cat == "Null"), "metric"])

wilcox.test(ocean_ec_spearman[which(ocean_ec_spearman$cat == "PAPRICA"), "metric"],
            ocean_ec_spearman[which(ocean_ec_spearman$cat == "Null"), "metric"])

wilcox.test(blueberry_ec_spearman[which(blueberry_ec_spearman$cat == "PAPRICA"), "metric"],
            blueberry_ec_spearman[which(blueberry_ec_spearman$cat == "Null"), "metric"])