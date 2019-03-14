### Plot 16S vs MGS validations from previously saved RDS files.

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

library(ggplot2)
library(reshape2)
library(cowplot)

# Plot KO validations 16S vs MGS.
hmp_ko_spearman <- readRDS("hmp_ko_spearman_df.rds")
mammal_ko_spearman <- readRDS("mammal_ko_spearman_df.rds")
ocean_ko_spearman <- readRDS("ocean_ko_spearman_df.rds")
blueberry_ko_spearman <- readRDS("blueberry_ko_spearman_df.rds")


colnames(hmp_ko_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
hmp_ko_spearman <- hmp_ko_spearman[with(hmp_ko_spearman, order(Category, Sample)),]
hmp_ko_spearman$Database <- "Other"
hmp_ko_spearman$Database[grep("NSTI" , hmp_ko_spearman$Category)] <- "PICRUSt2"
hmp_ko_spearman$Database[which(hmp_ko_spearman$Category=="Null")] <- "Null"
hmp_ko_spearman$Category <- factor(hmp_ko_spearman$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
hmp_ko_spearman_melt <- melt(hmp_ko_spearman)

hmp_spearman_boxplots <- ggplot(hmp_ko_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("HMP")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

colnames(mammal_ko_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
mammal_ko_spearman <- mammal_ko_spearman[with(mammal_ko_spearman, order(Category, Sample)),]
mammal_ko_spearman$Database <- "Other"
mammal_ko_spearman$Database[grep("NSTI" , mammal_ko_spearman$Category)] <- "PICRUSt2"
mammal_ko_spearman$Database[which(mammal_ko_spearman$Category=="Null")] <- "Null"
mammal_ko_spearman$Category <- factor(mammal_ko_spearman$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
mammal_ko_spearman_melt <- melt(mammal_ko_spearman)

mammal_spearman_boxplots <- ggplot(mammal_ko_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Mammal")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))



colnames(ocean_ko_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
ocean_ko_spearman <- ocean_ko_spearman[with(ocean_ko_spearman, order(Category, Sample)),]
ocean_ko_spearman$Database <- "Other"
ocean_ko_spearman$Database[grep("NSTI" , ocean_ko_spearman$Category)] <- "PICRUSt2"
ocean_ko_spearman$Database[which(ocean_ko_spearman$Category=="Null")] <- "Null"
ocean_ko_spearman$Category <- factor(ocean_ko_spearman$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
ocean_ko_spearman_melt <- melt(ocean_ko_spearman)

ocean_spearman_boxplots <- ggplot(ocean_ko_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Ocean")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

colnames(blueberry_ko_spearman) <- c("Spearman correlation coefficient", "Category", "metric_type", "Sample")
blueberry_ko_spearman <- blueberry_ko_spearman[with(blueberry_ko_spearman, order(Category, Sample)),]
blueberry_ko_spearman$Database <- "Other"
blueberry_ko_spearman$Database[grep("NSTI" , blueberry_ko_spearman$Category)] <- "PICRUSt2"
blueberry_ko_spearman$Database[which(blueberry_ko_spearman$Category=="Null")] <- "Null"
blueberry_ko_spearman$Category <- factor(blueberry_ko_spearman$Category, levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "NSTI=2 (GG)", "NSTI=2", "NSTI=1.5", "NSTI=1", "NSTI=0.5", "NSTI=0.25", "NSTI=0.1", "NSTI=0.05"))
blueberry_ko_spearman_melt <- melt(blueberry_ko_spearman)

blueberry_spearman_boxplots <- ggplot(blueberry_ko_spearman_melt, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  coord_flip() + ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Soil (Blueberry)")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4"))

plot_grid(hmp_spearman_boxplots,
          mammal_spearman_boxplots,
          ocean_spearman_boxplots,
          blueberry_spearman_boxplots,
          labels=c("A", "B", "C", "D"))

# Also plot without NSTI cut-offs specified.
hmp_ko_spearman_melt_simplified <- hmp_ko_spearman_melt
hmp_ko_spearman_melt_simplified$Category <- as.character(hmp_ko_spearman_melt_simplified$Category)
hmp_ko_spearman_melt_simplified$Category[which(hmp_ko_spearman_melt_simplified$Category=="NSTI=2")] <- "PICRUSt2"
hmp_ko_spearman_melt_simplified$Category[which(hmp_ko_spearman_melt_simplified$Category=="NSTI=2 (GG)")] <- "PICRUSt2 (GG)"
hmp_ko_spearman_melt_simplified <- hmp_ko_spearman_melt_simplified[-grep("NSTI", hmp_ko_spearman_melt_simplified$Category), ]
hmp_ko_spearman_melt_simplified$Category <- factor(hmp_ko_spearman_melt_simplified$Category,
                                                   levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2 (GG)", "PICRUSt2"))
hmp_simplified_plot <- ggplot(hmp_ko_spearman_melt_simplified, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
                          ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("HMP")  + guides(fill=FALSE) +
                           scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("")

mammal_ko_spearman_melt_simplified <- mammal_ko_spearman_melt
mammal_ko_spearman_melt_simplified$Category <- as.character(mammal_ko_spearman_melt_simplified$Category)
mammal_ko_spearman_melt_simplified$Category[which(mammal_ko_spearman_melt_simplified$Category=="NSTI=2")] <- "PICRUSt2"
mammal_ko_spearman_melt_simplified$Category[which(mammal_ko_spearman_melt_simplified$Category=="NSTI=2 (GG)")] <- "PICRUSt2 (GG)"
mammal_ko_spearman_melt_simplified <- mammal_ko_spearman_melt_simplified[-grep("NSTI", mammal_ko_spearman_melt_simplified$Category), ]
mammal_ko_spearman_melt_simplified$Category <- factor(mammal_ko_spearman_melt_simplified$Category,
                                                   levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2 (GG)", "PICRUSt2"))
mammal_simplified_plot <- ggplot(mammal_ko_spearman_melt_simplified, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Mammal")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("")

ocean_ko_spearman_melt_simplified <- ocean_ko_spearman_melt
ocean_ko_spearman_melt_simplified$Category <- as.character(ocean_ko_spearman_melt_simplified$Category)
ocean_ko_spearman_melt_simplified$Category[which(ocean_ko_spearman_melt_simplified$Category=="NSTI=2")] <- "PICRUSt2"
ocean_ko_spearman_melt_simplified$Category[which(ocean_ko_spearman_melt_simplified$Category=="NSTI=2 (GG)")] <- "PICRUSt2 (GG)"
ocean_ko_spearman_melt_simplified <- ocean_ko_spearman_melt_simplified[-grep("NSTI", ocean_ko_spearman_melt_simplified$Category), ]
ocean_ko_spearman_melt_simplified$Category <- factor(ocean_ko_spearman_melt_simplified$Category,
                                                   levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2 (GG)", "PICRUSt2"))
ocean_simplified_plot <- ggplot(ocean_ko_spearman_melt_simplified, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Ocean")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("")

blueberry_ko_spearman_melt_simplified <- blueberry_ko_spearman_melt
blueberry_ko_spearman_melt_simplified$Category <- as.character(blueberry_ko_spearman_melt_simplified$Category)
blueberry_ko_spearman_melt_simplified$Category[which(blueberry_ko_spearman_melt_simplified$Category=="NSTI=2")] <- "PICRUSt2"
blueberry_ko_spearman_melt_simplified$Category[which(blueberry_ko_spearman_melt_simplified$Category=="NSTI=2 (GG)")] <- "PICRUSt2 (GG)"
blueberry_ko_spearman_melt_simplified <- blueberry_ko_spearman_melt_simplified[-grep("NSTI", blueberry_ko_spearman_melt_simplified$Category), ]
blueberry_ko_spearman_melt_simplified$Category <- factor(blueberry_ko_spearman_melt_simplified$Category,
                                                   levels=c("Null", "Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2 (GG)", "PICRUSt2"))
blueberry_simplified_plot <- ggplot(blueberry_ko_spearman_melt_simplified, aes(x=Category, y=value, fill=Database)) + geom_boxplot() +
  ylim(c(0.5, 1.0)) + ylab(c("Spearman correlation coefficient")) + ggtitle("Soil (Blueberry)")  + guides(fill=FALSE) +
  scale_fill_manual(values=c("light grey", "#F8766D", "#00BFC4")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + xlab("")
  

plot_grid(hmp_simplified_plot,
          mammal_simplified_plot,
          ocean_simplified_plot,
          blueberry_simplified_plot,
          labels=c("A", "B", "C", "D"))

hmp_spearman_nsti2_vs_null <- wilcox.test(hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2")],
                                                    hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="Null")], paired=TRUE)

hmp_spearman_nsti2_vs_tax4fun <- wilcox.test(hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2")],
                                                    hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="Tax4Fun")], paired=TRUE)

hmp_spearman_nsti2_vs_panfp <- wilcox.test(hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2")],
                                                    hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="PanFP")], paired=TRUE)

hmp_spearman_nsti2_vs_piphillin <- wilcox.test(hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2")],
                                                    hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="Piphillin")], paired=TRUE)


hmp_spearman_nsti2_vs_picrust1 <- wilcox.test(hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2")],
                                                    hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="PICRUSt1")], paired=TRUE)


hmp_spearman_nsti2_vs_nsti2gg <- wilcox.test(hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2")],
                                                    hmp_ko_spearman$`Spearman correlation coefficient`[which(hmp_ko_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

hmp_spearman_wilcox <- list(nsti2_vs_null=hmp_spearman_nsti2_vs_null,
                            nsti2_vs_tax4fun=hmp_spearman_nsti2_vs_tax4fun,
                            nsti2_vs_panfp=hmp_spearman_nsti2_vs_panfp,
                            nsti2_vs_piphillin=hmp_spearman_nsti2_vs_piphillin,
                            nsti2_vs_picrust1=hmp_spearman_nsti2_vs_picrust1,
                            nsti2_vs_nsti2gg=hmp_spearman_nsti2_vs_nsti2gg)

mammal_spearman_nsti2_vs_null <- wilcox.test(mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2")],
                                          mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="Null")], paired=TRUE)

mammal_spearman_nsti2_vs_tax4fun <- wilcox.test(mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2")],
                                             mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="Tax4Fun")], paired=TRUE)

mammal_spearman_nsti2_vs_panfp <- wilcox.test(mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2")],
                                           mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="PanFP")], paired=TRUE)

mammal_spearman_nsti2_vs_piphillin <- wilcox.test(mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2")],
                                               mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="Piphillin")], paired=TRUE)


mammal_spearman_nsti2_vs_picrust1 <- wilcox.test(mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2")],
                                              mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="PICRUSt1")], paired=TRUE)


mammal_spearman_nsti2_vs_nsti2gg <- wilcox.test(mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2")],
                                             mammal_ko_spearman$`Spearman correlation coefficient`[which(mammal_ko_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

mammal_spearman_wilcox <- list(nsti2_vs_null=mammal_spearman_nsti2_vs_null,
                            nsti2_vs_tax4fun=mammal_spearman_nsti2_vs_tax4fun,
                            nsti2_vs_panfp=mammal_spearman_nsti2_vs_panfp,
                            nsti2_vs_piphillin=mammal_spearman_nsti2_vs_piphillin,
                            nsti2_vs_picrust1=mammal_spearman_nsti2_vs_picrust1,
                            nsti2_vs_nsti2gg=mammal_spearman_nsti2_vs_nsti2gg)


ocean_spearman_nsti2_vs_null <- wilcox.test(ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2")],
                                          ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="Null")], paired=TRUE)

ocean_spearman_nsti2_vs_tax4fun <- wilcox.test(ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2")],
                                             ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="Tax4Fun")], paired=TRUE)

ocean_spearman_nsti2_vs_panfp <- wilcox.test(ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2")],
                                           ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="PanFP")], paired=TRUE)

ocean_spearman_nsti2_vs_piphillin <- wilcox.test(ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2")],
                                               ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="Piphillin")], paired=TRUE)


ocean_spearman_nsti2_vs_picrust1 <- wilcox.test(ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2")],
                                              ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="PICRUSt1")], paired=TRUE)


ocean_spearman_nsti2_vs_nsti2gg <- wilcox.test(ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2")],
                                             ocean_ko_spearman$`Spearman correlation coefficient`[which(ocean_ko_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

ocean_spearman_wilcox <- list(nsti2_vs_null=ocean_spearman_nsti2_vs_null,
                            nsti2_vs_tax4fun=ocean_spearman_nsti2_vs_tax4fun,
                            nsti2_vs_panfp=ocean_spearman_nsti2_vs_panfp,
                            nsti2_vs_piphillin=ocean_spearman_nsti2_vs_piphillin,
                            nsti2_vs_picrust1=ocean_spearman_nsti2_vs_picrust1,
                            nsti2_vs_nsti2gg=ocean_spearman_nsti2_vs_nsti2gg)


blueberry_spearman_nsti2_vs_null <- wilcox.test(blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2")],
                                          blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="Null")], paired=TRUE)

blueberry_spearman_nsti2_vs_tax4fun <- wilcox.test(blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2")],
                                             blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="Tax4Fun")], paired=TRUE)

blueberry_spearman_nsti2_vs_panfp <- wilcox.test(blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2")],
                                           blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="PanFP")], paired=TRUE)

blueberry_spearman_nsti2_vs_piphillin <- wilcox.test(blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2")],
                                               blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="Piphillin")], paired=TRUE)


blueberry_spearman_nsti2_vs_picrust1 <- wilcox.test(blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2")],
                                              blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="PICRUSt1")], paired=TRUE)


blueberry_spearman_nsti2_vs_nsti2gg <- wilcox.test(blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2")],
                                             blueberry_ko_spearman$`Spearman correlation coefficient`[which(blueberry_ko_spearman$Category=="NSTI=2 (GG)")], paired=TRUE)

blueberry_spearman_wilcox <- list(nsti2_vs_null=blueberry_spearman_nsti2_vs_null,
                            nsti2_vs_tax4fun=blueberry_spearman_nsti2_vs_tax4fun,
                            nsti2_vs_panfp=blueberry_spearman_nsti2_vs_panfp,
                            nsti2_vs_piphillin=blueberry_spearman_nsti2_vs_piphillin,
                            nsti2_vs_picrust1=blueberry_spearman_nsti2_vs_picrust1,
                            nsti2_vs_nsti2gg=blueberry_spearman_nsti2_vs_nsti2gg)

combined_spearman_wilcox <- list(hmp=hmp_spearman_wilcox,
                                 mammal=mammal_spearman_wilcox,
                                 ocean=ocean_spearman_wilcox,
                                 blueberry=blueberry_spearman_wilcox)

saveRDS(object = combined_spearman_wilcox, file="combined_spearman_wilcox_ko.rds")

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")
combined_spearman_wilcox <- readRDS("combined_spearman_wilcox_ko.rds")

sapply(combined_spearman_wilcox$hmp, function(x) { x$p.value })
sapply(combined_spearman_wilcox$mammal, function(x) { x$p.value })
sapply(combined_spearman_wilcox$ocean, function(x) { x$p.value })
sapply(combined_spearman_wilcox$blueberry, function(x) { x$p.value })


# Mean and sd values:
mean(hmp_ko_spearman[hmp_ko_spearman$cat=="NSTI=2", "metric"])
sd(hmp_ko_spearman[hmp_ko_spearman$cat=="NSTI=2", "metric"])

mean(mammal_ko_spearman[mammal_ko_spearman$cat=="NSTI=2", "metric"])
sd(mammal_ko_spearman[mammal_ko_spearman$cat=="NSTI=2", "metric"])

mean(ocean_ko_spearman[ocean_ko_spearman$cat=="NSTI=2", "metric"])
sd(ocean_ko_spearman[ocean_ko_spearman$cat=="NSTI=2", "metric"])

mean(blueberry_ko_spearman[blueberry_ko_spearman$cat=="NSTI=2", "metric"])
sd(blueberry_ko_spearman[blueberry_ko_spearman$cat=="NSTI=2", "metric"])
