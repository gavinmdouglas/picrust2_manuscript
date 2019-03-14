
rm(list=ls())

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data")

wilcoxon_df <- data.frame(
  group1="Assembly Null",
  group2="Assembly",
  y.position = 1.05,
  p_symbol= "**")

LOOCV_out_18S <- readRDS("saved_RDS/18S_LOOCV_metrics.rds")
new_levels <- rev(levels(LOOCV_out_18S$level))
LOOCV_out_18S$level <- as.character(LOOCV_out_18S$level)
LOOCV_out_18S$level <- factor(LOOCV_out_18S$level, levels=new_levels)
LOOCV_out_18S_plot <- ggplot(LOOCV_out_18S, aes(x=level, y=mean_rho, fill=c("coral3"))) +
  geom_boxplot() +
  ggtitle("Fungi 18S EC Numbers") +
  guides(fill=FALSE) +
  xlab("") +
  ylim(c(0, 1.1)) +
  ylab(c("Mean Spearman's Correlation Coefficient")) +
  scale_fill_manual(values=c("coral3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(wilcoxon_df, label = "p_symbol")


LOOCV_out_ITS <- readRDS("saved_RDS/ITS_LOOCV_metrics.rds")
new_levels <- rev(levels(LOOCV_out_ITS$level))
LOOCV_out_ITS$level <- as.character(LOOCV_out_ITS$level)
LOOCV_out_ITS$level <- factor(LOOCV_out_ITS$level, levels=new_levels)
LOOCV_out_ITS_plot <- ggplot(LOOCV_out_ITS, aes(x=level, y=mean_rho, fill=c("dodgerblue3"))) +
  geom_boxplot() +
  xlab("") +
  ggtitle("Fungi ITS EC Numbers") +
  guides(fill=FALSE) +
  ylim(c(0, 1.05)) +
  ylab(c("Mean Spearman's Correlation Coefficient")) +
  scale_fill_manual(values=c("dodgerblue3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(wilcoxon_df, label = "p_symbol")


LOOCV_out_16S <- readRDS("saved_RDS/16S_LOOCV_metrics.rds")
new_levels <- rev(levels(LOOCV_out_16S$level))
LOOCV_out_16S$level <- as.character(LOOCV_out_16S$level)
LOOCV_out_16S$level <- factor(LOOCV_out_16S$level, levels=new_levels)
LOOCV_out_16S_plot <- ggplot(LOOCV_out_16S, aes(x=level, y=mean_rho, fill=c("dark orange"))) +
  geom_boxplot() +
  xlab("") +
  ggtitle("Prokaryotic 16S KOs") +
  guides(fill=FALSE) +
  ylim(c(0, 1.05)) +
  ylab(c("Mean Spearman's Correlation Coefficient")) +
  scale_fill_manual(values=c("dark orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(wilcoxon_df, label = "p_symbol")



plot_grid(LOOCV_out_18S_plot,
          LOOCV_out_ITS_plot,
          LOOCV_out_16S_plot,
          labels=c("A", "B", "C"),
          ncol=3,
          nrow=1)

wilcox.test(LOOCV_out_18S$mean_rho[which(LOOCV_out_18S$level=="Assembly")], LOOCV_out_18S$mean_rho[which(LOOCV_out_18S$level=="Assembly Null")])$p.value
wilcox.test(LOOCV_out_ITS$mean_rho[which(LOOCV_out_ITS$level=="Assembly")], LOOCV_out_ITS$mean_rho[which(LOOCV_out_ITS$level=="Assembly Null")])$p.value
wilcox.test(LOOCV_out_16S$mean_rho[which(LOOCV_out_16S$level=="Assembly")], LOOCV_out_16S$mean_rho[which(LOOCV_out_16S$level=="Assembly Null")])$p.value

