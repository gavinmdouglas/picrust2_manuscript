### Plot number of genera contributing to each pathway in each 16S-MGS paired dataset.

rm(list=ls(all.names=TRUE))

library(ggplot2)
library(reshape2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

all_16S_vs_mgs_contrib <- readRDS("16S_vs_MGS_strat_contrib.rds")

hmp_num_contrib_genera_mgs_vs_16S <- ggplot(all_16S_vs_mgs_contrib$hmp$num_contrib, aes(x=mgs_mean, y=picrust2_mean)) +
  geom_point(size=2) +
  scale_x_continuous(limits=c(0, 30), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, 20), expand = c(0, 0)) +
  ggtitle("Mean # contributing genera (HMP)") +
  ylab("PICRUSt2") +
  xlab("Metagenomics") +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

mammal_num_contrib_genera_mgs_vs_16S <- ggplot(all_16S_vs_mgs_contrib$mammal$num_contrib, aes(x=mgs_mean, y=picrust2_mean)) +
  geom_point(size=2) +
  scale_x_continuous(limits=c(0, 15), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, 15), expand = c(0, 0)) +
  ggtitle("Mean # contributing genera (Mammalian stool)") +
  ylab("PICRUSt2") +
  xlab("Metagenomics") +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ocean_num_contrib_genera_mgs_vs_16S <- ggplot(all_16S_vs_mgs_contrib$ocean$num_contrib, aes(x=mgs_mean, y=picrust2_mean)) +
  geom_point(size=2) +
  scale_x_continuous(limits=c(0, 15), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, 55), expand = c(0, 0)) +
  ggtitle("Mean # contributing genera (Ocean)") +
  ylab("PICRUSt2") +
  xlab("Metagenomics") +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

blueberry_num_contrib_genera_mgs_vs_16S <- ggplot(all_16S_vs_mgs_contrib$blueberry$num_contrib, aes(x=mgs_mean, y=picrust2_mean)) +
  geom_point(size=2) +
  scale_x_continuous(limits=c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, 40), expand = c(0, 0)) +
  ggtitle("Mean # contributing genera (Blueberry soil)") +
  ylab("PICRUSt2") +
  xlab("Metagenomics") +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


pdf(file = "../../figures/Supp_num_contrib_genera.pdf", width=9, height=9)

plot_grid(hmp_num_contrib_genera_mgs_vs_16S,
          mammal_num_contrib_genera_mgs_vs_16S,
          ocean_num_contrib_genera_mgs_vs_16S,
          blueberry_num_contrib_genera_mgs_vs_16S,
          labels=c("a", "b", "c", "d"),
          nrow=2,
          ncol=2)

dev.off()