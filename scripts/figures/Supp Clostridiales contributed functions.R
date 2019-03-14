### Plot functions contributed by Clostridiales.

rm(list=ls(all=TRUE))

library(ggplot2)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

nonibd_sig_higher_ratio_prep_melt <- readRDS("results_out/nonibd_sig_higher_ratio_prep_melt.rds")

nonibd_sig_higher_ratio_prep_melt$variable <- as.character(nonibd_sig_higher_ratio_prep_melt$variable)

ggplot(nonibd_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=diagnosis)) +
  geom_boxplot(width=0.75, outlier.shape = NA) +
  coord_flip() +
  scale_fill_manual(values=c("black", "grey")) +
  xlab("") +
  ylab("log2((Contributed by Significant Clostridiales + 1)/(Contributed by Other + 1))") +
  labs(fill="Disease State") +
  theme(legend.position = c(0.8, 0.9),
        legend.background = element_rect(color = "black", 
                                         fill = "white", size = 0.2, linetype = "solid"))
