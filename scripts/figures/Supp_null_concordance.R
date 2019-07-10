rm(list=ls(all=TRUE))

library(ggplot2)
library(ggbeeswarm)

combined_reps_melt <- readRDS(file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/null_concordance.rds")

combined_reps_melt <- combined_reps_melt[-which(combined_reps_melt$num_genomes %in% c("2", "3", "4", "6", "7", "8", "9", "10", "75")), ]

ggplot(combined_reps_melt, aes(x=num_genomes, y=value, fill=db)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  scale_fill_manual(values=c("#a6cee3" ,  "#1f78b4" , "#b2df8a" , "#33a02c" , "#fb9a99")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Spearman correlation coefficient") +
  xlab("Number of genomes") +
  #scale_x_discrete(labels=gsub("n.", "", levels(combined_reps_melt$variable))) +
  guides(fill=guide_legend(title="Database")) +
  facet_grid(. ~ num_genomes, scales = "free", space = "free", switch="x") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())