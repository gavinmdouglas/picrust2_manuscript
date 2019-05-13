### Plot example results from HMP2 dataset.

rm(list=ls(all=TRUE))

library(ggplot2)
library(reshape2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

# Panel A.
# First plot pathways enriched in Proteobacteria (CD vs nonIBD).
cd_sig_higher_ratio_prep_melt <- readRDS("results_out/cd_sig_higher_ratio_prep_melt.rds")

cd_sig_higher_ratio_prep_melt$variable <- as.character(cd_sig_higher_ratio_prep_melt$variable)
cd_sig_higher_ratio_prep_melt$variable <- factor(cd_sig_higher_ratio_prep_melt$variable, levels=c("PWY0-42", "PWY-5189", "PWY-5188"))

cd_sig_higher_ratio_plot <- ggplot(cd_sig_higher_ratio_prep_melt, aes(x=variable, y=log2ratio, fill=diagnosis)) +
  geom_boxplot(width=0.75, outlier.shape = NA) +
  coord_flip() +
  scale_fill_manual(values=c("black", "grey")) +
  xlab("") +
  ylab(expression('log'[2]*'((Contributed by Proteobacteria + 1)/(Contributed by Other + 1))')) +
  scale_y_continuous(limits=c(-5, 8)) +
  labs(fill="Phenotype") +
  theme(legend.position = c(0.7, 0.8),
        legend.background = element_rect(color = "black", 
                                         fill = "white", size = 0.2, linetype = "solid"))

# Panel B
# Plot scatterplot of # contributors in MGS stoopl vs 16S ileal samples.
num_contrib_genera_out <- readRDS("results_out/num_contrib_genera_out.rds")
num_contrib_genera_mgs_vs_16S <- ggplot(num_contrib_genera_out, aes(x=mgs_mean, y=picrust2_mean)) +
  geom_point(size=2) +
  scale_x_continuous(limits=c(0, 25), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("# Genera Contributing to each Pathway") +
  ylab("PICRUSt2 (Ileum)") +
  xlab("Metagenomics (Stool)") +
  coord_cartesian(clip = 'off')

# Panel C
# Plot scatterplot of top 10 genera contributors to tetrapyrolle I in PICRUSt2 / MGS.
PWY_5188_breakdown <- readRDS("results_out/PWY_5188_breakdown.rds")

# Identify top contributors in either PICRUSt2 or MGS and set everything else to be "Other".
PWY_5188_breakdown$Genus <- PWY_5188_breakdown$Row.names

picrust2_cutoff_abun <- sort(PWY_5188_breakdown$picrust2_mean, decreasing = TRUE)[10]
mgs_cutoff_abun <- sort(PWY_5188_breakdown$mgs_mean, decreasing = TRUE)[10]
rows2keep <- which(PWY_5188_breakdown$picrust2_mean >= picrust2_cutoff_abun | PWY_5188_breakdown$mgs_mean >= mgs_cutoff_abun)
PWY_5188_breakdown$Genus[-rows2keep] <- "Other"

genera_ordered <- c("Blautia", "Clostridium", "Coprococcus", "Dialister", "Roseburia", "Ruminococcus", 
                    "Veillonella", "Fusobacterium", "Bilophila", "Campylobacter", "Escherichia", 
                    "Haemophilus", "Sutterella", "Akkermansia", "Other")

# Colours chosen based on phyla of all genera.
manual_col <- c("firebrick", "firebrick2", "deeppink", "deeppink3", "indianred1", "indianred3", "red", 
                "yellow3", "steelblue", "steelblue2", "royalblue2", "lightslateblue", "dodgerblue", 
                "black", "grey")

PWY_5188_breakdown$Genus <- gsub("g__", "", PWY_5188_breakdown$Genus)

PWY_5188_breakdown$Genus <- factor(PWY_5188_breakdown$Genus, levels=(genera_ordered))

PWY_5188_breakdown_scatterplot <- ggplot(data = PWY_5188_breakdown, aes(x = mgs_mean, y = picrust2_mean, colour=Genus)) + 
  geom_point(size=4) +
  scale_colour_manual(values=manual_col, labels = c(substitute(italic(genus), env=list(genus=genera_ordered[1])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[2])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[3])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[4])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[5])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[6])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[7])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[8])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[9])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[10])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[11])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[12])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[13])),
                                                    substitute(italic(genus), env=list(genus=genera_ordered[14])),
                                                    "Other")) +
  scale_x_continuous(limits=c(0, 0.07), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, 0.04), expand = c(0, 0)) +
  coord_cartesian(clip = 'off') +
  ggtitle("PWY-5188: tetrapyrrole biosynthesis I\n(from glutamate)") +
  xlab("Metagenomics (Stool) Contributing %") +
  ylab("PICRUSt2 (Ileum) Contributing %") +
  theme(legend.text.align = 0)

# Plot stacked barcharts of main contributing to PWY-6572 and PWY0-1533

hmp2_16S_pathabun_strat_genus_sum <- readRDS("results_out/hmp2_16S_pathabun_strat_full_genus_sum.rds")


# Remove "Bacteria" from genus string:
hmp2_16S_pathabun_strat_genus_sum$genus <- gsub("k__Bacteria; ", "", hmp2_16S_pathabun_strat_genus_sum$genus)

# Identify genera contributing most abundance across both pathways.
hmp2_16S_pathabun_strat_genus_sum_PWY_6572 <- hmp2_16S_pathabun_strat_genus_sum[which(hmp2_16S_pathabun_strat_genus_sum$pathway == "PWY-6572"), ]
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt <- melt(hmp2_16S_pathabun_strat_genus_sum_PWY_6572)
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_tmp <- hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_tmp <- hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_tmp[, -3]
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_genera <- aggregate(value ~ genus + pathway, data=hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_tmp, FUN=sum)
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_genera <- hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_genera[with(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_genera, order(value, decreasing = TRUE)),]
PWY_6572_top_genera <- head(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_genera$genus, 10)

hmp2_16S_pathabun_strat_genus_sum_PWY0_1533 <- hmp2_16S_pathabun_strat_genus_sum[which(hmp2_16S_pathabun_strat_genus_sum$pathway == "PWY0-1533"), ]
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt <- melt(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533)
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_tmp <- hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_tmp <- hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_tmp[, -3]
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_genera <- aggregate(value ~ genus + pathway, data=hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_tmp, FUN=sum)
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_genera <- hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_genera[with(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_genera, order(value, decreasing = TRUE)),]
PWY0_1533_top_genera <- head(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_genera$genus, 10)

overall_top_genera <- sort(unique(c(PWY0_1533_top_genera, PWY_6572_top_genera)))

#qual_col <- c('#e6194b', '#3cb44b', 'yellow', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 'greenyellow', '#fabebe', '#008080', '#e6beff',
#              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', 'royalblue1', 'grey')

qual_col <- c('#e6194b', '#3cb44b', 'yellow', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 'greenyellow', '#fabebe', '#008080', 'grey')

hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus_clean <- hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus_clean <- hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus

hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt[which(! hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus_clean %in%  PWY_6572_top_genera), "genus_clean"] <- "Other"
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt[which(! hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus_clean %in%  PWY0_1533_top_genera), "genus_clean"] <- "Other"

hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus_clean <- factor(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus_clean, levels=c(overall_top_genera, "Other"))
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus_clean <- factor(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus_clean, levels=c(overall_top_genera, "Other"))


# First get plot of only a shared legend.
hmp2_16S_pathabun_strat_genus_sum_TMP <- melt(hmp2_16S_pathabun_strat_genus_sum)
hmp2_16S_pathabun_strat_genus_sum_TMP <- hmp2_16S_pathabun_strat_genus_sum_TMP[which(hmp2_16S_pathabun_strat_genus_sum_TMP$pathway %in% c("PWY-6572", "PWY0-1533")),]
hmp2_16S_pathabun_strat_genus_sum_TMP[which(! hmp2_16S_pathabun_strat_genus_sum_TMP$genus %in%  overall_top_genera), "genus"] <- "Other"
hmp2_16S_pathabun_strat_genus_sum_TMP$genus <- factor(hmp2_16S_pathabun_strat_genus_sum_TMP$genus, levels=c(overall_top_genera, "Other"))
tmp_plot <- ggplot(hmp2_16S_pathabun_strat_genus_sum_TMP, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=qual_col) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_blank(),
        legend.text=element_text(size=8)) +
  labs(fill="Genus")
grobs <- ggplotGrob(tmp_plot)$grobs
stacked_legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


# Get samples ordered by total relative abundance.
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_sample <- aggregate(value ~ variable, data=hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt, FUN=sum)
hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_sample <- hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_sample[with(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_sample, order(value, decreasing = FALSE)),]

hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$variable <- factor(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$variable,
                                                                   levels=hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt_by_sample$variable)

PWY_6572_col <- qual_col[which(levels(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus_clean) %in% hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt$genus_clean)]

PWY_6572_stacked <- ggplot(hmp2_16S_pathabun_strat_genus_sum_PWY_6572_melt, aes(x=variable, y=value, fill=genus_clean)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=PWY_6572_col) + 
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.text=element_text(size=12),
                     axis.title=element_text(size=14),
                     plot.title = element_text(hjust=0.2, vjust=-10, face="bold")) +
  ylab("Relative Abundance (%)") +
  xlab("Sample") +
  ggtitle("PWY-6572: Chondroitin sulfate degradation I (bacterial)") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=FALSE)


  
# Get samples ordered by total relative abundance.
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_sample <- aggregate(value ~ variable, data=hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt, FUN=sum)
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_sample <- hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_sample[with(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_sample, order(value, decreasing = FALSE)),]
  
hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$variable <- factor(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$variable,
                                                                    levels=hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt_by_sample$variable)
  
PWY0_1533_col <- qual_col[which(levels(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus_clean) %in% hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt$genus_clean)]

PWY0_1533_stacked <- ggplot(hmp2_16S_pathabun_strat_genus_sum_PWY0_1533_melt, aes(x=variable, y=value, fill=genus_clean)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=PWY0_1533_col) + 
  theme_bw() + theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x=element_blank(),
                     axis.text=element_text(size=12),
                     axis.title=element_text(size=14),
                     plot.title = element_text(hjust=0.2, vjust=-10, face="bold")) +
  ylab("Relative Abundance (%)") +
  xlab("Sample") +
  ggtitle("PWY0-1533: Methylphosphonate degradation I") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=FALSE)


#20x12
### Plot final figure.
plot_grid(cd_sig_higher_ratio_plot,
          num_contrib_genera_mgs_vs_16S,
          PWY_5188_breakdown_scatterplot,
          PWY_6572_stacked,
          PWY0_1533_stacked,
          stacked_legend,
          labels=c("A", "B", "C", "D", "E", ""),
          nrow=2,
          ncol=3)
