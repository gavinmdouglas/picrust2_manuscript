### Supplementary figure comparing # pathways identified as a function of phylogenetic diversity of the samples.

rm(list=ls(all=TRUE))

library(ggplot2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/")

hmp_picrust1_KEGG <- read.table("working_tables/hmp_pathway_comparison/picrust1_kegg_pathways.biom.tsv",
                                skip=1, comment.char = "", quote="", stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_picrust2_KEGG <- read.table("working_tables/hmp_pathway_comparison/picrust2_kegg_pathways.tsv",
                                stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_picrust2_metacyc <- read.table("16S_validation/picrust2_out/hmp_picrust2_path_nsti2.0.tsv",
                                   stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_faith_pd <- read.table("working_tables/hmp_pathway_comparison/hmp_faiths_diversity.tsv",
                           header=T, sep='\t', row.names=1)

# Make sure that the subset of samples used elsewhere in the manuscript is used here:
setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")
hmp_ko_predictions <- read_in_ko_predictions("hmp")

samples2keep <- colnames(hmp_ko_predictions$all_kos_overlap$picrust2_ko_nsti2)

hmp_picrust1_KEGG <- hmp_picrust1_KEGG[, samples2keep]

hmp_picrust2_KEGG <- hmp_picrust2_KEGG[, samples2keep]

hmp_picrust2_metacyc <- hmp_picrust2_metacyc[, samples2keep]

hmp_faith_pd <- hmp_faith_pd[samples2keep, , drop=FALSE]

# Remove all PICRUSt1 pathways that are not possible in PICRUSt2.
kegg_pathway_descrip <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz"),
                                   header=FALSE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")

rownames(kegg_pathway_descrip) <- kegg_pathway_descrip$V2

hmp_picrust1_KEGG_subset <- hmp_picrust1_KEGG[-which(rowSums(hmp_picrust1_KEGG) == 0), ]
hmp_picrust1_KEGG_subset <- hmp_picrust1_KEGG_subset[which(rownames(hmp_picrust1_KEGG_subset) %in% rownames(kegg_pathway_descrip)), ]

hmp_picrust2_KEGG_subset <- hmp_picrust2_KEGG[-which(rowSums(hmp_picrust2_KEGG) == 0), ]


hmp_pathway_richness_df <- data.frame(sample=samples2keep,
                                      faith_pd=hmp_faith_pd$faith_pd,
                                      picrust1_kegg_richness=colSums(hmp_picrust1_KEGG_subset > 0),
                                      picrust2_kegg_richness=colSums(hmp_picrust2_KEGG_subset > 0),
                                      picrust2_metacyc_richness=colSums(hmp_picrust2_metacyc > 0))

hmp_pathway_richness_melt <- melt(hmp_pathway_richness_df, id.vars=c("sample", "faith_pd"))


picrust1_kegg_richness_plot <- ggplot(hmp_pathway_richness_melt[which(hmp_pathway_richness_melt$variable == "picrust1_kegg_richness"), ],
                                      aes(x=faith_pd, y=value)) +
                                      geom_point(size=1.5) + ylab("Number of Pathways") +
                                      xlab("Faith's Phylogenetic Diversity") +
                                      ggtitle("PICRUSt1 (KEGG)") +
                                      ylim(0, 350)

picrust2_kegg_richness_plot <- ggplot(hmp_pathway_richness_melt[which(hmp_pathway_richness_melt$variable == "picrust2_kegg_richness"), ],
                                      aes(x=faith_pd, y=value)) +
  geom_point(size=1.5) + ylab("Number of Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ggtitle("PICRUSt2 (KEGG)") +
  ylim(0, 350)


picrust2_metacyc_richness_plot <- ggplot(hmp_pathway_richness_melt[which(hmp_pathway_richness_melt$variable == "picrust2_metacyc_richness"), ],
                                      aes(x=faith_pd, y=value)) +
  geom_point(size=1.5) + ylab("Number of Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ggtitle("PICRUSt2 (MetaCyc)") +
  ylim(0, 350)

pdf(file = "../../figures/Supp_pathway_richness.pdf", width=12, height=4)

plot_grid(picrust1_kegg_richness_plot, picrust2_kegg_richness_plot, picrust2_metacyc_richness_plot, labels = c("a", "b", "c"), nrow=1)

dev.off()


# Calculation of how many more pathways are called by PICRUSt1 vs PICRUSt2.
mean(colSums(hmp_picrust1_KEGG_subset > 0) / colSums(hmp_picrust2_KEGG_subset > 0))

