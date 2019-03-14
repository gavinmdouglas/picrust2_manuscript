### Comparing # pathways identified as a function of phylogenetic diversity of the samples.

rm(list=ls())

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/path_richness_v_faithdiv")

library(ggplot2)
library(cowplot)

# Read in input tables.
hmp_picrust1_KEGG_pathways <- read.table("otu_table_gg-only_filt_norm_ko_L3.biom.txt", 
                                skip=1, comment.char = "", quote="", stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

colnames(hmp_picrust1_KEGG_pathways) <- gsub(".nonchimera.fasta", "", colnames(hmp_picrust1_KEGG_pathways))

hmp_picrust1_KEGG_pathways <- hmp_picrust1_KEGG_pathways[, -which(colnames(hmp_picrust1_KEGG_pathways) == "KEGG_Pathways")]

kegg_pathway_descrip <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz"),
                                   header=FALSE, sep="\t", stringsAsFactors = FALSE, quote="")

# Remove any pathways in the PICRUSt1 output that aren't listed in the database used with PICRUSt2.
hmp_picrust1_KEGG_pathways <- hmp_picrust1_KEGG_pathways[-which(! rownames(hmp_picrust1_KEGG_pathways) %in% kegg_pathway_descrip$V2), ]

hmp_metacyc_pathways <- read.table("metacyc_path_abun_unstrat.tsv", stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_picrust2_KEGG_pathways <- read.table("KEGG_path_abun_unstrat.tsv", stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_faith_pd <- read.table("alpha-diversity.tsv", header=T, sep='\t', row.names=1)


# Subset tables to overlapping samples.
overlapping_col <- colnames(hmp_picrust1_KEGG_pathways)[which(colnames(hmp_picrust1_KEGG_pathways) %in% colnames(hmp_metacyc_pathways))]

hmp_faith_pd_subset <- hmp_faith_pd[overlapping_col,, drop=FALSE]
hmp_picrust1_KEGG_pathways <- hmp_picrust1_KEGG_pathways[, overlapping_col]
hmp_picrust2_KEGG_pathways <- hmp_picrust2_KEGG_pathways[, overlapping_col]
hmp_metacyc_pathways <- hmp_metacyc_pathways[, overlapping_col]

hmp_faith_pd_pathways <- data.frame(faith=hmp_faith_pd_subset$faith_pd,
                                    KEGG_picrust1=colSums(hmp_picrust1_KEGG_pathways > 0),
                                    KEGG_picrust2=colSums(hmp_picrust2_KEGG_pathways > 0),
                                    metacyc=colSums(hmp_metacyc_pathways > 0))


kegg_picrust1_pathways <- ggplot(data=hmp_faith_pd_pathways, aes(x=faith, y=KEGG_picrust1)) +
  geom_point(size=2) +
  ylab("Number of Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ggtitle("PICRUSt1 (KEGG)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350))

kegg_picrust2_pathways <- ggplot(data=hmp_faith_pd_pathways, aes(x=faith, y=KEGG_picrust2)) +
  geom_point(size=2) +
  ylab("Number of Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ggtitle("PICRUSt2 (KEGG)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350))


metacyc_pathways <- ggplot(data=hmp_faith_pd_pathways, aes(x=faith, y=metacyc)) +
  geom_point(size=2) +
  ylab("Number of Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ggtitle("PICRUSt2 (MetaCyc)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 350))

plot_grid(kegg_picrust1_pathways, kegg_picrust2_pathways, metacyc_pathways, ncol = 3, labels=c("A", "B", "C"))
