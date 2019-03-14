### Supplementary figure comparing # pathways identified as a function of phylogenetic diversity of the samples.

hmp_KEGG_pathways <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime_pipeline/final_otu_tables/otu_table_gg-only_filt_norm_ko_L3.biom.txt",
                                skip=1, comment.char = "", quote="", stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

colnames(hmp_KEGG_pathways) <- gsub(".nonchimera.fasta", "", colnames(hmp_KEGG_pathways))

hmp_KEGG_pathways <- hmp_KEGG_pathways[, -which(colnames(hmp_KEGG_pathways) == "KEGG_Pathways")]

hmp_metacyc_pathways <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/pathways_out/path_abun_unstrat.tsv",
                                   stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_faith_pd <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/diversity/faith_pd_diversity/exported/alpha-diversity.tsv",
                           header=T, sep='\t', row.names=1)


overlapping_col <- colnames(hmp_KEGG_pathways)[which(colnames(hmp_KEGG_pathways) %in% colnames(hmp_metacyc_pathways))]

hmp_faith_pd_subset <- hmp_faith_pd[overlapping_col,, drop=FALSE]

hmp_KEGG_pathways <- hmp_KEGG_pathways[, overlapping_col]
hmp_metacyc_pathways <- hmp_metacyc_pathways[, overlapping_col]

hmp_KEGG_pathways_relab <- data.frame(sweep(hmp_KEGG_pathways, 2, colSums(hmp_KEGG_pathways), '/')) * 100
hmp_metacyc_pathways_relab <- data.frame(sweep(hmp_metacyc_pathways, 2, colSums(hmp_metacyc_pathways), '/')) * 100

hmp_KEGG_pathways_relab <- hmp_KEGG_pathways_relab[-which(rowSums(hmp_KEGG_pathways_relab) == 0),]

hmp_faith_pd_vs_KEGG_richness <- data.frame(sample=overlapping_col,
                                            faith_pd=hmp_faith_pd_subset$faith_pd,
                                            kegg_richness=colSums(hmp_KEGG_pathways_relab > 0))

hmp_faith_pd_vs_KEGG_richness_melt <- melt(hmp_faith_pd_vs_KEGG_richness, id.vars=c("sample", "faith_pd"))

kegg_richness_plot <- ggplot(hmp_faith_pd_vs_KEGG_richness_melt, aes(x=faith_pd, y=value)) +
  geom_point(size=1.5) + ylab("Number of KEGG Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ylim(0, 350)

hmp_faith_pd_vs_MetaCyc_richness <- data.frame(sample=overlapping_col,
                                               faith_pd=hmp_faith_pd_subset$faith_pd,
                                               MetaCyc_richness=colSums(hmp_metacyc_pathways_relab > 0))

hmp_faith_pd_vs_MetaCyc_richness_melt <- melt(hmp_faith_pd_vs_MetaCyc_richness, id.vars=c("sample", "faith_pd"))

metacyc_richness_plot <- ggplot(hmp_faith_pd_vs_MetaCyc_richness_melt, aes(x=faith_pd, y=value)) +
  geom_point(size=1.5) + ylab("Number of MetaCyc Pathways") +
  xlab("Faith's Phylogenetic Diversity") +
  ylim(0, 350)
