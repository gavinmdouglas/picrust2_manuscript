library(ggplot2)
library(reshape2)
library(cowplot)
library(ggforce)

### Figure 1A will be flowchart figure (made in powerpoint).
### Figure 1B will be barplot with error bars showing % ASVs that match 100% with reference 16S sequences.
### Figure 1C and 1D will be scatterplots showing how pathway diversity changes with taxonomic diversity.

### Panels B-D will be saved as a single PDF.

# First make Figure 1B.
hmp_asv_blast <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/hmp_16S_rep_seqs_reference_align.txt",
                               header=FALSE, sep="\t", stringsAsFactors = FALSE)
hmp_otu_blast <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime_pipeline/clustering/rep_set_final_filt_revcomp.align.txt",
                           header=FALSE, sep="\t", stringsAsFactors = FALSE)

# Get % ASVs/OTUs 100% identical to reference sequences.
hmp_asv_abun <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/hmp_16S.biom.tsv",
                           header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")

hmp_otu_abun <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime_pipeline/final_otu_tables/otu_table_gg_10_2006.biom.tsv",
                           header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")
colnames(hmp_otu_abun) <- gsub(".nonchimera.fasta", "", colnames(hmp_otu_abun))

hmp_overlapping_samples <- colnames(hmp_asv_abun)[which(colnames(hmp_asv_abun) %in% colnames(hmp_otu_abun))]

hmp_asv_abun <- hmp_asv_abun[, hmp_overlapping_samples]
hmp_otu_abun <- hmp_otu_abun[, hmp_overlapping_samples]

hmp_asv_perfect_match <- c()
hmp_otu_perfect_match <- c()

for(hmp_sample in hmp_overlapping_samples) {
  
    asv_subset <- rownames(hmp_asv_abun)[which(hmp_asv_abun[, hmp_sample] > 0)]
    hmp_asv_blast_subset <- hmp_asv_blast[which(hmp_asv_blast$V1 %in% asv_subset),]
    hmp_asv_perfect_match <- c(hmp_asv_perfect_match, (length(which(hmp_asv_blast_subset$V3 == 100))/length(asv_subset)) * 100)

    otu_subset <- rownames(hmp_otu_abun)[which(hmp_otu_abun[, hmp_sample] > 0)]
    hmp_otu_blast_subset <- hmp_otu_blast[which(hmp_otu_blast$V1 %in% otu_subset),]
    hmp_otu_perfect_match <- c(hmp_otu_perfect_match, (length(which(hmp_otu_blast_subset$V3 == 100))/length(otu_subset)) * 100)

}

hmp_perfect_match <- data.frame(sample=hmp_overlapping_samples,
                                ASVs=hmp_asv_perfect_match,
                                OTUs=hmp_otu_perfect_match)

hmp_perfect_match_melt <- melt(hmp_perfect_match)

hmp_perfect_match_melt$variable <- factor(hmp_perfect_match_melt$variable, levels=c("OTUs", "ASVs"))

percent_matching_hmp_plot <- ggplot(hmp_perfect_match_melt, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA, fill="light grey") +
  geom_sina() + ylab("% Perfect Matches with Reference") + xlab("") +
  guides(fill=FALSE, colour=FALSE) + scale_fill_manual(values=c("light grey"))

### Figure 1C and 1D
### Comparing # pathways identified as a function of phylogenetic diversity of the samples.
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

plot_grid(percent_matching_hmp_plot, kegg_richness_plot, metacyc_richness_plot, labels = c("B", "C", "D"), ncol=3, nrow=1)