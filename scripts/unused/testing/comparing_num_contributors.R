ocean_denovo_kos <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_predicted.tsv",
                             header=T, sep="\t", row.names=1)

ocean_gg_kos <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/KO_predicted.tsv",
                             header=T, sep="\t", row.names=1)

tmp1 <- colSums(ocean_denovo_kos > 0)
tmp2 <- colSums(ocean_gg_kos > 0)

row_sum1 <- rowSums(ocean_denovo_kos)
row_sum2 <- rowSums(ocean_gg_kos)

ocean_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

ocean_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

ocean_denovo_marker_filt <- ocean_denovo_marker[-which(ocean_denovo_marker$metadata_NSTI > 2),]
ocean_gg_marker_filt <- ocean_gg_marker[-which(ocean_gg_marker$metadata_NSTI > 2),]

ocean_denovo_kos_filt <- ocean_denovo_kos[rownames(ocean_denovo_marker_filt),]
ocean_gg_kos_filt <- ocean_gg_kos[rownames(ocean_gg_marker_filt),]

boxplot()

ocean_ratio_num_contrib <- (colSums(ocean_denovo_kos_filt > 0) + 1)/(colSums(ocean_gg_kos_filt > 0) + 1)

hmp_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                                   header=T, sep="\t", row.names=1)

hmp_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                               header=T, sep="\t", row.names=1)

mammal_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

mammal_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                               header=T, sep="\t", row.names=1)

boxplot(hmp_denovo_marker$metadata_NSTI, hmp_gg_marker$metadata_NSTI)

boxplot(mammal_denovo_marker$metadata_NSTI, mammal_gg_marker$metadata_NSTI)

boxplot(ocean_denovo_marker$metadata_NSTI, ocean_gg_marker$metadata_NSTI, ylim=c(0, 2))



blueberry_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_pipeline_2.1.0-b/marker_predicted_and_nsti.tsv",
                                   header=T, sep="\t", row.names=1)

blueberry_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                               header=T, sep="\t", row.names=1)

boxplot(ocean_denovo_marker$metadata_NSTI, ocean_gg_marker$metadata_NSTI, ylim=c(0, 2))
boxplot(blueberry_denovo_marker$metadata_NSTI, blueberry_gg_marker$metadata_NSTI, ylim=c(0, 2))


setwd("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/")

mammal_denovo_pathabun_strat <- read.table("picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                              header=T, sep="\t", stringsAsFactors = FALSE)

mammal_gg_pathabun_strat <- read.table("picrust2_full_output_picrust1-prep_2.1.0-b/pathways_out_per_seq/path_abun_strat.tsv",
                                              header=T, sep="\t", stringsAsFactors = FALSE)


mammal_denovo_pathabun_strat_subset <- mammal_denovo_pathabun_strat[, c("pathway", "sequence", "A1")]
mammal_gg_pathabun_strat_subset <- mammal_gg_pathabun_strat[, c("pathway", "sequence", "A1")]

mammal_denovo_pathabun_strat_subset <- mammal_denovo_pathabun_strat_subset[-which(mammal_denovo_pathabun_strat_subset$A1 == 0),]
mammal_gg_pathabun_strat_subset <- mammal_gg_pathabun_strat_subset[-which(mammal_gg_pathabun_strat_subset$A1 == 0),]


mammal_denovo_pathabun_strat_subset_path_table <- table(mammal_denovo_pathabun_strat_subset$pathway)

mammal_gg_pathabun_strat_subset_path_table <- table(mammal_gg_pathabun_strat_subset$pathway)

all_path <- unique(c(mammal_denovo_pathabun_strat_subset$pathway, mammal_gg_pathabun_strat_subset$pathway))

pathway_contributors <- data.frame(matrix(NA, nrow=length(all_path), ncol=2))
colnames(pathway_contributors) <- c("denovo", "gg")
rownames(pathway_contributors) <- all_path

for(pathway in all_path) {
  pathway_contributors[pathway,] <- c(length(which(mammal_denovo_pathabun_strat_subset$pathway==pathway)),
                                      length(which(mammal_gg_pathabun_strat_subset$pathway==pathway)))
}



hmp_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

hmp_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                              header=T, sep="\t", row.names=1)



mammal_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

mammal_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                              header=T, sep="\t", row.names=1)



ocean_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

ocean_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                              header=T, sep="\t", row.names=1)

blueberry_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_pipeline_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)

blueberry_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                              header=T, sep="\t", row.names=1)

boxplot(hmp_denovo_marker$metadata_NSTI, hmp_gg_marker$metadata_NSTI,
        mammal_denovo_marker$metadata_NSTI, mammal_gg_marker$metadata_NSTI,
        ocean_denovo_marker$metadata_NSTI, ocean_gg_marker$metadata_NSTI,
        blueberry_denovo_marker$metadata_NSTI, blueberry_gg_marker$metadata_NSTI,
        ylim=c(0,2),
        names=c("HMP Denovo", "HMP GG", "Mammal Denovo", "Mammal GG",
                "Ocean Denovo", "Ocean GG", "Blueberry Denovo", "Blueberry GG"))


mammal_denovo_percent <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/deblur_output_final/iGEM_16S_rep_seqs_reference_align.txt",
                                header=F, sep="\t", stringsAsFactors = FALSE)

mammal_gg_percent <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/qiime_pipeline/clustering/rep_set_final_filt_ref_align.txt",
                                header=F, sep="\t", stringsAsFactors = FALSE)

rownames(mammal_gg_percent) <- mammal_gg_percent$V1
rownames(mammal_denovo_percent) <- mammal_denovo_percent$V1

mammal_gg_marker$percent_id <- mammal_gg_percent[rownames(mammal_gg_marker), "V3"]
mammal_denovo_marker$percent_id <- mammal_denovo_percent[rownames(mammal_denovo_marker), "V3"]

par(mfrow=c(2,1))
plot(mammal_gg_marker$metadata_NSTI, mammal_gg_marker$percent_id)
plot(mammal_denovo_marker$metadata_NSTI, mammal_denovo_marker$percent_id)


hmp_gg_tax_biom <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime_pipeline/final_otu_tables/otu_table_gg_10_2006.biom.tsv",
                              header=T, stringsAsFactors = FALSE, skip=1, comment.char="", sep="\t")

hmp_gg_tax_biom_species <- hmp_gg_tax_biom[grep("s__", hmp_gg_tax_biom$taxonomy),]
hmp_gg_tax_biom_species <- hmp_gg_tax_biom_species[-grep("s__$", hmp_gg_tax_biom_species$taxonomy),]


hmp_gg_tax_biom_species[which(hmp_gg_tax_biom_species$taxonomy == "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Dorea; s__formicigenerans"),,drop=FALSE]