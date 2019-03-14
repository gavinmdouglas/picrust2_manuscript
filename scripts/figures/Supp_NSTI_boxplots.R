### Exploring NSTI values in each 16S validation dataset (and in the 18S and ITS datasets).

setwd("/home/gavin/projects/picrust_pipeline/data/validation/")

library(ggplot2)
library(cowplot)
library(Biostrings)

# Read in HMP2 NSTI values, which were originally used to select the cut-off of 2.
hmp2_val_nsti <- read.table("/home/gavin/projects/hmp2_ibd_working/16S/April2018_redownload/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

boxplot(hmp2_val_nsti$metadata_NSTI, ylab="NSTI",  col="grey")
abline(h=2, lty=2, lwd=2)

# Read in NSTI values.
hmp_nsti <- read.table("hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

mammal_nsti <- read.table("iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blueberry_nsti <- read.table("blueberry/16S/picrust2_pipeline/picrust2_full_output_pipeline_2.1.0-b/marker_predicted_and_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

ocean_nsti <- read.table("ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blueberry_18S_nsti <- read.table("blueberry/18S/picrust2_full_output/marker_predicted_and_nsti.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_ITS_nsti <- read.table("wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/marker_predicted_and_nsti.tsv",
                             header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

# Read in weighted NSTI values.
hmp_nsti_weighted <- read.table("hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

mammal_nsti_weighted <- read.table("iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blueberry_nsti_weighted <- read.table("blueberry/16S/picrust2_pipeline/picrust2_full_output_pipeline_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

ocean_nsti_weighted <- read.table("ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blueberry_18S_nsti_weighted <- read.table("blueberry//18S/picrust2_full_output/ec_18S_counts_metagenome_out/weighted_nsti.tsv",
                                 header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

wine_ITS_nsti_weighted <- read.table("wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/ec_ITS_counts_metagenome_out//weighted_nsti.tsv",
                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

# Read in percent identity values.
hmp_percent_id <- read.table("hmp/16S/qiime2_artifacts/hmp_16S_rep_seqs_reference_align.txt",
                                header=F, sep="\t", stringsAsFactors = FALSE)

mammal_percent_id <- read.table("iGEM/16S/deblur_output_final/iGEM_16S_rep_seqs_reference_align.txt",
                             header=F, sep="\t", stringsAsFactors = FALSE)

blueberry_percent_id <- read.table("blueberry/16S/deblur_output_exported/blueberry_16S_rep_seqs_reference_align.txt",
                               header=F, sep="\t", stringsAsFactors = FALSE)

ocean_percent_id <- read.table("ocean/16S/deblur_output_final/ocean_16S_rep_seqs_reference_align.txt",
                             header=F, sep="\t", stringsAsFactors = FALSE)

blueberry_18S_percent_id <- read.table("blueberry/18S/deblur_output_exported/blueberry_18S_rep_seqs_align.txt",
                                       header=F, sep="\t", stringsAsFactors = FALSE)

wine_ITS_percent_id <- read.table("wine_fungi/ITS/ALL.BIT.derep.hash.n10.swarm.2.representative_renamed.align.txt",
                                  header=F, sep="\t", stringsAsFactors = FALSE)
missing_wine_seqs <- readDNAStringSet("wine_fungi/ITS/ALL.BIT.derep.hash.n10.swarm.2.representative_renamed.aligned.UNMATCHED.fasta")
unmatched_wine_percent_id <- data.frame(matrix(NA, nrow=length(missing_wine_seqs), ncol=12))
colnames(unmatched_wine_percent_id) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
unmatched_wine_percent_id$V1 <- names(missing_wine_seqs)
unmatched_wine_percent_id$V3 <- 0
wine_ITS_percent_id_all <- rbind(wine_ITS_percent_id, unmatched_wine_percent_id)

# Get distribution of wine ITS sequence lengths.
all_wine_seqs <- readDNAStringSet("wine_fungi/ITS/ALL.BIT.derep.hash.n10.swarm.2.representative_renamed.fasta")
all_wine_seq_lengths <- sapply(all_wine_seqs, length)

### Old code used to prep Soil data.
#Soil - 195 could not be aligned.
# soil_percent_id <- read.table("soil_crossbiome/16S/qiime2_artifacts/soil_16S_rep_seqs_reference_align.txt",
#                                header=F, sep="\t", stringsAsFactors = FALSE)
# 
# # Read in unmatched seqs to get their names.
# missing_soil_seqs <- readDNAStringSet("soil_crossbiome/16S/qiime2_artifacts/soil_16S_rep_seqs_ref_unmatched.fasta")
# 
# # Make empty dataframe for missing ASVs.
# unmatched_soil_percent_id <- data.frame(matrix(NA, nrow=length(missing_soil_seqs), ncol=12))
# colnames(unmatched_soil_percent_id) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
# unmatched_soil_percent_id$V1 <- names(missing_soil_seqs)
# unmatched_soil_percent_id$V3 <- 0
# 
# soil_percent_id_all <- rbind(soil_percent_id, unmatched_soil_percent_id)

# Read in biom tables to get abundance of excluded ASVs.
hmp_biom <- read.table("hmp/16S/qiime2_artifacts/hmp_16S_taxa.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

mammal_biom <- read.table("iGEM/16S/deblur_output_final/iGEM_16S_TAXA.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

ocean_biom <- read.table("ocean/16S/deblur_output_final/ocean_16S_TAXA.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

blueberry_biom <- read.table("blueberry/16S/deblur_output_exported/blueberry_16S_taxa.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

blueberry_18S_biom <- read.table("blueberry/18S/deblur_output_exported/blueberry_18S.biom.tsv",
                             header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

wine_ITS_biom <- read.table("wine_fungi/ITS/2014_swarm.2.otutable_renamed.txt",
                                 header=TRUE, sep="\t", stringsAsFactors = FALSE, comment.char="", row.names=1)

wine_ITS_biom <- wine_ITS_biom[, c("T1_D1.A", "T2_D1.A", "T2_D3.B", "Y1_D1.A", "Y2_D1.A", "Y2_D3.B", "Y3_D1.A", "Y3_D3.B")]

# Remove "ConsensusLineage" column.
hmp_biom <- hmp_biom[, -which(colnames(hmp_biom) == "ConsensusLineage")]
mammal_biom <- mammal_biom[, -which(colnames(mammal_biom) == "ConsensusLineage")]
ocean_biom <- ocean_biom[, -which(colnames(ocean_biom) == "ConsensusLineage")]
blueberry_biom <- blueberry_biom[, -which(colnames(blueberry_biom) == "ConsensusLineage")]

# Convert to relative abundances.
hmp_biom_relab <- data.frame(sweep(hmp_biom, 2, colSums(hmp_biom), FUN="/")) * 100
mammal_biom_relab <- data.frame(sweep(mammal_biom, 2, colSums(mammal_biom), FUN="/")) * 100
ocean_biom_relab <- data.frame(sweep(ocean_biom, 2, colSums(ocean_biom), FUN="/")) * 100
blueberry_biom_relab <- data.frame(sweep(blueberry_biom, 2, colSums(blueberry_biom), FUN="/")) * 100
blueberry_18S_biom_relab <- data.frame(sweep(blueberry_18S_biom, 2, colSums(blueberry_18S_biom), FUN="/")) * 100
wine_ITS_biom_relab <- data.frame(sweep(wine_ITS_biom, 2, colSums(wine_ITS_biom), FUN="/")) * 100

### Make combined dataframes of NSTI and percent identity per dataset.
# First get ASVs in same order between both tables.
hmp_nsti <- hmp_nsti[hmp_percent_id$V1,]
mammal_nsti <- mammal_nsti[mammal_percent_id$V1,]
ocean_nsti <- ocean_nsti[ocean_percent_id$V1,]
blueberry_nsti <- blueberry_nsti[blueberry_percent_id$V1,]
blueberry_18S_nsti <- blueberry_18S_nsti[blueberry_18S_percent_id$V1,]
wine_ITS_nsti <- wine_ITS_nsti[wine_ITS_percent_id_all$V1,]


hmp_nsti_id <- data.frame(asv=rownames(hmp_nsti), nsti=hmp_nsti$metadata_NSTI, percent_id=hmp_percent_id$V3, dataset="HMP", stringsAsFactors = FALSE)
mammal_nsti_id <- data.frame(asv=rownames(mammal_nsti), nsti=mammal_nsti$metadata_NSTI, percent_id=mammal_percent_id$V3, dataset="Mammal", stringsAsFactors = FALSE)
ocean_nsti_id <- data.frame(asv=rownames(ocean_nsti), nsti=ocean_nsti$metadata_NSTI, percent_id=ocean_percent_id$V3, dataset="Ocean", stringsAsFactors = FALSE)
blueberry_nsti_id <- data.frame(asv=rownames(blueberry_nsti), nsti=blueberry_nsti$metadata_NSTI, percent_id=blueberry_percent_id$V3, dataset="Soil ", stringsAsFactors = FALSE)
blueberry_18S_nsti_id <- data.frame(asv=rownames(blueberry_18S_nsti), nsti=blueberry_18S_nsti$metadata_NSTI, percent_id=blueberry_18S_percent_id$V3, dataset="Soil 18S ", stringsAsFactors = FALSE)
wine_ITS_nsti_id_all <- data.frame(asv=rownames(wine_ITS_nsti), nsti=wine_ITS_nsti$metadata_NSTI, percent_id=wine_ITS_percent_id_all$V3, dataset="Wine ITS", stringsAsFactors = FALSE)

combined_nsti_id <- rbind(hmp_nsti_id, mammal_nsti_id, ocean_nsti_id, blueberry_nsti_id, blueberry_18S_nsti_id, wine_ITS_nsti_id_all)

# Determine how many ASVs were excluded due to NSTI cut-off and what % of total relative abundance these correspond to.
hmp_excluded_asvs <- hmp_nsti_id[which(hmp_nsti_id$nsti > 2), "asv"]
hmp_excluded_asvs_summed_per <- sum(rowSums(hmp_biom_relab[hmp_excluded_asvs,]))
(length(hmp_excluded_asvs)/nrow(hmp_nsti_id))*100
(hmp_excluded_asvs_summed_per/sum(hmp_biom_relab))*100
length(hmp_excluded_asvs)

ocean_excluded_asvs <- ocean_nsti_id[which(ocean_nsti_id$nsti > 2), "asv"]
ocean_excluded_asvs_summed_per <- sum(rowSums(ocean_biom_relab[ocean_excluded_asvs,]))
(length(ocean_excluded_asvs)/nrow(ocean_nsti_id))*100
(ocean_excluded_asvs_summed_per/sum(ocean_biom_relab))*100
length(ocean_excluded_asvs)

length(which(mammal_nsti_id$nsti > 2))

blueberry_excluded_asvs <- blueberry_nsti_id[which(blueberry_nsti_id$nsti > 2), "asv"]
blueberry_excluded_asvs_summed_per <- sum(rowSums(blueberry_biom_relab[blueberry_excluded_asvs,]))
(length(blueberry_excluded_asvs)/nrow(blueberry_nsti_id))*100
(blueberry_excluded_asvs_summed_per/sum(blueberry_biom_relab))*100
length(blueberry_excluded_asvs)

blueberry_18S_excluded_asvs <- blueberry_18S_nsti_id[which(blueberry_18S_nsti_id$nsti > 2), "asv"]
blueberry_18S_excluded_asvs_summed_per <- sum(rowSums(blueberry_18S_biom_relab[blueberry_18S_excluded_asvs,]))
(length(blueberry_18S_excluded_asvs)/nrow(blueberry_18S_nsti_id))*100
(blueberry_18S_excluded_asvs_summed_per/sum(blueberry_18S_biom_relab))*100
length(blueberry_18S_excluded_asvs)

wine_ITS_excluded_asvs <- wine_ITS_nsti_id_all[which(wine_ITS_nsti_id_all$nsti > 2), "asv"]
wine_ITS_excluded_asvs_summed_per <- sum(rowSums(wine_ITS_biom_relab[wine_ITS_excluded_asvs,]))
(length(wine_ITS_excluded_asvs)/nrow(wine_ITS_nsti_id_all))*100
(wine_ITS_excluded_asvs_summed_per/sum(wine_ITS_biom_relab))*100
length(wine_ITS_excluded_asvs)

# HMP: 2 ASVs, corresponding to 0.107% of ASVs and 0.002% of all relative abundance.
# Mammal: 0 ASVs.
# Ocean: 3 ASVs, corresponding to 0.261% of ASVs and 0.025% of all relative abundance.
# Blueberry: 0 ASVs
# Blueberry 18S: 0 ASVs
# Wine ITS: 0 ASVs

hmp_nsti_weighted$dataset <- "HMP"
mammal_nsti_weighted$dataset <- "Mammal"
ocean_nsti_weighted$dataset <- "Ocean"
blueberry_nsti_weighted$dataset <- "Soil "
blueberry_18S_nsti_weighted$dataset <- "Soil 18S "
wine_ITS_nsti_weighted$dataset <- "Wine ITS"

combined_nsti_weighted <- rbind(hmp_nsti_weighted, mammal_nsti_weighted, ocean_nsti_weighted,
                                blueberry_nsti_weighted, blueberry_18S_nsti_weighted, wine_ITS_nsti_weighted)

combined_nsti_weighted$dataset <- factor(combined_nsti_weighted$dataset, levels=c("HMP", "Mammal", "Ocean",
                                                                                  "Soil ", "Soil 18S ",
                                                                                  "Wine ITS"))

combined_nsti_id$dataset  <- factor(combined_nsti_id$dataset, levels=c("HMP", "Mammal", "Ocean",
                                                                             "Soil ", "Soil 18S ",
                                                                             "Wine ITS"))

percent_id_boxplots <- ggplot(combined_nsti_id, aes(dataset, 100 - percent_id)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("100% - (% Identity)") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110)) 

full_nsti_boxplots <- ggplot(combined_nsti_id, aes(dataset, nsti)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("Nearest Sequenced Taxon Index") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + geom_hline(yintercept=c(2), linetype="dotted")

cropped_nsti_boxplots <- ggplot(combined_nsti_id, aes(dataset, nsti)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("Nearest Sequenced Taxon Index") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4)) + geom_hline(yintercept=c(2), linetype="dotted")


weighted_nsti_boxplots <- ggplot(combined_nsti_weighted, aes(dataset, weighted_NSTI)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("Weighted Nearest Sequenced Taxon Index") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

plot_grid(full_nsti_boxplots, cropped_nsti_boxplots, weighted_nsti_boxplots, percent_id_boxplots, labels = c("A", "B", "C", "D"), align="h", axis="b")






### Testing significance:

# Kruskal-Wallis test for significant differences in NSTI values.
combined_nsti_id_16S <- combined_nsti_id
combined_nsti_id_16S$dataset <- as.character(combined_nsti_id_16S$dataset)
combined_nsti_id_16S <- combined_nsti_id_16S[-which(combined_nsti_id_16S$dataset %in% c("Soil 18S ", "Wine ITS")), ]
combined_nsti_id_16S$dataset <- factor(combined_nsti_id_16S$dataset)

kruskal.test(dataset ~ nsti, data=combined_nsti_id_16S)
# Kruskal-Wallis chi-squared = 6504.3, df = 5615, p-value = 6.662e-16

# By % id
kruskal.test(dataset ~ percent_id, data=combined_nsti_id)
# Kruskal-Wallis chi-squared = 4040.2, df = 237, p-value < 2.2e-16

# By weighted NSTI
weighted_NSTI_16S <- combined_nsti_weighted
weighted_NSTI_16S$dataset <- as.character(weighted_NSTI_16S$dataset)
weighted_NSTI_16S <- weighted_NSTI_16S[-which(weighted_NSTI_16S$dataset %in% c("Soil 18S ", "Wine ITS")), ]
weighted_NSTI_16S$dataset <- factor(weighted_NSTI_16S$dataset)
kruskal.test(dataset ~ weighted_NSTI, data=weighted_NSTI_16S)
# Kruskal-Wallis chi-squared = 229, df = 229, p-value = 0.4876

# mean and sd NSTI values for each dataset:
# HMP - 0.1145611,  0.4894008
mean(hmp_nsti$metadata_NSTI)
sd(hmp_nsti$metadata_NSTI)

# iGEM - 0.1961959, 0.2161644
mean(mammal_nsti$metadata_NSTI)
sd(mammal_nsti$metadata_NSTI)

# Ocean - 0.5068952, 2.060702
mean(ocean_nsti$metadata_NSTI)
sd(ocean_nsti$metadata_NSTI)

# Blueberry Soil - 0.2760852, 0.2971092
mean(blueberry_nsti$metadata_NSTI)
sd(blueberry_nsti$metadata_NSTI)


# mean and sd % id values for each dataset:
# HMP - 97.25469,  4.358125
mean(hmp_percent_id$V3)
sd(hmp_percent_id$V3)

# Mammal - 91.7839, 4.56758
mean(mammal_percent_id$V3)
sd(mammal_percent_id$V3)

# Ocean - 90.20322, 6.797014
mean(ocean_percent_id$V3)
sd(ocean_percent_id$V3)

# Blueberry - 90.31872, 5.566207
mean(blueberry_percent_id$V3)
sd(blueberry_percent_id$V3)
