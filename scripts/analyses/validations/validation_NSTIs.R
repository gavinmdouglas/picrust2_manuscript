### Exploring NSTI values in each 16S validation dataset (and in the 18S and ITS datasets).

rm(list=ls(all=TRUE))

setwd("/home/gavin/projects/picrust_pipeline/data/validation/")

library(ggplot2)
library(cowplot)
library(Biostrings)

# Read in NSTI values.
cameroon_nsti <- read.table("cameroon/16S_workflow/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

hmp_nsti <- read.table("hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

indian_nsti <- read.table("indian/16S_workflow/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

mammal_nsti <- read.table("iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blueberry_nsti <- read.table("blueberry/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

ocean_nsti <- read.table("ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

primate_nsti <- read.table("primate/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/marker_predicted_and_nsti.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)


# Read in weighted NSTI values.
cameroon_nsti_weighted <- read.table("cameroon/16S_workflow/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

indian_nsti_weighted <- read.table("indian/16S_workflow/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                          header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

hmp_nsti_weighted <- read.table("hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

mammal_nsti_weighted <- read.table("iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blueberry_nsti_weighted <- read.table("blueberry/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

ocean_nsti_weighted <- read.table("ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

primate_nsti_weighted <- read.table("primate/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/KO_metagenome_out/weighted_nsti.tsv",
                           header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

# Read in percent identity values.
cameroon_percent_id <- read.table("cameroon/16S_workflow/deblur_output_final/cameroon_16S_rep_seqs_reference_align.txt",
                             header=F, sep="\t", stringsAsFactors = FALSE)

hmp_percent_id <- read.table("hmp/16S/qiime2_artifacts/hmp_16S_rep_seqs_reference_align.txt",
                                header=F, sep="\t", stringsAsFactors = FALSE)

indian_percent_id <- read.table("indian/16S_workflow/deblur_output_final/indian_16S_rep_seqs_reference_align.txt",
                                  header=F, sep="\t", stringsAsFactors = FALSE)

mammal_percent_id <- read.table("iGEM/16S/deblur_output_final/iGEM_16S_rep_seqs_reference_align.txt",
                             header=F, sep="\t", stringsAsFactors = FALSE)

blueberry_percent_id <- read.table("blueberry/16S/deblur_output_exported/blueberry_16S_rep_seqs_reference_align.txt",
                               header=F, sep="\t", stringsAsFactors = FALSE)

ocean_percent_id <- read.table("ocean/16S/deblur_output_final/ocean_16S_rep_seqs_reference_align.txt",
                             header=F, sep="\t", stringsAsFactors = FALSE)

primate_percent_id <- read.table("primate/16S/final_files/primate_16S_rep_seqs_reference_align.txt",
                                header=F, sep="\t", stringsAsFactors = FALSE)

missing_primate_seqs <- readDNAStringSet("primate/16S/final_files/primate_asvs.fna")
unmatched_primate_percent_id <- data.frame(matrix(NA, nrow=length(missing_primate_seqs), ncol=12))
colnames(unmatched_primate_percent_id) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
unmatched_primate_percent_id$V1 <- names(missing_primate_seqs)
unmatched_primate_percent_id$V3 <- 0
primate_percent_id_all <- rbind(primate_percent_id, unmatched_primate_percent_id)

# Read in biom tables to get abundance of excluded ASVs.
cameroon_biom <- read.table("cameroon/16S_workflow/deblur_output_final/cameroon_16S.biom.tsv",
                            header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

indian_biom <- read.table("indian/16S_workflow/deblur_output_final/indian_16S.biom.tsv",
                          header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

hmp_biom <- read.table("hmp/16S/qiime2_artifacts/hmp_16S.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

mammal_biom <- read.table("iGEM/16S/deblur_output_final/iGEM_16S.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

ocean_biom <- read.table("ocean/16S/deblur_output_final/ocean_16S.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

blueberry_biom <- read.table("blueberry/16S/deblur_output_exported/blueberry_16S.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

primate_biom <- read.table("primate/16S/final_files/primate_16S.biom.tsv",
                           header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

# Convert to relative abundances.
cameroon_biom_relab <- data.frame(sweep(cameroon_biom, 2, colSums(cameroon_biom), FUN="/")) * 100
indian_biom_relab <- data.frame(sweep(indian_biom, 2, colSums(indian_biom), FUN="/")) * 100
hmp_biom_relab <- data.frame(sweep(hmp_biom, 2, colSums(hmp_biom), FUN="/")) * 100
mammal_biom_relab <- data.frame(sweep(mammal_biom, 2, colSums(mammal_biom), FUN="/")) * 100
ocean_biom_relab <- data.frame(sweep(ocean_biom, 2, colSums(ocean_biom), FUN="/")) * 100
blueberry_biom_relab <- data.frame(sweep(blueberry_biom, 2, colSums(blueberry_biom), FUN="/")) * 100
primate_biom_relab <- data.frame(sweep(primate_biom, 2, colSums(primate_biom), FUN="/")) * 100

### Make combined dataframes of NSTI and percent identity per dataset.
# First get ASVs in same order between both tables.
cameroon_nsti <- cameroon_nsti[cameroon_percent_id$V1,]
indian_nsti <- indian_nsti[indian_percent_id$V1,]
hmp_nsti <- hmp_nsti[hmp_percent_id$V1,]
mammal_nsti <- mammal_nsti[mammal_percent_id$V1,]
ocean_nsti <- ocean_nsti[ocean_percent_id$V1,]
blueberry_nsti <- blueberry_nsti[blueberry_percent_id$V1,]
primate_nsti <- primate_nsti[primate_percent_id$V1,]

cameroon_nsti_id <- data.frame(asv=rownames(cameroon_nsti), nsti=cameroon_nsti$metadata_NSTI, percent_id=cameroon_percent_id$V3, dataset="Cameroonian", stringsAsFactors = FALSE)
indian_nsti_id <- data.frame(asv=rownames(indian_nsti), nsti=indian_nsti$metadata_NSTI, percent_id=indian_percent_id$V3, dataset="Indian", stringsAsFactors = FALSE)
hmp_nsti_id <- data.frame(asv=rownames(hmp_nsti), nsti=hmp_nsti$metadata_NSTI, percent_id=hmp_percent_id$V3, dataset="HMP", stringsAsFactors = FALSE)
mammal_nsti_id <- data.frame(asv=rownames(mammal_nsti), nsti=mammal_nsti$metadata_NSTI, percent_id=mammal_percent_id$V3, dataset="Mammal", stringsAsFactors = FALSE)
ocean_nsti_id <- data.frame(asv=rownames(ocean_nsti), nsti=ocean_nsti$metadata_NSTI, percent_id=ocean_percent_id$V3, dataset="Ocean", stringsAsFactors = FALSE)
blueberry_nsti_id <- data.frame(asv=rownames(blueberry_nsti), nsti=blueberry_nsti$metadata_NSTI, percent_id=blueberry_percent_id$V3, dataset="Soil (Blueberry)", stringsAsFactors = FALSE)
primate_nsti_id <- data.frame(asv=rownames(primate_nsti), nsti=primate_nsti$metadata_NSTI, percent_id=primate_percent_id$V3, dataset="Primate", stringsAsFactors = FALSE)

combined_nsti_id <- rbind(cameroon_nsti_id, indian_nsti_id, hmp_nsti_id, mammal_nsti_id, ocean_nsti_id, blueberry_nsti_id, primate_nsti_id)

# Determine how many ASVs were excluded due to NSTI cut-off and what % of total relative abundance these correspond to.
cameroon_excluded_asvs <- cameroon_nsti_id[which(cameroon_nsti_id$nsti > 2), "asv"]
cameroon_excluded_asvs_summed_per <- sum(rowSums(cameroon_biom_relab[cameroon_excluded_asvs,]))
(length(cameroon_excluded_asvs)/nrow(cameroon_nsti_id))*100
(cameroon_excluded_asvs_summed_per/sum(cameroon_biom_relab))*100
length(cameroon_excluded_asvs)


indian_excluded_asvs <- indian_nsti_id[which(indian_nsti_id$nsti > 2), "asv"]
indian_excluded_asvs_summed_per <- sum(rowSums(indian_biom_relab[indian_excluded_asvs,]))
(length(indian_excluded_asvs)/nrow(indian_nsti_id))*100
(indian_excluded_asvs_summed_per/sum(indian_biom_relab))*100
length(indian_excluded_asvs)


hmp_excluded_asvs <- hmp_nsti_id[which(hmp_nsti_id$nsti > 2), "asv"]
hmp_excluded_asvs_summed_per <- sum(rowSums(hmp_biom_relab[hmp_excluded_asvs,]))
(length(hmp_excluded_asvs)/nrow(hmp_nsti_id))*100
(hmp_excluded_asvs_summed_per/sum(hmp_biom_relab))*100
length(hmp_excluded_asvs)

mammal_excluded_asvs <- mammal_nsti_id[which(mammal_nsti_id$nsti > 2), "asv"]
mammal_excluded_asvs_summed_per <- sum(rowSums(mammal_biom_relab[mammal_excluded_asvs,]))
(length(mammal_excluded_asvs)/nrow(mammal_nsti_id))*100
(mammal_excluded_asvs_summed_per/sum(mammal_biom_relab))*100
length(mammal_excluded_asvs)

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

primate_excluded_asvs <- primate_nsti_id[which(primate_nsti_id$nsti > 2), "asv"]
primate_excluded_asvs_summed_per <- sum(rowSums(primate_biom_relab[primate_excluded_asvs,]))
(length(primate_excluded_asvs)/nrow(primate_nsti_id))*100
(primate_excluded_asvs_summed_per/sum(primate_biom_relab))*100
length(primate_excluded_asvs)

mean((length(cameroon_excluded_asvs)/nrow(cameroon_nsti_id))*100,
     (length(indian_excluded_asvs)/nrow(indian_nsti_id))*100,
     (length(hmp_excluded_asvs)/nrow(hmp_nsti_id))*100,
     (length(mammal_excluded_asvs)/nrow(mammal_nsti_id))*100,
     (length(ocean_excluded_asvs)/nrow(ocean_nsti_id))*100,
     (length(blueberry_excluded_asvs)/nrow(blueberry_nsti_id))*100,
     (length(primate_excluded_asvs)/nrow(primate_nsti_id))*100)
     
mean((hmp_excluded_asvs_summed_per/sum(hmp_biom_relab))*100,
     (mammal_excluded_asvs_summed_per/sum(mammal_biom_relab))*100,
     (ocean_excluded_asvs_summed_per/sum(ocean_biom_relab))*100,
     (blueberry_excluded_asvs_summed_per/sum(blueberry_biom_relab))*100)

cameroon_nsti_weighted$dataset <- "Cameroonian"
indian_nsti_weighted$dataset <- "Indian"
hmp_nsti_weighted$dataset <- "HMP"
mammal_nsti_weighted$dataset <- "Mammal"
ocean_nsti_weighted$dataset <- "Ocean"
blueberry_nsti_weighted$dataset <- "Soil (Blueberry)"
primate_nsti_weighted$dataset <- "Primate"

combined_nsti_weighted <- rbind(cameroon_nsti_weighted, indian_nsti_weighted, hmp_nsti_weighted, mammal_nsti_weighted, ocean_nsti_weighted,
                                blueberry_nsti_weighted, primate_nsti_weighted)

combined_nsti_weighted$dataset <- factor(combined_nsti_weighted$dataset, levels=c("Cameroonian", "HMP",  "Indian",
                                                                                  "Mammal", "Ocean", "Primate" "Soil (Blueberry)" ))

combined_nsti_id$dataset  <- factor(combined_nsti_id$dataset, levels=c("Cameroonian", "HMP", "Indian", 
                                                                       "Mammal", "Ocean", "Primate", "Soil (Blueberry)" ))

percent_id_boxplots <- ggplot(combined_nsti_id, aes(dataset, 100 - percent_id)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("100% - (% Identity)") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110))

full_nsti_boxplots <- ggplot(combined_nsti_id, aes(dataset, nsti)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("Nearest Sequenced Taxon Index") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) + geom_hline(yintercept=c(2), linetype="dotted")

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

pdf(file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/figures/Supp_NSTI_boxplots.pdf", width=12, height=8)

plot_grid(full_nsti_boxplots, cropped_nsti_boxplots, weighted_nsti_boxplots, percent_id_boxplots,
          labels = c("a", "b", "c", "d"), align="h", axis="b")

dev.off()

# Kruskal-Wallis test for significant differences in NSTI values.
kruskal.test(dataset ~ nsti, data=combined_nsti_id)
# Kruskal-Wallis chi-squared = 20804, df = 16148, p-value < 2.2e-16

# By % id
kruskal.test(dataset ~ percent_id, data=combined_nsti_id)
# Kruskal-Wallis chi-squared = 9553.1, df = 334, p-value < 2.2e-16

# By weighted NSTI
kruskal.test(dataset ~ weighted_NSTI, data=combined_nsti_weighted)
# Kruskal-Wallis chi-squared = 560, df = 560, p-value = 0.4921

# mean and sd NSTI values for each dataset:

# Cameroon - 0.3713435, 0.3483202
mean(cameroon_nsti$metadata_NSTI)
sd(cameroon_nsti$metadata_NSTI)

# Indian - 0.09970479, 0.1117562
mean(indian_nsti$metadata_NSTI)
sd(indian_nsti$metadata_NSTI)

# HMP - 0.1145611,  0.4894008
mean(hmp_nsti$metadata_NSTI)
sd(hmp_nsti$metadata_NSTI)

# iGEM - 0.1961959, 0.2161644
mean(mammal_nsti$metadata_NSTI)
sd(mammal_nsti$metadata_NSTI)

# Ocean - 0.5068952, 2.060702
mean(ocean_nsti$metadata_NSTI)
sd(ocean_nsti$metadata_NSTI)

# Blueberry Soil - 0.3193135, 0.2254971
mean(blueberry_nsti$metadata_NSTI)
sd(blueberry_nsti$metadata_NSTI)

# Primate - 0.265337, 1.85833
mean(primate_nsti$metadata_NSTI)
sd(primate_nsti$metadata_NSTI)

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
