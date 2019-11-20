### Exploring NSTI values in each 16S validation dataset.

rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_datasets/")

library(ggplot2)
library(cowplot)
library(Biostrings)

datasets <- c("cameroon", "primate", "hmp", "mammal", "ocean", "blueberry", "indian")

dataset2name <- list("cameroon"="Cameroonian", "indian"="Indian", "hmp"="HMP", "mammal"="Mammal",
                     "ocean"="Ocean", "blueberry"="Soil (Blueberry)", "primate"="Primate")

nsti <- list()
weighted <- list()
percentID <- list()
biom <- list()
nsti_percentID <- list()

# Also read in KO spearman output to make sure that only ASVs in the samples actually used
# for validations (i.e. that overlap with MGS data) are used.
ko_spearman_out <- list()

for(d in datasets) {

  ko_spearman_rds <- paste("../saved_RDS/16S_vs_MGS_metrics/", d, "_ko_spearman_df.rds", sep="")
  ko_spearman_out[[d]] <- readRDS(ko_spearman_rds)
  
  sample_names <- levels(ko_spearman_out[[d]]$sample_names)
  
  biom_infile <- paste(d, "/", d, "_16S.biom.tsv", sep="")
  nsti_infile <- paste(d, "marker_predicted_and_nsti.tsv", sep="/")
  weighted_infile <- paste(d, "weighted_nsti.tsv", sep="/")
  percentID_infile <- paste(d, "/", d, "_16S_rep_seqs_reference_align.txt", sep="")
  
  nsti[[d]] <- read.table(nsti_infile, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
  weighted[[d]] <- read.table(weighted_infile, header=TRUE, sep="\t", stringsAsFactors = FALSE, row.names=1)
  percentID[[d]] <- read.table(percentID_infile, header=FALSE, sep="\t", stringsAsFactors = FALSE)
  biom[[d]] <- read.table(biom_infile, skip=1, comment.char="", header=TRUE, sep="\t", row.names=1)
  
  # First figure out if any ASVs should be excluded, because they are not found in any samples actually used for validations.
  biom[[d]] <- biom[[d]][, sample_names]
  
  zero_rows <- which(rowSums(biom[[d]]) == 0)
  
  if(length(zero_rows) > 0) {
    biom[[d]] <- biom[[d]][-zero_rows, ]
  }
  
  asv2keep <- rownames(biom[[d]])
  
  # Subset NSTI and weighted NSTI tables to the ASVs and samples used for validations.
  nsti[[d]] <- nsti[[d]][asv2keep, ]
  
  # Convert weighted table to matrix and back to make sure that sample names are interpreted the as in RDS
  # (e.g. add X if starts with number and convert many characters to ".")
  weighted[[d]] <- data.frame(t(data.frame(t(weighted[[d]]))))

  weighted[[d]] <- weighted[[d]][sample_names, , drop=FALSE]
  weighted[[d]]$dataset <- dataset2name[[d]]
  
  
  # Subset BLAST6OUT table as well, but also add in any ASVs that might be missing (set as NA).
  missing_seqs <- asv2keep[which(! asv2keep %in% percentID[[d]]$V1)]
  
  if(length(missing_seqs) > 0) {
    unmatched_percent_id <- data.frame(matrix(NA, nrow=length(missing_seqs), ncol=12))
    colnames(unmatched_percent_id) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
    unmatched_percent_id$V1 <- missing_seqs
    unmatched_percent_id$V3 <- NA
    percentID[[d]] <- rbind(percentID[[d]], unmatched_percent_id)
  }
  
  rownames(percentID[[d]]) <- percentID[[d]]$V1
  
  percentID[[d]] <- percentID[[d]][asv2keep, ]
  
  nsti_percentID[[d]] <- data.frame(asv=asv2keep, nsti=nsti[[d]]$metadata_NSTI, percent_id=percentID[[d]]$V3, dataset=dataset2name[[d]], stringsAsFactors = FALSE)

}


combined_nsti_id <- do.call("rbind", nsti_percentID)

combined_nsti_weighted <- do.call("rbind", weighted)

combined_nsti_id$dataset  <- factor(combined_nsti_id$dataset, levels=c("Cameroonian", "HMP", "Indian", 
                                                                       "Mammal", "Ocean", "Primate", "Soil (Blueberry)"))

combined_nsti_weighted$dataset <- factor(combined_nsti_weighted$dataset, levels=c("Cameroonian", "HMP",  "Indian",
                                                                                  "Mammal", "Ocean", "Primate", "Soil (Blueberry)"))

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

percent_id_boxplots <- ggplot(combined_nsti_id, aes(dataset, 100 - percent_id)) + geom_boxplot(fill="light grey") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ylab("100% - (% Identity)") + xlab("Dataset") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

pdf(file = "../../figures/Supp_NSTI_boxplots.pdf", width=12, height=8)

plot_grid(full_nsti_boxplots, cropped_nsti_boxplots, weighted_nsti_boxplots, percent_id_boxplots,
          labels = c("a", "b", "c", "d"), align="h", axis="b")

dev.off()

# Kruskal-Wallis test for significant differences in NSTI values.
kruskal.test(dataset ~ nsti, data=combined_nsti_id)
# Kruskal-Wallis chi-squared = 19499, df = 15203, p-value < 2.2e-16

# By % id
kruskal.test(dataset ~ percent_id, data=combined_nsti_id)
# Kruskal-Wallis chi-squared = 9284.1, df = 320, p-value < 2.2e-16

# By weighted NSTI
kruskal.test(dataset ~ weighted_NSTI, data=combined_nsti_weighted)
# Kruskal-Wallis chi-squared = 400, df = 400, p-value = 0.4906


# Calculate basic summary statistics.
dataset_nsti_stats <- data.frame(matrix(NA, nrow=7, ncol=5))
colnames(dataset_nsti_stats) <- c("nsti_mean", "nsti_sd", "num_asvs", "num_asvs_excluded", "percent_asvs_excluded")
rownames(dataset_nsti_stats) <- datasets

for(d in datasets) {
  
  d_name <- dataset2name[[d]]
  
  dataset_nsti_stats[d, ] <- c(mean(nsti[[d]]$metadata_NSTI),
                               sd(nsti[[d]]$metadata_NSTI),
                               nrow(nsti[[d]]),
                               length(which(nsti[[d]]$metadata_NSTI > 2)),
                               (length(which(nsti[[d]]$metadata_NSTI > 2)) / nrow(nsti[[d]])) * 100)
  
}

# Get mean and sd of % ASVs excluded across all datasets:
mean(dataset_nsti_stats$percent_asvs_excluded)
sd(dataset_nsti_stats$percent_asvs_excluded)
