### Get basic summaries of all 16S validation datasets to report after removing all samples that do not overlap with MGS data.
### Note that the crossbiome and microbial mat datasets were included for curiousity, but we were not included in the manuscript.

rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

datasets <- c("hmp", "mammal", "ocean", "blueberry", "crossbiome", "mat", "indian", "cameroon", "primate")

ko_infiles <- list()

for(dataset in datasets) {
  ko_infiles[[dataset]] <- read_in_ko_predictions(dataset)
}

setwd("/home/gavin/projects/picrust_pipeline/data/validation/")

# Read in biom tables
biom_tables <- list()
biom_tables[["cameroon"]] <- read.table("cameroon/16S_workflow/deblur_output_final/cameroon_16S.renamed.biom.tsv",
                                    header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["indian"]] <- read.table("indian/16S_workflow/deblur_output_final/indian_16S.renamed.biom.tsv",
                          header=TRUE, sep="\t", stringsAsFactors = FALSE, comment.char="", row.names=1)

biom_tables[["hmp"]] <- read.table("hmp/16S/qiime2_artifacts/hmp_16S.biom.tsv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["mammal"]] <- read.table("iGEM/16S/deblur_output_final/iGEM_16S.biom.tsv",
                          header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["ocean"]] <- read.table("ocean/16S/deblur_output_final/ocean_16S.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["blueberry"]] <- read.table("blueberry/16S/deblur_output_exported/blueberry_16S.biom.tsv",
                             header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["primate"]] <- read.table("primate/16S/final_files/primate_16S.biom.tsv",
                           header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["mat"]] <- read.table("microbial_mat/QIITA_16S/BIOM/47904/otu_table.biom.tsv",
                                       header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

biom_tables[["crossbiome"]] <- read.table("soil_crossbiome/16S/qiime2_artifacts/soil_16S.biom.tsv",
                                   header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)


# Subset to samples in final prediction tables (that have been intersected with MGS).

data_summary <- data.frame(matrix(NA, nrow=9, ncol=5))
rownames(data_summary) <- datasets
colnames(data_summary) <- c("n", "num_asvs", "mean_depth", "min_depth", "max_depth")

for(dataset in datasets) {

  biom_subset <- biom_tables[[dataset]][, colnames(ko_infiles[[dataset]]$all_kos_overlap$picrust2_ko_nsti2)]
  
  empty_rows <- which(rowSums(biom_subset) == 0)
  
  if(length(empty_rows) > 0) {
    biom_subset <- biom_subset[-empty_rows, ] 
  }

  data_summary[dataset, ] <- c(ncol(biom_subset), nrow(biom_subset), mean(colSums(biom_subset)), min(colSums(biom_subset)), max(colSums(biom_subset)))

}

