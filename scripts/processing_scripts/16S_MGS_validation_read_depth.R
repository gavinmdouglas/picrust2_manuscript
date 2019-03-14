### Commands to get mean sample read depth for each dataset.

setwd("/home/gavin/projects/picrust_pipeline/data/validation/")

# Get KO spearman validation output (to get which samples to subset to)
hmp_ko_spearman <- readRDS("/home/gavin/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/hmp_ko_spearman_df.rds")
mammal_ko_spearman <- readRDS("/home/gavin/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/mammal_ko_spearman_df.rds")
ocean_ko_spearman <- readRDS("/home/gavin/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/ocean_ko_spearman_df.rds")
soil_ko_spearman <- readRDS("/home/gavin/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/soil_ko_spearman_df.rds")

# Read in biom tables to get abundance of excluded ASVs.
hmp_biom <- read.table("hmp/16S/qiime2_artifacts/hmp_16S_taxa.biom.tsv",
                       header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

mammal_biom <- read.table("iGEM/16S/deblur_output_final/iGEM_16S_TAXA.biom.tsv",
                          header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

ocean_biom <- read.table("ocean/16S/deblur_output_final/ocean_16S_TAXA.biom.tsv",
                         header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

soil_biom <- read.table("soil_crossbiome/16S/qiime2_artifacts/crossbiome_16S_TAXA.biom.tsv",
                        header=TRUE, sep="\t", stringsAsFactors = FALSE, skip=1, comment.char="", row.names=1)

# Subset to overlapping samples and get mean depth.
mean(colSums(hmp_biom[,levels(hmp_ko_spearman$sample_names)]))
mean(colSums(mammal_biom[,levels(mammal_ko_spearman$sample_names)]))
mean(colSums(ocean_biom[,levels(ocean_ko_spearman$sample_names)]))
mean(colSums(soil_biom[,levels(soil_ko_spearman$sample_names)]))

# Get # ASVs:
nrow(hmp_biom)
nrow(mammal_biom)
nrow(ocean_biom)
nrow(soil_biom)

# Get mean MGS read depth:
