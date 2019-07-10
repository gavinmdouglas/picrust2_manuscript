rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

library(ggplot2)
library(stringi)
library(reshape2)
library(cowplot)
library(ggbeeswarm)

# Read in hmp2 metadata table.
hmp2_metadata <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_2018-08-20_col_subset.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE, check.names = TRUE)

# Read in stratified pathways (--per-seq-contrib option used).
pathabun_in_strat <- readRDS("hmp2_pathabun_strat_filt_count_Ileum.rds")

# Transpose stratified table and re-separate pathways and sequences as separate columns.
pathabun_in_strat <- data.frame(t(pathabun_in_strat), check.names = FALSE)
orig_col <- colnames(pathabun_in_strat)
pathabun_in_strat$pathway <- as.character(sapply(rownames(pathabun_in_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
pathabun_in_strat$sequence <- as.character(sapply(rownames(pathabun_in_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
pathabun_in_strat <- pathabun_in_strat[, c("pathway", "sequence", orig_col)]

# Unstrat pathways and phenotypes.
pathabun_in_unstrat <- readRDS("hmp2_pathabun_filt_count_Ileum.rds")
pathabun_in_unstrat <- data.frame(t(pathabun_in_unstrat), check.names = FALSE)
pheno_in_unstrat <- readRDS("hmp2_pheno_count_Ileum.rds")

# Read in ASVs:
asv_abun <- readRDS("hmp2_biom_count_Ileum.rds")

descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info_prokaryotes.txt.gz', 'rt')

path_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)
close(descrip_gzfile)

# Read in ASV taxonomy.
asv_tax_in <- read.table("../hmp2_asv_taxa.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)
asv_tax_in_levels <- add_tax_cols(asv_tax_in)
rownames(asv_tax_in_levels) <- asv_tax_in_levels$Feature.ID

# Get metadata subsetted to 16S samples.
biopsy_16S_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "biopsy_16S"), ]
biopsy_16S_meta <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Ileum"), ]
biopsy_16S_meta <- biopsy_16S_meta[-which(duplicated(biopsy_16S_meta$Participant.ID)), ]
rownames(biopsy_16S_meta) <- biopsy_16S_meta$Participant.ID

Ileum_CD_samples <- biopsy_16S_meta[which(biopsy_16S_meta$diagnosis == "CD"), "Participant.ID"]
Ileum_nonIBD_samples <- biopsy_16S_meta[which(biopsy_16S_meta$diagnosis == "nonIBD"), "Participant.ID"]
Ileum_CD_nonIBD_samples <- c(Ileum_CD_samples, Ileum_nonIBD_samples)


Ileum_CD_non_IBD_pathabun_in <- pathabun_in_unstrat[, Ileum_CD_nonIBD_samples[which(Ileum_CD_nonIBD_samples %in% colnames(pathabun_in_unstrat))]]
Ileum_CD_non_IBD_asv_abun <- asv_abun[, Ileum_CD_nonIBD_samples[which(Ileum_CD_nonIBD_samples %in% colnames(pathabun_in_unstrat))]]
Ileum_CD_non_IBD_asv_abun_taxa <- abun_by_taxa_level(asv_tax_in_levels, Ileum_CD_non_IBD_asv_abun)


Ileum_CD_non_IBD_pathabun_in_filt_relab <- data.frame(sweep(Ileum_CD_non_IBD_pathabun_in, 2, colSums(Ileum_CD_non_IBD_pathabun_in), FUN="/")) * 100
Ileum_CD_non_IBD_asv_abun_relab <- data.frame(sweep(Ileum_CD_non_IBD_asv_abun, 2, colSums(Ileum_CD_non_IBD_asv_abun), FUN="/")) * 100
Ileum_CD_non_IBD_asv_abun_Species_relab <- data.frame(sweep(Ileum_CD_non_IBD_asv_abun_taxa$Species, 2, colSums(Ileum_CD_non_IBD_asv_abun_taxa$Species), FUN="/")) * 100
Ileum_CD_non_IBD_asv_abun_Genus_relab <- data.frame(sweep(Ileum_CD_non_IBD_asv_abun_taxa$Genus, 2, colSums(Ileum_CD_non_IBD_asv_abun_taxa$Genus), FUN="/")) * 100
Ileum_CD_non_IBD_asv_abun_Phylum_relab <- data.frame(sweep(Ileum_CD_non_IBD_asv_abun_taxa$Phylum, 2, colSums(Ileum_CD_non_IBD_asv_abun_taxa$Phylum), FUN="/")) * 100

ASV_samples_meta <- biopsy_16S_meta[colnames(Ileum_CD_non_IBD_asv_abun), ]
ASV_CD_samples <- ASV_samples_meta[which(ASV_samples_meta$diagnosis == "CD"), "Participant.ID"]
ASV_nonIBD_samples <- ASV_samples_meta[which(ASV_samples_meta$diagnosis == "nonIBD"), "Participant.ID"]

sig_taxa_df <- data.frame(asv_sig=as.numeric(Ileum_CD_non_IBD_asv_abun_relab["2031d34eae50b711bfb7c1a7b9a22f39",]),
                          s_sig=as.numeric(Ileum_CD_non_IBD_asv_abun_Species_relab["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__",]),
                          g_sig=as.numeric(Ileum_CD_non_IBD_asv_abun_Genus_relab["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira",]),
                          p_sig=as.numeric(Ileum_CD_non_IBD_asv_abun_Phylum_relab["k__Bacteria; p__Proteobacteria",]),
                          diagnosis=biopsy_16S_meta[colnames(Ileum_CD_non_IBD_asv_abun_relab), "diagnosis"])


asv_boxplot <- ggplot(sig_taxa_df, aes(x=diagnosis, y=asv_sig, fill=diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("ASV for unclassified o__Clostridiales") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

species_boxplot <- ggplot(sig_taxa_df, aes(x=diagnosis, y=s_sig, fill=diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("g__Lachnospira; s__") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

genus_boxplot <- ggplot(sig_taxa_df, aes(x=diagnosis, y=g_sig, fill=diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("g__Lachnospira") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 


phylum_boxplot <- ggplot(sig_taxa_df, aes(x=diagnosis, y=p_sig, fill=diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size=1, dodge.width=0.8) +
  ylab("Relative Abundance (%)") +
  xlab("") +
  ggtitle("p__Proteobacteria") +
  scale_fill_manual(values=c("white", "grey")) +
  theme(legend.position = "none") 

# 10 x 8
plot_grid(asv_boxplot, species_boxplot, genus_boxplot, phylum_boxplot, labels=c("A", "B", "C", "D"))
