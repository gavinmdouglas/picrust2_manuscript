### Commands to prep all files for HMP2 analyses.
### i.e. transform, filter, and subset tables and write out to be used for subsequent analyses.

rm(list=ls(all=TRUE))

setwd("/home/gavin/projects/hmp2_ibd_working/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")

# Read in hmp2 metadata table.
hmp2_metadata <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_2018-08-20_col_subset.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE, check.names = TRUE)

### Prep PICRUSt2 pathway predictions ###
# Read in stratified per-seq-contrib pathway abundances.
hmp2_pathabun_strat <- read.table("/home/gavin/projects/hmp2_ibd_working/16S/April2018_redownload/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                  header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Get unstratified pathway abundances based on the per-seq-contrib table (so that contributors will sum to unstratified abundance).
hmp2_pathabun_strat_tmp <- hmp2_pathabun_strat[, -which(colnames(hmp2_pathabun_strat) == "sequence")]
hmp2_pathabun <- aggregate(. ~ pathway, data=hmp2_pathabun_strat_tmp, FUN=sum)
rm(hmp2_pathabun_strat_tmp)
rownames(hmp2_pathabun) <- hmp2_pathabun$pathway
hmp2_pathabun <- hmp2_pathabun[, -which(colnames(hmp2_pathabun) == "pathway")]

colnames(hmp2_pathabun) <- gsub("HMP\\.2\\.", "", colnames(hmp2_pathabun))

# Remove any pathways classified as "engineered" or "superpathways".
descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info_prokaryotes.txt.gz', 'rt')

pathway_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)

close(descrip_gzfile)

pathway_descrip_subset <- pathway_descrip[rownames(hmp2_pathabun),, drop=FALSE]

path2remove_i <- c(grep("superpathway", pathway_descrip_subset$V2), grep("engineered", pathway_descrip_subset$V2))

path2remove <- rownames(pathway_descrip_subset)[path2remove_i]

hmp2_pathabun_filt <- hmp2_pathabun[-which(rownames(hmp2_pathabun) %in% path2remove), ]

hmp2_pathabun_filt_count <- hmp2_pathabun_filt

# Convert to relative abundance and perform asin transformation.
hmp2_pathabun_filt <- data.frame(sweep(hmp2_pathabun_filt, 2, colSums(hmp2_pathabun_filt), '/'), check.names=FALSE)
hmp2_pathabun_filt_asin <- asinTransform(hmp2_pathabun_filt)

# Split into ileum and rectum samples and re-name columns to be participant ids.
# Sort mapfiles so that 1st week should be first (in case of duplicates from same patient) and the deduplicate.
# Also transpose final abundance tables.
biopsy_16S_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "biopsy_16S"), ]

biopsy_16S_meta_Ileum <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Ileum"), ]
biopsy_16S_meta_Ileum <- biopsy_16S_meta_Ileum[with(biopsy_16S_meta_Ileum, order(Participant.ID, week_num)),]
biopsy_16S_meta_Ileum <- biopsy_16S_meta_Ileum[-which(duplicated(biopsy_16S_meta_Ileum$Participant.ID)), ]
biopsy_16S_meta_Ileum <- biopsy_16S_meta_Ileum[which(biopsy_16S_meta_Ileum$External.ID %in% colnames(hmp2_pathabun_filt_asin)), ]

hmp2_pathabun_filt_asin_Ileum <- data.frame(t(hmp2_pathabun_filt_asin[, biopsy_16S_meta_Ileum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_filt_asin_Ileum) <- biopsy_16S_meta_Ileum$Participant.

biopsy_16S_meta_Rectum <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Rectum"), ]
biopsy_16S_meta_Rectum <- biopsy_16S_meta_Rectum[with(biopsy_16S_meta_Rectum, order(Participant.ID, week_num)),]
biopsy_16S_meta_Rectum <- biopsy_16S_meta_Rectum[-which(duplicated(biopsy_16S_meta_Rectum$Participant.ID)), ]
biopsy_16S_meta_Rectum <- biopsy_16S_meta_Rectum[which(biopsy_16S_meta_Rectum$External.ID %in% colnames(hmp2_pathabun_filt_asin)), ]

hmp2_pathabun_filt_asin_Rectum <- data.frame(t(hmp2_pathabun_filt_asin[, biopsy_16S_meta_Rectum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_filt_asin_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID

##### Also prep stratified pathway predictions. ##### 
# Note that this data wont be transformed since we're interested in the relative abundance of contributors.

colnames(hmp2_pathabun_strat) <- gsub("HMP\\.2\\.", "", colnames(hmp2_pathabun_strat))

hmp2_pathabun_strat_filt <- hmp2_pathabun_strat[-which(hmp2_pathabun_strat$pathway %in% path2remove), ]
rownames(hmp2_pathabun_strat_filt) <- paste(hmp2_pathabun_strat_filt$pathway, hmp2_pathabun_strat_filt$sequence, sep="|")
hmp2_pathabun_strat_filt <- hmp2_pathabun_strat_filt[, -which(colnames(hmp2_pathabun_strat_filt) %in% c("pathway", "sequence"))]

hmp2_pathabun_strat_filt_count <- hmp2_pathabun_strat_filt

# Convert to relative abundance
hmp2_pathabun_strat_filt <- data.frame(sweep(hmp2_pathabun_strat_filt, 2, colSums(hmp2_pathabun_strat_filt), '/'), check.names=FALSE)

# Subset to Ileum and Rectum samples.
hmp2_pathabun_strat_filt_Ileum <- data.frame(t(hmp2_pathabun_strat_filt[, biopsy_16S_meta_Ileum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_strat_filt_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_pathabun_strat_filt_Rectum <- data.frame(t(hmp2_pathabun_strat_filt[, biopsy_16S_meta_Rectum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_strat_filt_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID



####### Prep ASV relative abundances ####### 
hmp2_biom <- read.table("16S/April2018_redownload/deblur_output_final/hmp2_ibd_16S_4000reads.biom.tsv",
                        header=TRUE, skip=1, row.names=1, comment.char="", sep="\t")

colnames(hmp2_biom) <- gsub("HMP\\.2\\.", "", colnames(hmp2_biom))

hmp2_biom_count <- hmp2_biom

hmp2_biom <- data.frame(sweep(hmp2_biom, 2, colSums(hmp2_biom), '/'), check.names = FALSE)

# Split into Ileum and Rectum samples based on the above mapfiles used for processing pathway abundance.
hmp2_biom_Ileum <- hmp2_biom[, biopsy_16S_meta_Ileum$External.ID]
colnames(hmp2_biom_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_biom_Rectum <- hmp2_biom[, biopsy_16S_meta_Rectum$External.ID]
colnames(hmp2_biom_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID

# Also read in taxa classifications for this subset of ASVs.
hmp2_asv_taxa <- read.table("16S/April2018_redownload/deblur_output_final/hmp2_ibd_16S_taxonomy_unfilt.tsv",
                             header=TRUE, sep="\t", stringsAsFactors = FALSE)

hmp2_asv_taxa_subset <- hmp2_asv_taxa[which(hmp2_asv_taxa$Feature.ID %in% rownames(hmp2_biom)), ]

####### Prep IMG phenotype predictions ####### 
hmp2_pheno <- read.table("16S/April2018_redownload/picrust2_full_output_2.1.0-b/metagenome_out_pheno/pred_metagenome_unstrat.tsv",
                        header=TRUE,row.names=1, sep="\t")

colnames(hmp2_pheno) <- gsub("HMP\\.2\\.", "", colnames(hmp2_pheno))

hmp2_pheno_count <- hmp2_pheno

hmp2_pheno <- data.frame(sweep(hmp2_pheno, 2, colSums(hmp2_pheno), '/'), check.names = FALSE)

hmp2_pheno_asin <- asinTransform(hmp2_pheno)

# Split into Ileum and Rectum samples based on the above mapfiles used for processing pathway abundance.
hmp2_pheno_asin_Ileum <- hmp2_pheno_asin[, biopsy_16S_meta_Ileum$External.ID]
colnames(hmp2_pheno_asin_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_pheno_asin_Rectum <- hmp2_pheno_asin[, biopsy_16S_meta_Rectum$External.ID]
colnames(hmp2_pheno_asin_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID




##### Also get count tables formatted for pathways, ASVs, and phenos.
hmp2_pathabun_filt_count_Ileum <- data.frame(t(hmp2_pathabun_filt_count[, biopsy_16S_meta_Ileum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_filt_count_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_pathabun_strat_filt_count_Ileum <- data.frame(t(hmp2_pathabun_strat_filt_count[, biopsy_16S_meta_Ileum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_strat_filt_count_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_biom_count_Ileum <- hmp2_biom_count[, biopsy_16S_meta_Ileum$External.ID]
colnames(hmp2_biom_count_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID

hmp2_pheno_count_Ileum <- hmp2_pheno_count[, biopsy_16S_meta_Ileum$External.ID]
colnames(hmp2_pheno_count_Ileum) <- biopsy_16S_meta_Ileum$Participant.ID


hmp2_pathabun_filt_count_Rectum <- data.frame(t(hmp2_pathabun_filt_count[, biopsy_16S_meta_Rectum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_filt_count_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID

hmp2_pathabun_strat_filt_count_Rectum <- data.frame(t(hmp2_pathabun_strat_filt_count[, biopsy_16S_meta_Rectum$External.ID]), check.names=FALSE)
rownames(hmp2_pathabun_strat_filt_count_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID

hmp2_biom_count_Rectum <- hmp2_biom_count[, biopsy_16S_meta_Rectum$External.ID]
colnames(hmp2_biom_count_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID

hmp2_pheno_count_Rectum <- hmp2_pheno_count[, biopsy_16S_meta_Rectum$External.ID]
colnames(hmp2_pheno_count_Rectum) <- biopsy_16S_meta_Rectum$Participant.ID



####### Prep metabolome data #######
# Load metabolome data.
hmp2_metabolome <- read.csv("HMP2_metabolomics_proteomics/HMP2_metabolomics.csv",
                            header=T, sep=",", comment.char="", stringsAsFactors = FALSE)

# Write out map of compounds to metabolites.
compound_metabolite_map <- hmp2_metabolome[, c("Compound", "Metabolite")]

# Remove any compounds that are of unidentified metabolite.
rownames(hmp2_metabolome) <- hmp2_metabolome$Compound
hmp2_metabolome_filt <- hmp2_metabolome[-which(hmp2_metabolome$Metabolite == ""), ]

# Remove first 7 columns, which contain non-abundance data.
hmp2_metabolome_filt <- hmp2_metabolome_filt[, -c(1:7)]

# Set NA metabolite abundances to be 0.
hmp2_metabolome_filt[is.na(hmp2_metabolome_filt)] <- 0

# Subset metadata to metabolomics rows.
metabolome_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "metabolomics"), ]

# Sort by participant and week number so that first sample should always be first per participant.
metabolome_meta_sorted <- metabolome_meta[with(metabolome_meta, order(Participant.ID, week_num)),]

# Remove all samples except for 1st per patient.
metabolome_meta_sorted_dedup <- metabolome_meta_sorted[-which(duplicated(metabolome_meta_sorted$Participant.ID)), ]

# Subset the metabolome table to this subset of samples and rename to be participant ids.
hmp2_metabolome_filt_subset <- hmp2_metabolome_filt[, metabolome_meta_sorted_dedup$External.ID]
colnames(hmp2_metabolome_filt_subset) <- metabolome_meta_sorted_dedup$Participant.ID

# log transform after adding pseudocount of 1 and transpose.
hmp2_metabolome_filt_subset_log <- data.frame(t(log(hmp2_metabolome_filt_subset + 1)), check.names=FALSE)


#### Prep RNA-seq data ####
# Read in RNA table with genes as columns and samples as rows.
host_rna <- data.frame(t(read.table("host_rnaseq/host_tx_counts.biom.tsv",
                                    skip=1, comment="", sep="\t", header=T, row.names=1)))

# Identify samples that have fewer than 20,000,000 raw reads and remove them from rnaseq.
# Resulted in 5 samples being removed.
host_rna <- host_rna[-which(rowSums(host_rna) < 20e6), ]

# Normalize RNAseq table to be reads per million.
host_rna_rpm <- host_rna / (rowSums(host_rna) / 1e6)

# Log transform after adding pseudocount of 1.
host_rna_rpm_log <- log(host_rna_rpm + 1)

# Get RNA-seq subset of mapfile and subset RPM table to Ileum and Rectum biopsies
# and then re-name rows to be participant ids.
rnaseq_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "host_transcriptomics"), ]

rnaseq_meta_Ileum <- rnaseq_meta[which(rnaseq_meta$biopsy_location == "Ileum"), ]
rnaseq_meta_Ileum <- rnaseq_meta_Ileum[with(rnaseq_meta_Ileum, order(Participant.ID, week_num)),]
rnaseq_meta_Ileum <- rnaseq_meta_Ileum[-which(duplicated(rnaseq_meta_Ileum$Participant.ID)), ]
rnaseq_meta_Ileum <- rnaseq_meta_Ileum[which(rnaseq_meta_Ileum$External.ID %in% rownames(host_rna_rpm_log)), ]

host_rna_rpm_log_Ileum <- host_rna_rpm_log[rnaseq_meta_Ileum$External.ID, ]
rownames(host_rna_rpm_log_Ileum) <- rnaseq_meta_Ileum$Participant.ID

rnaseq_meta_Rectum <- rnaseq_meta[which(rnaseq_meta$biopsy_location == "Rectum"), ]
rnaseq_meta_Rectum <- rnaseq_meta_Rectum[with(rnaseq_meta_Rectum, order(Participant.ID, week_num)),]
rnaseq_meta_Rectum <- rnaseq_meta_Rectum[-which(duplicated(rnaseq_meta_Rectum$Participant.ID)), ]
rnaseq_meta_Rectum <- rnaseq_meta_Rectum[which(rnaseq_meta_Rectum$External.ID %in% rownames(host_rna_rpm_log)), ]

host_rna_rpm_log_Rectum <- host_rna_rpm_log[rnaseq_meta_Rectum$External.ID, ]
rownames(host_rna_rpm_log_Rectum) <- rnaseq_meta_Rectum$Participant.ID


####### Prep MGS - HUMAnN2 output table #######
hmp2_mgs <- read.table("mgs/pathabundance_relab.pcl.tsv",
                       header=TRUE, sep="\t", comment.char="", quote="", stringsAsFactors = FALSE, row.names=1)

# First 14 rows contain metadata, which can be removed.
hmp2_mgs <- hmp2_mgs[-c(1:14), ]

# Get sample names to keep and re-name columns as above.
mgs_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "metagenomics"), ]
mgs_meta_sorted <- mgs_meta[with(mgs_meta, order(Participant.ID, week_num)),]
mgs_meta_sorted_dedup <- mgs_meta_sorted[-which(duplicated(mgs_meta_sorted$Participant.ID)), ]
hmp2_mgs_subset <- hmp2_mgs[, mgs_meta_sorted_dedup$External.ID]
colnames(hmp2_mgs_subset) <- mgs_meta_sorted_dedup$Participant.ID
orig_rownames <- rownames(hmp2_mgs_subset)

# Convert columns to numeric.
hmp2_mgs_subset <- sapply(hmp2_mgs_subset, as.numeric)
rownames(hmp2_mgs_subset) <- orig_rownames

# Remove description from rownames:
rownames(hmp2_mgs_subset) <- gsub(":.*\\|", "|", rownames(hmp2_mgs_subset))
rownames(hmp2_mgs_subset) <- gsub(":.*$", "", rownames(hmp2_mgs_subset))

# Split into strat and unstratified and convert to relative abundance.
strat_rows <- grep("\\|", rownames(hmp2_mgs_subset))

hmp2_mgs_subset_strat <- hmp2_mgs_subset[strat_rows, ]
hmp2_mgs_subset_unstrat <- hmp2_mgs_subset[-strat_rows, ]

hmp2_mgs_subset_strat_relab <- data.frame(sweep(hmp2_mgs_subset_strat, 2, colSums(hmp2_mgs_subset_strat), '/'))
hmp2_mgs_subset_unstrat_relab <- data.frame(sweep(hmp2_mgs_subset_unstrat, 2, colSums(hmp2_mgs_subset_unstrat), '/'))

# For unstrat table, also perform asin transformation.
hmp2_mgs_subset_unstrat_relab_asin <- asinTransform(hmp2_mgs_subset_unstrat_relab)


####### Print out final tables for samples of interest. #######

### First, print out misc mapping tables:
# ASV to taxa classifications.
write.table(x=hmp2_asv_taxa_subset,
            file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/hmp2_asv_taxa.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

# Compound to metabolite mapping.
write.table(x=compound_metabolite_map,
            file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/compound_metabolite_map.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

# Write out full tables.
saveRDS(object = hmp2_pathabun_strat_filt_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_pathabun_strat_filt_Ileum.rds")
saveRDS(object = hmp2_pathabun_filt_asin_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_pathabun_filt_asin_Ileum.rds")
saveRDS(object = hmp2_biom_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_biom_Ileum.rds")
saveRDS(object = hmp2_pheno_asin_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_pheno_asin_Ileum.rds")
saveRDS(object = host_rna_rpm_log_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_rna_rpm_log_Ileum.rds")

saveRDS(object = hmp2_pathabun_strat_filt_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_pathabun_strat_filt_Rectum.rds")
saveRDS(object = hmp2_pathabun_filt_asin_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_pathabun_filt_asin_Rectum.rds")
saveRDS(object = hmp2_biom_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_biom_Rectum.rds")
saveRDS(object = hmp2_pheno_asin_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_pheno_asin_Rectum.rds")
saveRDS(object = host_rna_rpm_log_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_rna_rpm_log_Rectum.rds")


saveRDS(object = hmp2_metabolome_filt_subset_log,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_metabolome_filt_subset_log.rds")
saveRDS(object = hmp2_mgs_subset_strat_relab,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_mgs_subset_strat_relab.rds")
saveRDS(object = hmp2_mgs_subset_unstrat_relab_asin,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/prepped_tables/hmp2_mgs_subset_unstrat_relab_asin.rds")



# Also save tables with raw counts.
saveRDS(object = hmp2_pathabun_filt_count_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_pathabun_filt_count_Ileum.rds")

saveRDS(object = hmp2_pathabun_strat_filt_count_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_pathabun_strat_filt_count_Ileum.rds")

saveRDS(object = hmp2_biom_count_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_biom_count_Ileum.rds")

saveRDS(object = hmp2_pheno_count_Ileum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_pheno_count_Ileum.rds")

saveRDS(object = hmp2_pathabun_filt_count_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_pathabun_filt_count_Rectum.rds")

saveRDS(object = hmp2_pathabun_strat_filt_count_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_pathabun_strat_filt_count_Rectum.rds")

saveRDS(object = hmp2_biom_count_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_biom_count_Rectum.rds")

saveRDS(object = hmp2_pheno_count_Rectum,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/hmp2_pheno_count_Rectum.rds")
