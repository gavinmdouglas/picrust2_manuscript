### Compare taxonomic contributions to pathways between the 16S (biopsy) and MGS (stool) datarm(list=ls(all=TRUE)).

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

# Read in stratified HMP2 pathway abundances
hmp2_16S_pathabun_strat <- data.frame(t(readRDS("count_tables/hmp2_pathabun_strat_filt_count_Ileum.rds")), check.names=FALSE)

hmp2_mgs_pathabun_strat <- readRDS("prepped_tables/hmp2_mgs_subset_strat_relab.rds") * 100

# Subset to overlapping CD and nonIBD samples only.
hmp2_metadata <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_2018-08-20_col_subset.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE)
biopsy_16S_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "biopsy_16S"), ]
biopsy_16S_meta <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Ileum"), ]
biopsy_16S_meta <- biopsy_16S_meta[-which(duplicated(biopsy_16S_meta$Participant.ID)), ]
rownames(biopsy_16S_meta) <- biopsy_16S_meta$Participant.ID

raw_CD_nonIBD_Ileum_samples <- c(biopsy_16S_meta[which(biopsy_16S_meta$diagnosis == "CD"), "Participant.ID"],
                                 biopsy_16S_meta[which(biopsy_16S_meta$diagnosis == "nonIBD"), "Participant.ID"])

CD_nonIBD_Ileum_samples <- raw_CD_nonIBD_Ileum_samples[which(raw_CD_nonIBD_Ileum_samples %in% colnames(hmp2_16S_pathabun_strat))]
CD_nonIBD_Ileum_samples <- CD_nonIBD_Ileum_samples[which(CD_nonIBD_Ileum_samples %in% colnames(hmp2_mgs_pathabun_strat))]

hmp2_16S_pathabun_strat <- hmp2_16S_pathabun_strat[, CD_nonIBD_Ileum_samples]
hmp2_16S_pathabun_strat <- hmp2_16S_pathabun_strat[-which(rowSums(hmp2_16S_pathabun_strat) == 0), ]

hmp2_mgs_pathabun_strat <- hmp2_mgs_pathabun_strat[, CD_nonIBD_Ileum_samples]
hmp2_mgs_pathabun_strat <- hmp2_mgs_pathabun_strat[-which(rowSums(hmp2_mgs_pathabun_strat) == 0), ]

# Rearrange both tables to be contributions by genera.
asv_tax_in <- read.table("hmp2_asv_taxa.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
asv_tax_in_levels <- add_tax_cols(asv_tax_in)
rownames(asv_tax_in_levels) <- asv_tax_in_levels$Feature.ID

mgs_orig_col <- colnames(hmp2_16S_pathabun_strat$sequence)
picrust2_orig_col <- colnames(hmp2_16S_pathabun_strat)

hmp2_mgs_pathabun_strat$pathway <- as.character(sapply(rownames(hmp2_mgs_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
hmp2_mgs_pathabun_strat$genus <- as.character(sapply(rownames(hmp2_mgs_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
hmp2_mgs_pathabun_strat$genus <- gsub(".s__.*$", "", hmp2_mgs_pathabun_strat$genus)

hmp2_16S_pathabun_strat$pathway <- as.character(sapply(rownames(hmp2_16S_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
hmp2_16S_pathabun_strat$sequence <- as.character(sapply(rownames(hmp2_16S_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
hmp2_16S_pathabun_strat$genus <- gsub("^.*g__", "g__", asv_tax_in_levels[hmp2_16S_pathabun_strat$sequence, "Genus"])
hmp2_16S_pathabun_strat$genus[which(hmp2_16S_pathabun_strat$genus == "g__")] <- "unclassified"
hmp2_16S_pathabun_strat <- hmp2_16S_pathabun_strat[, -which(colnames(hmp2_16S_pathabun_strat) == "sequence")]

# Only keep pathways in PICRUSt2 mapfile (i.e. prokaryotic pathways only) that aren't superpathways or engineered.
descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 'rt')

pathway_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)

close(descrip_gzfile)

path2keep <- hmp2_mgs_pathabun_strat$pathway

if(length(which(duplicated(path2keep))) > 0) {
  path2keep <- path2keep[-which(duplicated(path2keep))]
}

pathway_descrip_subset <- pathway_descrip[path2keep, , drop=FALSE]




path2remove_i <- c(grep("superpathway", pathway_descrip_subset$V2), grep("engineered", pathway_descrip_subset$V2))

path2remove <- rownames(pathway_descrip_subset)[path2remove_i]

hmp2_mgs_pathabun_strat <- hmp2_mgs_pathabun_strat[-which(hmp2_mgs_pathabun_strat$pathway %in% path2remove), ]

# Re-convert both to relative abundance.
hmp2_16S_pathabun_strat[, 1:(ncol(hmp2_16S_pathabun_strat) - 2)] <- data.frame(sweep(hmp2_16S_pathabun_strat[, 1:(ncol(hmp2_16S_pathabun_strat) - 2)],
                                                                                     2, colSums(hmp2_16S_pathabun_strat[, 1:(ncol(hmp2_16S_pathabun_strat) - 2)]), '/'), check.names = FALSE) * 100
hmp2_mgs_pathabun_strat[, 1:(ncol(hmp2_mgs_pathabun_strat) - 2)] <- data.frame(sweep(hmp2_mgs_pathabun_strat[, 1:(ncol(hmp2_mgs_pathabun_strat) - 2)],
                                                                                     2, colSums(hmp2_mgs_pathabun_strat[, 1:(ncol(hmp2_mgs_pathabun_strat) - 2)]), '/'), check.names = FALSE) * 100

hmp2_mgs_pathabun_strat_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=hmp2_mgs_pathabun_strat)
hmp2_16S_pathabun_strat_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=hmp2_16S_pathabun_strat)

##################

# Also get collapsed genera without requiring that final genus level be classified and save for other analyses.
hmp2_16S_pathabun_strat_full <- data.frame(t(readRDS("count_tables/hmp2_pathabun_strat_filt_count_Ileum.rds")), check.names=FALSE)
hmp2_16S_pathabun_strat_full <- hmp2_16S_pathabun_strat_full[, CD_nonIBD_Ileum_samples]
hmp2_16S_pathabun_strat_full <- hmp2_16S_pathabun_strat_full[-which(rowSums(hmp2_16S_pathabun_strat_full) == 0), ]

hmp2_16S_pathabun_strat_full$pathway <- as.character(sapply(rownames(hmp2_16S_pathabun_strat_full), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
hmp2_16S_pathabun_strat_full$sequence <- as.character(sapply(rownames(hmp2_16S_pathabun_strat_full), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
hmp2_16S_pathabun_strat_full$genus <- asv_tax_in_levels[hmp2_16S_pathabun_strat_full$sequence, "Genus"]
hmp2_16S_pathabun_strat_full <- hmp2_16S_pathabun_strat_full[, -which(colnames(hmp2_16S_pathabun_strat_full) == "sequence")]
hmp2_16S_pathabun_strat_full[, 1:(ncol(hmp2_16S_pathabun_strat_full) - 2)] <- data.frame(sweep(hmp2_16S_pathabun_strat_full[, 1:(ncol(hmp2_16S_pathabun_strat_full) - 2)],
                                                 2, colSums(hmp2_16S_pathabun_strat_full[, 1:(ncol(hmp2_16S_pathabun_strat_full) - 2)]), '/'), check.names = FALSE) * 100
hmp2_16S_pathabun_strat_full_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=hmp2_16S_pathabun_strat_full)
#saveRDS(object = hmp2_16S_pathabun_strat_full_genus_sum, file = "results_out/hmp2_16S_pathabun_strat_full_genus_sum.rds")

###############

num_contrib_genera_out <- num_contrib_genera(pred_df = hmp2_16S_pathabun_strat_genus_sum,
                                             mgs_df = hmp2_mgs_pathabun_strat_genus_sum)

num_contrib_genera_out[is.na(num_contrib_genera_out)] <- 0

# Save output file.
saveRDS(object = num_contrib_genera_out, file = "results_out/num_contrib_genera_out.rds")


unique_path <- c(hmp2_16S_pathabun_strat_genus_sum$pathway, hmp2_mgs_pathabun_strat_genus_sum$pathway)
if(length(which(duplicated(unique_path))) > 0) {
  unique_path <- unique_path[-which(duplicated(unique_path))]
}

path_spearman <- rep(NA, length(unique_path))
names(path_spearman) <- unique_path

for(pathway in unique_path) {
  pathway_breakdown <- breakdown_mean_genera_contrib(hmp2_16S_pathabun_strat_genus_sum, hmp2_mgs_pathabun_strat_genus_sum, pathway)

  if(nrow(pathway_breakdown) <= 1) {
    next 
  }

  path_spearman[pathway] <- cor.test(pathway_breakdown$mgs_mean, pathway_breakdown$picrust2_mean, method="spearman")$estimate
}

saveRDS(object = path_spearman, file = "results_out/16S_vs_MGS_contrib_pathway_spearman.rds")

### TESTING:
# Next compare contributions to 1 pathway expected to be contributed mainly by Proteobacteria (PWY-5188) vs a different one.
# hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1 <- readRDS("results_out/hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.rds")
# 
# proteobacteria_enriched_pathways <- levels(cd_sig_higher_ratio_prep_melt$variable)
# 
# hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1[which(hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1$gene == "NAT8"), ]
# # PWY0-1533 NAT8 0.004757666 0.09840439
# 
# 
# PWY_5188_breakdown <- breakdown_mean_genera_contrib(hmp2_16S_pathabun_strat_genus_sum, hmp2_mgs_pathabun_strat_genus_sum, "PWY-5188")
# PWY0_1533_breakdown <- breakdown_mean_genera_contrib(hmp2_16S_pathabun_strat_genus_sum, hmp2_mgs_pathabun_strat_genus_sum, "PWY0-1533")
# 
# plot(PWY_5188_breakdown$mgs_mean, PWY_5188_breakdown$picrust2_mean)
# plot(PWY0_1533_breakdown$mgs_mean, PWY0_1533_breakdown$picrust2_mean)
# 
# ### PWY0-1533 is in MGS and at a similar relative abundance as 16S data.
# 
# 
# ### Stacked barchart of genera contributions.
# hmp2_16S_pathabun_strat_genus_sum_PWY_5188 <- hmp2_16S_pathabun_strat_genus_sum[which(hmp2_16S_pathabun_strat_genus_sum$pathway == "PWY0-1533"), ]
# hmp2_16S_pathabun_strat_genus_sum_PWY_5188_melt <- melt(hmp2_16S_pathabun_strat_genus_sum_PWY_5188)
# 
# ggplot(hmp2_16S_pathabun_strat_genus_sum_PWY_5188_melt, aes(x=variable, y=value, fill=genus)) +
#   geom_bar(stat="identity")
# 
