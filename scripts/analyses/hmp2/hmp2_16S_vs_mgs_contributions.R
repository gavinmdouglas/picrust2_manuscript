### Compare taxonomic contributions to pathways between the 16S (biopsy) and MGS (stool) datarm(list=ls(all=TRUE)).

rm(list=ls(all=TRUE))

# First compare numbers of genera contributing pathways (excluding "unclassified" in 16S biopsies and mgs stool).
num_contrib_genera <- function(pred_df, mgs_df) {
  
  # Remove unclassified rows.
  pred_df <- pred_df[-which(pred_df$genus == "unclassified"), ]
  mgs_df <- mgs_df[-which(mgs_df$genus == "unclassified"), ]
  
  all_path <- unique(c(pred_df$pathway, mgs_df$pathway))
  
  num_contrib_df <- data.frame(matrix(NA, nrow=length(all_path), ncol=4))
  colnames(num_contrib_df) <- c("picrust2_mean", "picrust2_sem", "mgs_mean", "mgs_sem")
  rownames(num_contrib_df) <- all_path
  
  for(pathway in all_path) {
    
    if(pathway %in% pred_df$pathway) {
      pred_df_subset <- pred_df[which(pred_df$pathway == pathway), ]
      
      num_contrib_df[pathway, "picrust2_mean"] <- mean(colSums(pred_df_subset > 0))
      num_contrib_df[pathway, "picrust2_sem"] <- sd(colSums(pred_df_subset > 0)) / sqrt(ncol(pred_df_subset))
    }
    
    if(pathway %in% mgs_df$pathway) {
      mgs_df_subset <- mgs_df[which(mgs_df$pathway == pathway), ]
      
      num_contrib_df[pathway, "mgs_mean"] <- mean(colSums(mgs_df_subset > 0))
      num_contrib_df[pathway, "mgs_sem"] <- sd(colSums(mgs_df_subset > 0)) / sqrt(ncol(mgs_df_subset))
    }
    
  }
  
  return(num_contrib_df)
}

pathway_mean_sem <- function(in_df, pathway, col_str) {
  
  pathway_subset_abun <- in_df[which(in_df$pathway == pathway), ]
  
  pathway_subset_abun <- pathway_subset_abun[-which(pathway_subset_abun$genus == "unclassified"), ]
  
  rownames(pathway_subset_abun) <- pathway_subset_abun$genus
  
  pathway_subset_abun <- pathway_subset_abun[, -which(colnames(pathway_subset_abun) %in% c("genus", "pathway"))]
  
  pathway_subset_abun <- data.frame(sweep(pathway_subset_abun, 2, colSums(pathway_subset_abun), '/'), check.names = FALSE) * 100
  
  if(length(which(is.na(pathway_subset_abun))) > 0) {
    pathway_subset_abun[is.na(pathway_subset_abun)] <- 0
  }
  pathway_subset_abun_genus_mean <- rowMeans(pathway_subset_abun)
  
  pathway_subset_abun_genus_sd <- apply(pathway_subset_abun, 1, sd)
  
  pathway_subset_abun_genus_sem <- pathway_subset_abun_genus_sd / sqrt(ncol(pathway_subset_abun))
  
  pathway_subset_abun_genus_upper <- pathway_subset_abun_genus_mean + pathway_subset_abun_genus_sem * 2
  pathway_subset_abun_genus_lower <- pathway_subset_abun_genus_mean - pathway_subset_abun_genus_sem * 2
  
  if(length(which(pathway_subset_abun_genus_lower < 0)) > 0) {
    pathway_subset_abun_genus_lower[which(pathway_subset_abun_genus_lower < 0)] <- 0
  }
  
  out_df <- data.frame(mean=pathway_subset_abun_genus_mean,
                       sem=pathway_subset_abun_genus_sem,
                       lower=pathway_subset_abun_genus_lower,
                       upper=pathway_subset_abun_genus_upper)
  rownames(out_df) <- names(pathway_subset_abun_genus_mean) 
  colnames(out_df) <- paste(col_str, colnames(out_df), sep="_")
  
  return(out_df)
}

breakdown_mean_genera_contrib <- function(pred_df, mgs_df, pathway) {
  
  pred_summary <- pathway_mean_sem(pred_df, pathway, "picrust2")
  mgs_summary <- pathway_mean_sem(mgs_df, pathway, "mgs")
  
  merged_out <- merge(pred_summary, mgs_summary, by="row.names", all.x=TRUE)
  
  if(length(which(is.na(merged_out))) > 0 ) {
    merged_out[is.na(merged_out)] <- 0
  }
  return(merged_out)
}


setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")

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
descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info_prokaryotes.txt.gz', 'rt')

pathway_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)

close(descrip_gzfile)

pathway_descrip_subset <- pathway_descrip[unique(hmp2_mgs_pathabun_strat$pathway),, drop=FALSE]

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


# Next compare contributions to 1 pathway expected to be contributed mainly by Proteobacteria (PWY-5188) vs a different one.
hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1 <- readRDS("results_out/hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.rds")

proteobacteria_enriched_pathways <- levels(cd_sig_higher_ratio_prep_melt$variable)

hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1[which(hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1$gene == "NAT8"), ]
# PWY0-1533 NAT8 0.004757666 0.09840439


PWY_5188_breakdown <- breakdown_mean_genera_contrib(hmp2_16S_pathabun_strat_genus_sum, hmp2_mgs_pathabun_strat_genus_sum, "PWY-5188")
PWY0_1533_breakdown <- breakdown_mean_genera_contrib(hmp2_16S_pathabun_strat_genus_sum, hmp2_mgs_pathabun_strat_genus_sum, "PWY0-1533")

plot(PWY_5188_breakdown$mgs_mean, PWY_5188_breakdown$picrust2_mean)
plot(PWY0_1533_breakdown$mgs_mean, PWY0_1533_breakdown$picrust2_mean)

### PWY0-1533 is in MGS and at a similar relative abundance as 16S data.


### Stacked barchart of genera contributions.
hmp2_16S_pathabun_strat_genus_sum_PWY_5188 <- hmp2_16S_pathabun_strat_genus_sum[which(hmp2_16S_pathabun_strat_genus_sum$pathway == "PWY0-1533"), ]
hmp2_16S_pathabun_strat_genus_sum_PWY_5188_melt <- melt(hmp2_16S_pathabun_strat_genus_sum_PWY_5188)

ggplot(hmp2_16S_pathabun_strat_genus_sum_PWY_5188_melt, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity")

