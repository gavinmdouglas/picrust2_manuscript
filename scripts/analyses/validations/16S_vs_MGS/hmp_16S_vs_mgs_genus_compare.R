rm(list=ls(all.names=TRUE))

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 'rt')

pathway_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)

close(descrip_gzfile)

path2remove_i <- c(grep("superpathway", pathway_descrip$V2), grep("engineered", pathway_descrip$V2))

path2remove <- rownames(pathway_descrip)[path2remove_i]


compare_strat_pred_vs_mgs <- function(picrust2_strat, mgs_strat, asv_taxa, mgs_col_str2replace, path2rm, sample_name_map=NULL, blueberry_input=FALSE) {
  
  pathabun_strat <- read.table(picrust2_strat, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  
  pathabun_strat_filt <- pathabun_strat[-which(pathabun_strat$pathway %in% path2rm), ]
  
  # Read in ASV taxonomy.
  tax_in <- read.table(asv_taxa, header=T, sep="\t", stringsAsFactors = FALSE)
  
  tax_in_levels <- add_tax_cols(tax_in)
  rownames(tax_in_levels) <- tax_in_levels$Feature.ID
  
  # Add genus as new column:
  pathabun_strat_filt$sequence <- tax_in_levels[pathabun_strat_filt$sequence, "Genus"]
  colnames(pathabun_strat_filt)[2] <- "genus"
   
  MGS_pathabun_strat <- read.table(mgs_strat, header=TRUE, sep="\t", stringsAsFactors = FALSE, comment.char="", row.names=1, quote="")
  
  colnames(MGS_pathabun_strat) <- gsub(mgs_col_str2replace, "", colnames(MGS_pathabun_strat))
  
  # Remove description from rownames:
  rownames(MGS_pathabun_strat) <- gsub(":.*\\|", "|", rownames(MGS_pathabun_strat))
  rownames(MGS_pathabun_strat) <- gsub(":.*$", "", rownames(MGS_pathabun_strat))
  
  # Subset to stratified rows only.
  strat_rows <- grep("\\|", rownames(MGS_pathabun_strat))
  
  MGS_pathabun_strat <- MGS_pathabun_strat[strat_rows, ]
  
  if(! is.null(sample_name_map)) {
    sample_name_map_in <- read.table(sample_name_map, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    rownames(sample_name_map_in) <- sample_name_map_in$runids_mgs
    colnames(MGS_pathabun_strat) <- sample_name_map_in[colnames(MGS_pathabun_strat), "runids_16S"]
  }
  
  if(blueberry_input) {
    colnames(MGS_pathabun_strat) <- gsub("^BB", "Bact", colnames(MGS_pathabun_strat))
    colnames(MGS_pathabun_strat) <- gsub("_", "-", colnames(MGS_pathabun_strat))
  }

  orig_MGS_col <- colnames(MGS_pathabun_strat)
  
  overlapping_col <- orig_MGS_col[which(orig_MGS_col %in% colnames(pathabun_strat_filt))]
  
  MGS_pathabun_strat$pathway <- as.character(sapply(rownames(MGS_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
  MGS_pathabun_strat$genus <- as.character(sapply(rownames(MGS_pathabun_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
  MGS_pathabun_strat$genus <- gsub(".s__.*$", "", MGS_pathabun_strat$genus)

  pathabun_strat_filt <- pathabun_strat_filt[, c("pathway", "genus", overlapping_col)]
  MGS_pathabun_strat <- MGS_pathabun_strat[, c("pathway", "genus", overlapping_col)]
  
  MGS_pathabun_strat_filt <- MGS_pathabun_strat[-which(MGS_pathabun_strat$pathway %in% path2remove), ]
  
  # Remove NA rows from 16S output (this was a bug in an earlier version of PICRUSt2).
  if(length(which(rowSums(is.na(pathabun_strat_filt)) == ncol(pathabun_strat_filt) - 2)) > 0) {
    pathabun_strat_filt <- pathabun_strat_filt[-which(rowSums(is.na(pathabun_strat_filt)) == ncol(pathabun_strat_filt) - 2), ]
  }
  # Convert both to relative abundance.
  pathabun_strat_filt[, 3:(ncol(pathabun_strat_filt))] <- data.frame(sweep(pathabun_strat_filt[, 3:(ncol(pathabun_strat_filt))],
                                                                                       2, colSums(pathabun_strat_filt[, 3:(ncol(pathabun_strat_filt))]), '/'), check.names = FALSE) * 100
  MGS_pathabun_strat_filt[, 3:(ncol(MGS_pathabun_strat_filt))] <- data.frame(sweep(MGS_pathabun_strat_filt[, 3:(ncol(MGS_pathabun_strat_filt))],
                                                                                   2, colSums(MGS_pathabun_strat_filt[, 3:(ncol(MGS_pathabun_strat_filt))]), '/'), check.names = FALSE) * 100
  
  pathabun_strat_filt$genus <- gsub("^.*g__", "g__", pathabun_strat_filt$genus)
  
  if(length(which(pathabun_strat_filt$genus == "g__")) > 0) {
    pathabun_strat_filt$genus[which(pathabun_strat_filt$genus == "g__")] <- "unclassified"
  }
  # Aggregate rows by matching pathways / genus.
  MGS_pathabun_strat_filt_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=MGS_pathabun_strat_filt)
  pathabun_strat_filt_genus_sum <- aggregate(. ~ genus + pathway, FUN=sum, data=pathabun_strat_filt)
  
  num_contrib_genera_out <- num_contrib_genera(pred_df = pathabun_strat_filt_genus_sum,
                                                   mgs_df = MGS_pathabun_strat_filt_genus_sum)
  
  num_contrib_genera_out[is.na(num_contrib_genera_out)] <- 0
  
  
  unique_path <- c(pathabun_strat_filt_genus_sum$pathway, MGS_pathabun_strat_filt_genus_sum$pathway)
  
  if(length(which(duplicated(unique_path))) > 0) {
    unique_path <- unique_path[-which(duplicated(unique_path))]
  }
  
  path_spearman <- rep(NA, length(unique_path))
  names(path_spearman) <- unique_path
  
  for(pathway in unique_path) {
    pathway_breakdown <- breakdown_mean_genera_contrib(pathabun_strat_filt_genus_sum, MGS_pathabun_strat_filt_genus_sum, pathway)
    
    if(nrow(pathway_breakdown) <= 1) {
      next 
    }
    
    path_spearman[pathway] <- cor.test(pathway_breakdown$mgs_mean, pathway_breakdown$picrust2_mean, method="spearman")$estimate
  }
  
  return(list(num_contrib=num_contrib_genera_out, pathway_spearman=path_spearman))

}


hmp_strat_pred_vs_mgs <- compare_strat_pred_vs_mgs(picrust2_strat = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                                   mgs_strat = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/mgs/humann2_final_out/humann2_pathabundance_stratified.tsv",
                                                   asv_taxa = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/ASVs_OTUs_taxa_classified/hmp_ASVs_taxa_GG/exported/taxonomy.tsv",
                                                   mgs_col_str2replace = "_R1_R2_cat_Abundance",
                                                   path2rm = path2remove)

mammal_strat_pred_vs_mgs <- compare_strat_pred_vs_mgs(picrust2_strat = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                                   mgs_strat = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/mgs/humann2_final_out/humann2_pathabundance_stratified.tsv",
                                                   asv_taxa = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/ASVs_OTUs_taxa_classified/mammal_ASVs_taxa_GG/exported/taxonomy.tsv",
                                                   mgs_col_str2replace = "_Abundance",
                                                   path2rm = path2remove)

ocean_strat_pred_vs_mgs <- compare_strat_pred_vs_mgs(picrust2_strat = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                                      mgs_strat = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/humann2_final_out/humann2_pathabundance_stratified.tsv",
                                                      asv_taxa = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/ASVs_OTUs_taxa_classified/ocean_ASVs_taxa_GG/exported/taxonomy.tsv",
                                                      mgs_col_str2replace = "_Abundance",
                                                      path2rm = path2remove,
                                                     sample_name_map="/home/gavin/projects/picrust_pipeline/data/validation/ocean/ocean_16S_mgs_sample_links.txt")

blue_strat_pred_vs_mgs <- compare_strat_pred_vs_mgs(picrust2_strat = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                                     mgs_strat = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/humann2_final_out/humann2_pathabundance_stratified.tsv",
                                                     asv_taxa = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/ASVs_OTUs_taxa_classified/blue_ASVs_taxa_GG/exported/taxonomy.tsv",
                                                     mgs_col_str2replace = "_Abundance",
                                                     path2rm = path2remove,
                                                     blueberry_input=TRUE)

strat_pred_vs_mgs <- list(hmp=hmp_strat_pred_vs_mgs,
                          mammal=mammal_strat_pred_vs_mgs,
                          ocean=ocean_strat_pred_vs_mgs,
                          blueberry=blue_strat_pred_vs_mgs)

saveRDS(object = strat_pred_vs_mgs, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/16S_vs_MGS_strat_contrib.rds")

