### Commands to parse the raw NCBI taxa classifications of all genomes.

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/raw_genome_taxa")

library(ape)

# Function to loop through identified taxa classes and return genome id, superkingdom, class, order, family, genus, species.
return_taxa_levels <- function(in_taxa_df, genome_id) {
  
  if(is.na(in_taxa_df) || dim(in_taxa_df)[1] == 0) {
    return(c(genome_id, rep(NA, 7)))
  }
  
  # Remove any rows that are for duplicated for the same rank.
  dup_ranks <- which(duplicated(in_taxa_df$rank))
  if(length(dup_ranks) > 0) {
    in_taxa_df <- in_taxa_df[-dup_ranks,]
  }
  
  taxa_out <- c()
  taxa_levels <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  
  for(taxa_level in taxa_levels) {
    if(taxa_level %in% in_taxa_df$rank) {
      taxa_out <- c(taxa_out, in_taxa_df[which(in_taxa_df$rank == taxa_level), "name"])
    } else {
      taxa_out <- c(taxa_out, NA)
    }
  }
  
  return(c(genome_id, taxa_out))
  
}


orig_picrust <- readRDS("IMG_16S_classifications_picrust1.rds")
new_picrust <- readRDS("IMG_16S_classifications.rds")
fungi_taxa <- readRDS("fungi_classifications.rds")



orig_picrust_levels <- lapply(names(orig_picrust), function(x) { return_taxa_levels(orig_picrust[[x]][[1]], x) })
orig_picrust_levels_df <- as.data.frame(do.call(rbind, orig_picrust_levels), stringsAsFactors = FALSE)
colnames(orig_picrust_levels_df) <- c("assembly", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

new_picrust_levels <- lapply(names(new_picrust), function(x) { return_taxa_levels(new_picrust[[x]][[1]], x) })
new_picrust_levels_df <- as.data.frame(do.call(rbind, new_picrust_levels), stringsAsFactors = FALSE)
colnames(new_picrust_levels_df) <- c("assembly", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fungi_taxa_levels <- lapply(names(fungi_taxa), function(x) { return_taxa_levels(fungi_taxa[[x]][[1]], x) })
fungi_taxa_levels_df <- as.data.frame(do.call(rbind, fungi_taxa_levels), stringsAsFactors = FALSE)
colnames(fungi_taxa_levels_df) <- c("assembly", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

rownames(fungi_taxa_levels_df) <- fungi_taxa_levels_df$assembly

# Get breakdown by ITS/18S databases separately.
fungi_18S_tree_labels <- read.tree("/home/gavin/projects/picrust_pipeline/fungal_genomes/fungi_alignments/raxml_18S/18S_ssu_align.eukarya.mask.derep.raxml.bestTree")$tip.label
fungi_ITS_tree_labels <- read.tree("/home/gavin/projects/picrust_pipeline/fungal_genomes/fungi_alignments/raxml_ITS/ITS_seqs_best_pass_derep.raxml.bestTree")$tip.label

fungi_18S_taxa_levels_df <- fungi_taxa_levels_df[fungi_18S_tree_labels,]
fungi_ITS_taxa_levels_df <- fungi_taxa_levels_df[fungi_ITS_tree_labels,]

# Get breakdowns of # taxa at each level.
orig_picrust_breakdown <- sapply(orig_picrust_levels_df, function(x) { length(unique(na.omit(x)))  })
new_picrust_breakdown <- sapply(new_picrust_levels_df, function(x) { length(unique(na.omit(x))) })
fungi_18S_breakdown <- sapply(fungi_18S_taxa_levels_df, function(x) { length(unique(na.omit(x))) })
fungi_ITS_breakdown <- sapply(fungi_ITS_taxa_levels_df, function(x) { length(unique(na.omit(x))) })
fungi_ALL_breakdown <- sapply(fungi_taxa_levels_df, function(x) { length(unique(na.omit(x))) })

# Write out tables
write.table(x = orig_picrust_levels_df, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/16S_PICRUSt1_taxa_levels.tsv",
            col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")

write.table(x = new_picrust_levels_df, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/16S_PICRUSt2_taxa_levels.tsv",
            col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")


write.table(x = fungi_18S_taxa_levels_df, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/18S_fungi_taxa_levels.tsv",
            col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")

write.table(x = fungi_ITS_taxa_levels_df, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/ITS_fungi_taxa_levels.tsv",
            col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")
