### Function to output a file of all the genome ids within a given taxa classification (on a single line).
### (given a set of genome ids and a taxa table).
### These tables will then be parsed for LOOCV analyses.

return_taxa_id_vector <- function(in_taxa) {
 
  taxa_id_vec <- c()
  
  tax_cols <- c("assembly", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for(col in tax_cols) {
    
    unique_categories <- unique(in_taxa[, col])
    
    for(category in na.omit(unique_categories)) {

      matching_assemblies <- in_taxa[which(in_taxa[, col] == category), "assembly"]
      
      category <- gsub(" ", "_", category)
      
      id_links <- paste(c(col, category, matching_assemblies), collapse="\t")
      
      taxa_id_vec <- c(taxa_id_vec, id_links)
    }
  }
  
  return(taxa_id_vec)
}


taxa_16S <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/16S_PICRUSt2_taxa_levels_NA_filled.tsv",
                             sep="\t", header=T, stringsAsFactors = FALSE, quote="")


taxa_16S_groupings <- return_taxa_id_vector(taxa_16S)

write.table(taxa_16S_groupings, "/home/gavin/projects/picrust_pipeline/16S_LOOCV/16S_taxa_groupings.tsv", col.names=F, row.names=F,
            sep="\t", quote=FALSE)
