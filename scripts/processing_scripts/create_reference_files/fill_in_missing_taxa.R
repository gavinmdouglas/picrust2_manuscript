### Commands to fill in placeholder taxa levels for all NA taxa that are identified otherwise to the phylum level.
### Also append "-cluster" to assemblies that represent multiple assemblies.

fill_in_NA_taxa <- function(tax_in) {
  
  missing_rows <- which(rowSums(is.na(tax_in)) > 0 )
  
  for(missing_row in missing_rows) {
    
    # Skip if Superkingdom or Phylum are missing.
    if(is.na(tax_in[missing_row, "Superkingdom"]) | is.na(tax_in[missing_row, "Phylum"])) {
      next
    }
    
    labels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    for(i in 2:6) {
      
      current_label = labels[i]
      past_label = labels[i-1]
      
      if(is.na(tax_in[missing_row, current_label])) {
        tax_in[missing_row, current_label] <- paste(tax_in[missing_row, past_label], "X", sep="_")
      }
      
    }
    
  }
  return(tax_in)
}

library(ape)

# Read in input files.
tree_16S <- read.tree("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/pro_ref/pro_ref.tre")

tax_16S <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/16S_PICRUSt2_taxa_levels.tsv", 
                      header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

tree_18S <- read.tree("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/fungi_18S/fungi_18S.tre")
tax_18S <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/18S_fungi_taxa_levels.tsv", 
                      header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

tree_ITS <- read.tree("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/fungi_ITS/fungi_ITS.tre")
tax_ITS <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/ITS_fungi_taxa_levels.tsv", 
                      header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

tax_16S_filled <- fill_in_NA_taxa(tax_16S)
tax_18S_filled <- fill_in_NA_taxa(tax_18S)
tax_ITS_filled <- fill_in_NA_taxa(tax_ITS)


# Add "-cluster" to assembly ids where needed and write output file.
missing_assemblies_i_16S <- which(! tax_16S_filled$assembly %in% tree_16S$tip.label)
tax_16S_filled[missing_assemblies_i_16S, "assembly"] <- gsub("$", "-cluster", tax_16S_filled[missing_assemblies_i_16S, "assembly"])
write.table(x = tax_16S_filled, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/16S_PICRUSt2_taxa_levels_NA_filled.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

missing_assemblies_i_18S <- which(! tax_18S_filled$assembly %in% tree_18S$tip.label)
tax_18S_filled[missing_assemblies_i_18S, "assembly"] <- gsub("$", "-cluster", tax_18S_filled[missing_assemblies_i_18S, "assembly"])
write.table(x = tax_18S_filled, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/18S_PICRUSt2_taxa_levels_NA_filled.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

missing_assemblies_i_ITS <- which(! tax_ITS_filled$assembly %in% tree_ITS$tip.label)
tax_ITS_filled[missing_assemblies_i_ITS, "assembly"] <- gsub("$", "-cluster", tax_ITS_filled[missing_assemblies_i_ITS, "assembly"])
write.table(x = tax_ITS_filled, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/ITS_PICRUSt2_taxa_levels_NA_filled.tsv",
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

