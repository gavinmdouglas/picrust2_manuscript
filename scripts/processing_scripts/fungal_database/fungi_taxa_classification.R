# Needed to run this locally since NCBI denied connection from vulcan.

# Function to try trimming end of genome name to see if that can be classified.
get_taxa_from_str <- function(str_in) {
  
  # Get number of spaces in input string and add 1.
  num_loop <- length(grep(" ", unlist(strsplit(str_in, "")))) + 1
  
  for(i in 1:num_loop) {
    taxa_out <- tryCatch(classification(c(str_in), db="ncbi", out_type="summary"), error=function(err) NA)
    # If got matching taxa then return.
    if(! is.na(taxa_out)) {
      return(taxa_out)
    }
    
    str_in <- gsub(" [^ ]*$", "", str_in)
    
  }
  
  return(NA)
  
}

setwd("/Users/Gavin/Dropbox/work/Langille/PICRUSt2/ms_code_for_github/fungi_taxa_classification")

library("taxizedb") # ‘0.1.7.9601’
library("ape")
library("parallel")

# Get final fungi ids from trees.
fungi_ITS_tree_labels <- read.tree("ITS_seqs_best_pass_derep.raxml.bestTree")$tip.label
fungi_18S_tree_labels <- read.tree("18S_ssu_align.eukarya.mask.derep.raxml.bestTree")$tip.label

fungi_ITS_tree_labels <- gsub("-cluster", "", fungi_ITS_tree_labels)
fungi_18S_tree_labels <- gsub("-cluster", "", fungi_18S_tree_labels)
label_union <- unique(c(fungi_ITS_tree_labels, fungi_18S_tree_labels))

fungi_species_to_id <- read.table("fungi_genome_assembly_id_links_no_mito_public.txt", header=F, sep="\t", stringsAsFactors = FALSE, comment.char="", quote="")
rownames(fungi_species_to_id) <- fungi_species_to_id$V2
fungi_species_to_id_subset <- fungi_species_to_id[label_union,]

# Change label for Naematella encephela UCDFST 68-887.2 v1.0 to be Naematelia encephaliformis
fungi_species_to_id_subset["Treen1", "V1"] <- "Naematelia encephaliformis"

# Change label for "Phlebia radiata Fr. (isolate 79, FBCC0043)" to be Phlebia radiata
fungi_species_to_id_subset["Phlrad1", "V1"] <- "Phlebia radiata"

fungi_classifications <- mclapply(fungi_species_to_id_subset$V1, get_taxa_from_str, mc.cores=4)

names(fungi_classifications) <- fungi_species_to_id_subset$V2

saveRDS(object = fungi_classifications, file = "fungi_classifications.rds")
