### Get full NCBI taxa levels (e.g. Phylum, Class...) based on species name.

### Do this for all protozoa/fungi in the 18S/ITS databases.

### Also do this for all 16S sequences!

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

setwd("/Users/Gavin/Dropbox/work/Langille/PICRUSt2/ms_code_for_github/")

library("taxizedb") # ‘0.1.7.9601’
library("ape")
library("parallel")

# Read in tip labels of interest.
tree_18S <- read.tree("18S_ssu_align.eukarya.mask.derep.raxml.bestTree")
tree_ITS <- read.tree("ITS_seqs_best_pass_derep.raxml.bestTree")

# Read in RefSeq assembly info, which contains taxid linked to ids.
protozoa_summary <- read.table("protozoa_assembly_summary_13Feb2018.txt",
                               header=T,
                               sep="\t",
                               skip=1,
                               comment.char="",
                               row.names=1,
                               quote="")

fungi_summary <- read.table("fungi_assembly_summary_13Feb2018.txt",
                               header=T,
                               sep="\t",
                               skip=1,
                               comment.char="",
                               row.names=1,
                            quote="")

# Combine these 2 tables.
euk_summary <- rbind(protozoa_summary, fungi_summary)

euk_tips <- unique(c(tree_18S$tip.label, tree_ITS$tip.label))

# Remove "_cluster" from the name of 2 tips
euk_tips <- gsub("_cluster", "", euk_tips)

# Subset to just ids in at least 1 database.
euk_summary_subset <- euk_summary[euk_tips,]

euk_taxa <- lapply(euk_summary_subset$taxid, classification, db="ncbi")

# Prokaryotic taxonomy is a little messier since sometimes the IMG genome name isn't in RefSeq.
# Need to keep trimming off the end until there is a taxa match.
IMG_species <- read.table("16S_IMG_species.tsv",
                          header=T,
                          sep="\t",
                          stringsAsFactors = FALSE,
                          quote="",
                          comment.char = "")

tree_16S <- read.tree("reference.tre")

# Remove "_cluster" from tip names.
tip_names <- gsub("_cluster", "", tree_16S$tip.label)

# Subset to only those genomes in 16S tree.
IMG_species_subset <- IMG_species[which(IMG_species$IMG_id %in% tip_names),]

IMG_16S_classifications <- mclapply(IMG_species_subset$Species, get_taxa_from_str, mc.cores=4)

names(IMG_16S_classifications) <- IMG_species_subset$IMG_id

### Also get taxa classifications for genomes just in PICRUSt1 tree (which includes some now deleted genomes).
orig_picrust_ids <- read.table("gg_13_5_img_fixed.txt", header=F, sep="\t")

# Read in table of deleted genome info downloaded from IMG.
img_deleted <- read.table("taxondelete43705_11-sep-2018.txt", header=T, sep="\t", stringsAsFactors = F,
                          comment.char = "", quote = "")

# ids overlapping with current table.
IMG_species_picrust1_subset <- IMG_species[which(IMG_species$IMG_id %in% orig_picrust_ids$V2),]

# Get subset for overlapping deleting genomes.
img_deleted_picrust1_overlap <- img_deleted[which(img_deleted$Genome.ID %in% orig_picrust_ids$V2), c("Genome.ID", "Genome.Name")]

colnames(img_deleted_picrust1_overlap) <- c("IMG_id", "Species")

### 9 genomes are mysteriously missing from new and deleted tables. Found the taxonomy of these genomes by hand (where possible - NA if not).
missing_ids <- orig_picrust_ids[which(! orig_picrust_ids$V2 %in% IMG_species_picrust1_subset_and_del$IMG_id), "V2"]
missing_id_taxa <- c("Campylobacter coli LMG 9860", "Corynebacterium casei UCMA 3821", "Elizabethkingia anophelis Ag1",
                     "NA", "Haemophilus haemolyticus M21621", "Mycoplasma ovipneumoniae SC01", "Citreicella sp. 357",
                     "Ruegeria conchae TW15", "Pseudomonas fluorescens R124")

p1_missing_ids <- data.frame(IMG_id = missing_ids, Species = missing_id_taxa)

IMG_species_picrust1_subset_and_del <- rbind(IMG_species_picrust1_subset, img_deleted_picrust1_overlap, p1_missing_ids)

IMG_16S_classifications_picrust1 <- mclapply(IMG_species_picrust1_subset_and_del$Species, get_taxa_from_str, mc.cores=4)

names(IMG_16S_classifications_picrust1) <- IMG_species_picrust1_subset_and_del$IMG_id

# Save these taxa lists to RDS files.
saveRDS(object = euk_taxa, file = "euk_taxa.rds")
saveRDS(object = IMG_16S_classifications, file = "IMG_16S_classifications.rds")
saveRDS(object = IMG_16S_classifications_picrust1, file = "IMG_16S_classifications_picrust1.rds")

