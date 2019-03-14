### Identify eukaryota to download in RefSeq (all eukaryota except for plants and animals) and download the protein
### and genome sequences for each one.

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/raw_genome_taxa")

# Function to loop through identified taxa classes and return genome id, superkingdom, class, order, family, genus, species.
return_taxa_levels_refseq_modified <- function(in_taxa_df, genome_id) {
  
  if(is.na(in_taxa_df)) { return(c(genome_id, rep(NA, 8))) }
  
  # Remove any rows that are "no_rank".
  no_rank_rows <- which(in_taxa_df$rank == "no_rank")
  if(length(no_rank_rows) > 0) {
    in_taxa_df <- in_taxa_df[-no_rank_rows,,drop=FALSE]
  }
  
  # Return NAs if nothing left now.
  if(dim(in_taxa_df)[1] == 0) { return(c(genome_id, rep(NA, 8))) }
  
  # Remove any rows that are for duplicated for the same rank.
  dup_ranks <- which(duplicated(in_taxa_df$rank))
  if(length(dup_ranks) > 0) {
    in_taxa_df <- in_taxa_df[-dup_ranks,,drop=FALSE]
  }
  
  taxa_out <- c()
  taxa_levels <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  for(taxa_level in taxa_levels) {
    if(taxa_level %in% in_taxa_df$rank) {
      taxa_out <- c(taxa_out, in_taxa_df[which(in_taxa_df$rank == taxa_level), "name"])
    } else {
      taxa_out <- c(taxa_out, NA)
    }
  }
  
  return(c(genome_id, taxa_out))
  
}

all_refseq <- readRDS("refseq_summary_taxa.rds")

all_refseq_levels <- lapply(names(all_refseq), function(x) { return_taxa_levels_refseq_modified(all_refseq[[x]], x) })


# Creating dataframe manually since it is extremely slow otherwise (e.g. when rbind.data.frame is used).
all_refseq_levels_df <- data.frame(matrix(NA, nrow=142091, ncol=9))
for(i in 1:nrow(all_refseq_levels_df)) {
  all_refseq_levels_df[i,] <- all_refseq_levels[[i]]
}
colnames(all_refseq_levels_df) <- c("taxid", "Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#write.table(x = all_refseq_levels_df, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/all_RefSeq_taxa_levels_30Oct2018.tsv",
#            col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")

all_refseq_levels_df_eukaryota <- all_refseq_levels_df[which(all_refseq_levels_df$Superkingdom=="Eukaryota"),]

# remove Metazoa and Streptophyta
all_refseq_levels_df_eukaryota_filt <- all_refseq_levels_df_eukaryota[-which(all_refseq_levels_df_eukaryota$Kingdom == "Metazoa"),]
all_refseq_levels_df_eukaryota_filt <- all_refseq_levels_df_eukaryota_filt[-which(all_refseq_levels_df_eukaryota_filt$Phylum == "Streptophyta"),]

#write.table(x = all_refseq_levels_df_eukaryota_filt, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/all_RefSeq_taxa_levels_30Oct2018_eukaryota_microbes.tsv",
#            col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# Read in assembly summary from RefSeq.
refseq_assembly <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/assembly_summary_refseq_30Oct2018.txt",
                              header=T, sep="\t", comment.char="", skip=1, quote="")

refseq_assembly_eukaryota_microbes <- refseq_assembly[which(refseq_assembly$taxid %in% all_refseq_levels_df_eukaryota_filt$taxid),]

