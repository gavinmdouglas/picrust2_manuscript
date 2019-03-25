### Process tables of BUSCO genome completeness / duplication values to identify fungi reference
### genomes that should be excluded 

convert_fungi_filename_to_id <- function(ids, file_w_ids, file_w_maps) {
 
  # Parse through vector of ids / filenames and return only valid fungi genome ids.
  
   ids_df <- read.table(file_w_ids, header=F, sep="\t", stringsAsFactors = FALSE)

   nonstandard_names_df <- read.table(file_w_maps, header=F, sep="\t", stringsAsFactors = FALSE)
   rownames(nonstandard_names_df) <- nonstandard_names_df$V1
  
   matching_ids_i <- which(ids %in% ids_df$V2)
   non_matching_ids_i <- which(! ids %in% ids_df$V2)

   for(i in non_matching_ids_i) {
     #orig <- ids[i]
     genome_id <- identify_possible_id(id_in=ids[i], possible_ids=ids_df$V2)
     if(! isFALSE(genome_id)) {
       ids[i] <- genome_id
     } else {
       alt_name <- identify_possible_id(id_in=ids[i], possible_ids=nonstandard_names_df$V1)
       if( isFALSE(alt_name)) {
        print(paste("Warning - no possible id found for input", ids[i], sep=" "))
        ids[i] <- NA
       } else {
         ids[i] <- nonstandard_names_df[alt_name, "V2"]
       }
     }

   }
   
   return(ids)
   
}

identify_possible_id <- function(id_in, possible_ids) {
  
  id_split <- strsplit(id_in, "_|\\.")[[1]]
  
  current <- id_split[1]
  id_split <- id_split[-1]
  
  for(i in 1:(length(id_split) + 1)) {
    if(current %in% possible_ids) {
      return(current) 
    }
    current <- paste(current, id_split[i], sep="_")
  }
  return(FALSE)
}

fungi_1000_busco <- read.table("/home/gavin/projects/picrust_pipeline/fungal_genomes/busco_out.tsv",
                               header=T, sep="\t", stringsAsFactors = FALSE)


# Convert completeness and duplication to percentages.
fungi_1000_busco$C_per <- (fungi_1000_busco$C / fungi_1000_busco$Total) * 100
fungi_1000_busco$D_per <- (fungi_1000_busco$D / fungi_1000_busco$Total) * 100

fungi_1000_busco$genome <- gsub("short_summary_", "", fungi_1000_busco$File)
fungi_1000_busco$genome <- gsub("_GeneCatalog_proteins_.*$", "", fungi_1000_busco$genome)

par(mfrow=c(1, 2))
hist(fungi_1000_busco$C_per, breaks=20, ylim=c(0, 500), main="Euk completeness", xlab="Percent Complete BUSCOs", col="grey")
abline(v=70, lwd=2, lty=2)
hist(fungi_1000_busco$D_per, breaks=20, ylim=c(0, 500), xlim=c(0, 100), main="Euk duplicated", xlab="Percent Duplicated BUSCOs", col="grey")
abline(v=10, lwd=2, lty=2)

rows2remove <- unique(c(which(fungi_1000_busco$C_per < 70), which(fungi_1000_busco$D_per > 10)))

genomes2keep <- fungi_1000_busco$genome[-rows2remove]

# Label the genomes by consistent ids.
genomes2keep_renamed <- convert_fungi_filename_to_id(ids=genomes2keep,
                                    file_w_ids="/home/gavin/projects/picrust_pipeline/fungal_genomes/fungi_genome_assembly_id_links_no_mito_public.txt",
                                    file_w_maps="/home/gavin/projects/picrust_pipeline/fungal_genomes/1000_fungi_unique_names_to_ids.txt")

# Remove single NA genome (this one had multiple versions so only latest one used).
genomes2keep_renamed <- genomes2keep_renamed[-which(is.na(genomes2keep_renamed))]

