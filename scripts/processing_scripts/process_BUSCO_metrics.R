### Process tables of BUSCO genome completeness / duplication values to identify eukaryotic reference
<<<<<<< HEAD
### genomes that should be excluded based on eukaryota_odb9 database.

euk_busco <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/busco_out.tsv",
=======
### genomes that should be excluded.
### Ran against 3 BUSCO lineage database: protists_ensembl fungi_odb9 eukaryota_odb9    

# Read in protozoa genome metrics based on protist lineage genes.
protist_busco <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/RefSeq_downloads/protozoa_BUSCO_protist_output.tsv",
                            header=T,
                            sep="\t",
                            stringsAsFactors = FALSE)

fungi_busco <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/RefSeq_downloads/fungi_BUSCO_fungi_output.tsv",
                            header=T,
                            sep="\t",
                            stringsAsFactors = FALSE)

euk_busco <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/RefSeq_downloads/all_BUSCO_euk_output.tsv",
>>>>>>> 42b754e4586a977fbf18bdad273e65279b0d021b
                          header=T,
                          sep="\t",
                          stringsAsFactors = FALSE)


# Completeness based on all euk genes looks much noisier than each individual lineage:
par(mfrow=c(3,1))
hist(euk_busco$C, breaks=50)
<<<<<<< HEAD
=======
hist(protist_busco$C, breaks=50)
hist(fungi_busco$C, breaks=50)
# Fungal genomes seem to be quite complete
>>>>>>> 42b754e4586a977fbf18bdad273e65279b0d021b



# Convert completeness and duplication to percentages.
<<<<<<< HEAD
euk_busco$C_per <- (euk_busco$C / euk_busco$Total) * 100
euk_busco$D_per <- (euk_busco$D / euk_busco$Total) * 100

euk_busco$genome <- gsub("short_summary_GCF_", "", euk_busco$File)
euk_busco$genome <- gsub("_.*$", "", euk_busco$genome)
euk_busco$genome <- paste("GCF", euk_busco$genome, sep="_")

par(mfrow=c(2,1))
hist(euk_busco$C_per, breaks=20, ylim=c(0, 300), main="Euk completeness", xlab="Percent Complete BUSCOs", col="grey")
abline(v=70, lwd=2, lty=2)
hist(euk_busco$D_per, breaks=20, ylim=c(0, 400), main="Euk duplicated", xlab="Percent Duplicated BUSCOs", col="grey")
abline(v=10, lwd=2, lty=2)

rows2remove <- unique(c(which(euk_busco$C_per < 70), which(euk_busco$D_per > 10)))

genomes2keep <- euk_busco$genome[-rows2remove]


# Retained 321 genomes.

write.table(x = genomes2keep,
            file = "/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/passing_euk_genomes.txt",
            col.names = FALSE, row.names = FALSE, quote=FALSE)





######## Processing commands for 1000 fungi database ########
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
      #print(c(orig, ids[i]))
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

=======
fungi_busco$C_per <- (fungi_busco$C / fungi_busco$Total) * 100
fungi_busco$D_per <- (fungi_busco$D / fungi_busco$Total) * 100

protist_busco$C_per <- (protist_busco$C / protist_busco$Total) * 100
protist_busco$D_per <- (protist_busco$D / protist_busco$Total) * 100

fungi_busco$genome <- gsub("short_summary_GCF_", "", fungi_busco$File)
fungi_busco$genome <- gsub("_.*$", "", fungi_busco$genome)
fungi_busco$genome <- paste("GCF", fungi_busco$genome, sep="_")

protist_busco$genome <- gsub("short_summary_GCF_", "", protist_busco$File)
protist_busco$genome <- gsub("_.*$", "", protist_busco$genome)
protist_busco$genome <- paste("GCF", protist_busco$genome, sep="_")


par(mfrow=c(2,2))
hist(fungi_busco$C_per, breaks=20, ylim=c(0, 200), main="Fungi completeness", xlab="Percent Complete BUSCOs", col="grey")
abline(v=80, lwd=2, lty=2)
hist(protist_busco$C_per, breaks=20, ylim=c(0, 12), main="Protozoa completeness", xlab="Percent Complete BUSCOs", col="grey")
abline(v=80, lwd=2, lty=2)
hist(fungi_busco$D_per, breaks=20, ylim=c(0, 200), main="Fungi duplicated", xlab="Percent Duplicated BUSCOs", col="grey")
abline(v=10, lwd=2, lty=2)
hist(protist_busco$D_per, breaks=20, ylim=c(0, 60), main="Protozoa duplicated", xlab="Percent Duplicated BUSCOs", col="grey")
abline(v=10, lwd=2, lty=2)

# In addition to the above cut-offs, also identify genomes with very low numbers of annotated KOs.
euk_ko <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/func_tables/protozoa_fungi_faa_uniref90_hits_ko_regrouped.txt",
                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

par(mfrow=c(1,1))

# Identify number of different KOs per genome.
euk_ko_nozero <- colSums(euk_ko > 0)

hist(euk_ko_nozero, main="", ylim=c(0, 35), xlim=c(0, 3000), col="grey", xlab="Number of different KOs per eukaryotic genome", breaks=50)
abline(v=200, lwd=2, lty=2)

# Identify genomes with at least 200 KOs.
euk_ko_sufficient <- names(euk_ko_nozero)[which(euk_ko_nozero >= 200)]
euk_ko_sufficient <- gsub("GCF_", "", euk_ko_sufficient)
euk_ko_sufficient <- gsub("_.*$", "", euk_ko_sufficient)
euk_ko_sufficient <- paste("GCF", euk_ko_sufficient, sep="_")

# Keep genomes within the above cut-offs.
genomes2keep <- fungi_busco[which(fungi_busco$C_per >= 80 & fungi_busco$D_per <= 10 ), "genome"]
genomes2keep <- c(genomes2keep, protist_busco[which(protist_busco$C_per >= 80 & protist_busco$D_per <= 10 ), "genome"])

genomes2keep <- genomes2keep[which(genomes2keep %in% euk_ko_sufficient)]

# Retained 205 genomes after this point.

write.table(x = genomes2keep,
            file = "/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/passing_euk_genomes.txt",
            col.names = FALSE, row.names = FALSE, quote=FALSE)
>>>>>>> 42b754e4586a977fbf18bdad273e65279b0d021b
