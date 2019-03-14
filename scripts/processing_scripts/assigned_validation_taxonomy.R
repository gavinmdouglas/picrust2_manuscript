### Commands to assign taxonomy to 16S sequences in validation datasets to be used with PanFP.

library("dada2")
library("Biostrings")

### Assign taxa to HMP.
hmp_fasta_in <- as.character(readDNAStringSet("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/hmp_16S_rep_seqs.fasta"))
hmp_tax <- assignTaxonomy(hmp_fasta_in, refFasta = "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", minBoot=0, multithread=TRUE)
hmp_tax_species <- addSpecies(taxtab = hmp_tax, refFasta = "/scratch/db/dada2_ref_db/rdp_species_assignment_16.fa.gz")

hmp_tax_species[,"Kingdom"] <- paste("k__", hmp_tax_species[,"Kingdom"], sep="")
hmp_tax_species[,"Phylum"] <- paste("p__", hmp_tax_species[,"Phylum"], sep="")
hmp_tax_species[,"Class"] <- paste("c__", hmp_tax_species[,"Class"], sep="")
hmp_tax_species[,"Order"] <- paste("o__", hmp_tax_species[,"Order"], sep="")
hmp_tax_species[,"Family"] <- paste("f__", hmp_tax_species[,"Family"], sep="")
hmp_tax_species[,"Genus"] <- paste("g__", hmp_tax_species[,"Genus"], sep="")
hmp_tax_species[,"Species"] <- paste("s__", hmp_tax_species[,"Species"], sep="")

# Remove all NA values.
hmp_tax_species <- gsub("NA", "", hmp_tax_species)
hmp_taxa_combined <- apply(hmp_tax_species, 1, function(x) paste(x, collapse="; "))

hmp_taxa_out <- data.frame(names(hmp_fasta_in), hmp_taxa_combined)
colnames(hmp_taxa_out) <- c("#OTU ID", "ConsensusLineage")

write.table(x = hmp_taxa_out,
            file = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/hmp_16S_rep_seqs_tax_metadata.tsv", 
            quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

# Assign taxa to iGEM.
mammal_fasta_in <- as.character(readDNAStringSet("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/deblur_output_final/iGEM_16S_rep_seqs.fasta"))
mammal_tax <- assignTaxonomy(mammal_fasta_in, refFasta = "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", minBoot=0, multithread=TRUE)
mammal_tax_species <- addSpecies(taxtab = mammal_tax, refFasta = "/scratch/db/dada2_ref_db/rdp_species_assignment_16.fa.gz")

mammal_tax_species[,"Kingdom"] <- paste("k__", mammal_tax_species[,"Kingdom"], sep="")
mammal_tax_species[,"Phylum"] <- paste("p__", mammal_tax_species[,"Phylum"], sep="")
mammal_tax_species[,"Class"] <- paste("c__", mammal_tax_species[,"Class"], sep="")
mammal_tax_species[,"Order"] <- paste("o__", mammal_tax_species[,"Order"], sep="")
mammal_tax_species[,"Family"] <- paste("f__", mammal_tax_species[,"Family"], sep="")
mammal_tax_species[,"Genus"] <- paste("g__", mammal_tax_species[,"Genus"], sep="")
mammal_tax_species[,"Species"] <- paste("s__", mammal_tax_species[,"Species"], sep="")

# Remove all NA values.
mammal_tax_species <- gsub("NA", "", mammal_tax_species)
mammal_taxa_combined <- apply(mammal_tax_species, 1, function(x) paste(x, collapse="; "))

mammal_taxa_out <- data.frame(names(mammal_fasta_in), mammal_taxa_combined)
colnames(mammal_taxa_out) <- c("#OTU ID", "ConsensusLineage")

write.table(x = mammal_taxa_out,
            file = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/deblur_output_final/iGEM_16S_rep_seqs_tax_metadata.tsv", 
            quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

### soil ###
soil_fasta_in <- as.character(readDNAStringSet("/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S/qiime2_artifacts/soil_16S_rep_seqs.fasta"))
soil_tax <- assignTaxonomy(soil_fasta_in, refFasta = "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", minBoot=0, multithread=TRUE)
soil_tax_species <- addSpecies(taxtab = soil_tax, refFasta = "/scratch/db/dada2_ref_db/rdp_species_assignment_16.fa.gz")

soil_tax_species[,"Kingdom"] <- paste("k__", soil_tax_species[,"Kingdom"], sep="")
soil_tax_species[,"Phylum"] <- paste("p__", soil_tax_species[,"Phylum"], sep="")
soil_tax_species[,"Class"] <- paste("c__", soil_tax_species[,"Class"], sep="")
soil_tax_species[,"Order"] <- paste("o__", soil_tax_species[,"Order"], sep="")
soil_tax_species[,"Family"] <- paste("f__", soil_tax_species[,"Family"], sep="")
soil_tax_species[,"Genus"] <- paste("g__", soil_tax_species[,"Genus"], sep="")
soil_tax_species[,"Species"] <- paste("s__", soil_tax_species[,"Species"], sep="")

# Remove all NA values.
soil_tax_species <- gsub("NA", "", soil_tax_species)
soil_taxa_combined <- apply(soil_tax_species, 1, function(x) paste(x, collapse="; "))

soil_taxa_out <- data.frame(names(soil_fasta_in), soil_taxa_combined)
colnames(soil_taxa_out) <- c("#OTU ID", "ConsensusLineage")

write.table(x = soil_taxa_out,
            file = "/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S/qiime2_artifacts/soil_16S_rep_seqs_metadata.tsv", 
            quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

### ocean ###
ocean_fasta_in <- as.character(readDNAStringSet("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/deblur_output_final/ocean_16S_rep_seqs.fasta"))
ocean_tax <- assignTaxonomy(ocean_fasta_in, refFasta = "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", minBoot=0, multithread=TRUE)
ocean_tax_species <- addSpecies(taxtab = ocean_tax, refFasta = "/scratch/db/dada2_ref_db/rdp_species_assignment_16.fa.gz")

ocean_tax_species[,"Kingdom"] <- paste("k__", ocean_tax_species[,"Kingdom"], sep="")
ocean_tax_species[,"Phylum"] <- paste("p__", ocean_tax_species[,"Phylum"], sep="")
ocean_tax_species[,"Class"] <- paste("c__", ocean_tax_species[,"Class"], sep="")
ocean_tax_species[,"Order"] <- paste("o__", ocean_tax_species[,"Order"], sep="")
ocean_tax_species[,"Family"] <- paste("f__", ocean_tax_species[,"Family"], sep="")
ocean_tax_species[,"Genus"] <- paste("g__", ocean_tax_species[,"Genus"], sep="")
ocean_tax_species[,"Species"] <- paste("s__", ocean_tax_species[,"Species"], sep="")

# Remove all NA values.
ocean_tax_species <- gsub("NA", "", ocean_tax_species)
ocean_taxa_combined <- apply(ocean_tax_species, 1, function(x) paste(x, collapse="; "))

ocean_taxa_out <- data.frame(names(ocean_fasta_in), ocean_taxa_combined)
colnames(ocean_taxa_out) <- c("#OTU ID", "ConsensusLineage")

write.table(x = ocean_taxa_out,
            file = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/deblur_output_final/ocean_16S_rep_seqs_metadata.tsv", 
            quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)


# Assign taxa to blueberry.
blueberry_fasta_in <- as.character(readDNAStringSet("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/deblur_output_exported/blueberry_16S_rep_seqs.fna"))
blueberry_tax <- assignTaxonomy(blueberry_fasta_in, refFasta = "/scratch/db/dada2_ref_db/rdp_train_set_16.fa.gz", minBoot=0, multithread=TRUE)
blueberry_tax_species <- addSpecies(taxtab = blueberry_tax, refFasta = "/scratch/db/dada2_ref_db/rdp_species_assignment_16.fa.gz")

blueberry_tax_species[,"Kingdom"] <- paste("k__", blueberry_tax_species[,"Kingdom"], sep="")
blueberry_tax_species[,"Phylum"] <- paste("p__", blueberry_tax_species[,"Phylum"], sep="")
blueberry_tax_species[,"Class"] <- paste("c__", blueberry_tax_species[,"Class"], sep="")
blueberry_tax_species[,"Order"] <- paste("o__", blueberry_tax_species[,"Order"], sep="")
blueberry_tax_species[,"Family"] <- paste("f__", blueberry_tax_species[,"Family"], sep="")
blueberry_tax_species[,"Genus"] <- paste("g__", blueberry_tax_species[,"Genus"], sep="")
blueberry_tax_species[,"Species"] <- paste("s__", blueberry_tax_species[,"Species"], sep="")

# Remove all NA values.
blueberry_tax_species <- gsub("NA", "", blueberry_tax_species)
blueberry_taxa_combined <- apply(blueberry_tax_species, 1, function(x) paste(x, collapse="; "))

blueberry_taxa_out <- data.frame(names(blueberry_fasta_in), blueberry_taxa_combined)
colnames(blueberry_taxa_out) <- c("#OTU ID", "ConsensusLineage")

write.table(x = blueberry_taxa_out,
            file = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/deblur_output_exported/blueberry_16S_rep_seqs_tax_metadata.tsv", 
            quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)


