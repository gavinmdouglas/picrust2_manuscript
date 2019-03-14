setwd("/home/gavin/projects/picrust_pipeline/IMG_phenotypes/raw_phenotype_tables/")

# Read in IMG ids to make table of.
IMG_ids <- read.table("/home/gavin/projects/picrust_pipeline/IMG_pipeline_prep/centroid_ids.txt")$V1

# Read in IMG phenotypes (and which ones are focal and which represent "0").
IMG_pheno_links <- read.table("../IMG_phenotype_links.txt", header=T, sep="\t", stringsAsFactors = FALSE)
IMG_pheno_links$Phenotype <- gsub(" ", "_", IMG_pheno_links$Phenotype)
IMG_pheno_links$Non_phenotype <- gsub(" ", "_", IMG_pheno_links$Non_phenotype)

# Inititalize output table.
pheno_tab <- data.frame(matrix(NA, nrow=length(IMG_ids), ncol=nrow(IMG_pheno_links) + 1))
colnames(pheno_tab) <- c("genome", IMG_pheno_links$Phenotype)
pheno_tab$genome <- IMG_ids
rownames(pheno_tab) <- IMG_ids

# Read in raw files per phenotype that contains genome ids that are predicted
# to have phenotype.
raw_files <- list.files(".")
infiles <- list()
for(x in raw_files) {
  
  infile <- read.table(x, header=TRUE, quote="", comment.char="", stringsAsFactors = FALSE, sep="\t")
  
  expected_phenotype <- gsub(".txt", "", x)
  pheno_name <- expected_phenotype
  expected_phenotype <- gsub("_", " ", expected_phenotype)
  
  if(! all(infile$Phenotype == expected_phenotype)) {
   print(expected_phenotype)
   print(infile$Phenotype[1])
   print("Error - filename phenotype does not match file content in:")
   print(x) 
  }
  
  infiles[[pheno_name]] <- infile
}

# Loop through all primary phenotypes and set all genomes predicted to exhibit them to have setting of 1.
# If there is a secondary phenotype then set all these genomes to be 0 and leave rest as NA.
# Otherwise set all genomes that don't exhibit the phenotype to be 0.
for(i in 1:nrow(IMG_pheno_links)) {
  pheno <- IMG_pheno_links$Phenotype[i]
  secondary <- IMG_pheno_links$Non_phenotype[i]
  
  pheno_genomes <- as.character(infiles[[pheno]]$Genome.ID)

  pheno_genomes <- pheno_genomes[which(pheno_genomes %in% rownames(pheno_tab))]
  
  pheno_tab[pheno_genomes, pheno] <- 1
  
  # If no "negative" phenotype then set all other genomes to be missing.
  if(secondary == "") {
    other_genomes <- which(! rownames(pheno_tab) %in% pheno_genomes)
    pheno_tab[other_genomes, pheno] <- 0
  } else {
    secondary_genomes <- as.character(infiles[[secondary]]$Genome.ID)
    secondary_genomes <- secondary_genomes[which(secondary_genomes %in% rownames(pheno_tab))]
    
    pheno_tab[secondary_genomes, pheno] <- 0
  }
}

# Summarize table to compare counts of all types.
summary_df <- data.frame(
col_na = colSums(is.na(pheno_tab)),
col_1 = colSums(pheno_tab == 1, na.rm = TRUE),
col_0 = colSums(pheno_tab == 0, na.rm = TRUE))
table(rowSums(summary_df))

# Write out table.
write.table(x = pheno_tab, file = "../combined_table/IMG_phenotypes_raw.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)



#### Tests of mean phenotype table.
mean_pheno_table <- read.table("../combined_table/IMG_phenotypes_mean.txt", header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names = FALSE)

pheno_tab_tmp <- pheno_tab[,-1]

pheno_tab_tmp <- data.frame(sapply(pheno_tab_tmp, as.integer), check.names = FALSE)
rownames(pheno_tab_tmp) <- rownames(pheno_tab)

# First check that all non-cluster genomes match the original table exactly.
non_clusters <- grep("cluster", rownames(mean_pheno_table), value=TRUE, invert=TRUE)
mean_pheno_table_non_cluster <- mean_pheno_table[non_clusters, ]
pheno_tab_tmp_non_cluster <- pheno_tab_tmp[non_clusters,]
identical(colnames(pheno_tab_tmp_non_cluster), colnames(mean_pheno_table_non_cluster))
identical(pheno_tab_tmp_non_cluster, mean_pheno_table_non_cluster)

# Here are 5 random test clusters:
test_clusters <- c("650716092-cluster", "2654588122-cluster",
                   "2671180686-cluster", "2551306220-cluster",
                   "648276732-cluster")

#650716092 651053073
#2654588122 2684622647
#2671180686 2698536419 2703718879 2703718893 2703718898
#2551306220 2551306250 2551306265
#648276732 650716091

mean_pheno_table_test <- mean_pheno_table[test_clusters,]

# Now check a couple of examples where the cluster is NA - does this make sense given the constituents?
na_648276732_cluster <- colnames(mean_pheno_table_test)[which(is.na(mean_pheno_table_test["648276732-cluster",]))]

pheno_tab_tmp[c("648276732", "650716091"), na_648276732_cluster]

pheno_tab_tmp[c("2654588122", "2684622647"),]

pheno_tab_tmp[c("2671180686", "2698536419", "2703718879", "2703718893", "2703718898"),]
mean_pheno_table_test["2671180686-cluster",]
rbind(pheno_tab_tmp[c("2671180686", "2698536419", "2703718879", "2703718893", "2703718898"),], mean_pheno_table_test["2671180686-cluster",])
