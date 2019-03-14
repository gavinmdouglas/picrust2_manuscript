### Comparing number of genera represented by ASVs and OTUs in each 16S validation dataset.
### Basic proof that ASVs enable improved taxonomic diversity.

setwd("/home/gavin/projects/picrust_pipeline/data/validation/")


get_unclassified_genus <- function(taxa_vec) {
  
  taxa_level <- c("D_0__", ";D_1__", ";D_2__", ";D_3__", ";D_4__", ";D_5__")
  
  taxa_vec <- gsub(";D_6__.*$", "", taxa_vec)

  for(tax_label in taxa_level) {
    
    matching_i <- grep(tax_label, taxa_vec, invert=TRUE)
    
    missing2add <- paste(tax_label, "Unclassified", sep="")
    
    taxa_vec[matching_i] <- paste(taxa_vec[matching_i], missing2add, sep="")
    
  }
  
  taxa_vec <- gsub("uncultured bacterium", "Unclassified", taxa_vec)
  taxa_vec <- gsub("uncultured", "Unclassified", taxa_vec)
  taxa_vec <- gsub("gut metagenome", "Unclassified", taxa_vec)
  taxa_vec <- gsub("metagenome", "Unclassified", taxa_vec)
  taxa_vec <- gsub("Unclassified organism", "Unclassified", taxa_vec)
  
  return(taxa_vec)
  
}



blue_ASV_taxa <- read.table("blueberry/16S/ASVs_OTUs_taxa_classified/blue_ASVs_taxa/exported/taxonomy.tsv", header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
blue_OTU_taxa <- read.table("blueberry/16S/ASVs_OTUs_taxa_classified/blue_OTUs_taxa/exported/taxonomy.tsv", header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

blue_ASV_taxa$genus <- get_unclassified_genus(blue_ASV_taxa$Taxon)
blue_OTU_taxa$genus <- get_unclassified_genus(blue_OTU_taxa$Taxon)

hmp_ASV_taxa <- read.table("hmp/16S/ASVs_OTUs_taxa_classified/hmp_ASVs_taxa/exported/taxonomy.tsv", header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
hmp_OTU_taxa <- read.table("hmp/16S/ASVs_OTUs_taxa_classified/hmp_OTUs_taxa/exported/taxonomy.tsv", header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

hmp_ASV_taxa$genus <- get_unclassified_genus(hmp_ASV_taxa$Taxon)
hmp_OTU_taxa$genus <- get_unclassified_genus(hmp_OTU_taxa$Taxon)

rownames(hmp_ASV_taxa)[grep("Pseudomonas", hmp_ASV_taxa$genus)]
rownames(hmp_OTU_taxa)[grep("Pseudomonas", hmp_OTU_taxa$genus)]


blue_OTU_taxa_single <- names(which(table(blue_OTU_taxa$genus) == 1))

"D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales;D_4__Rhodomicrobiaceae;D_5__Rhodomicrobium"
rownames(blue_ASV_taxa)[which(blue_ASV_taxa$genus == "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales;D_4__Rhodomicrobiaceae;D_5__Rhodomicrobium")]
rownames(blue_OTU_taxa)[which(blue_OTU_taxa$genus == "D_0__Bacteria;D_1__Proteobacteria;D_2__Alphaproteobacteria;D_3__Rhizobiales;D_4__Rhodomicrobiaceae;D_5__Rhodomicrobium")]

"b08a057837436b2ee5b89eaa6f8015b7" "d5a21822a3e5dc64b9bd258eb22df501"
"1108989"


library(ape)

blue_ASV_tree <- read.tree(file = "blueberry/16S/picrust2_pipeline/picrust2_full_output_pipeline_2.1.0-b/out.tre")
blue_OTU_tree <- read.tree(file = "blueberry/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/out.tre")


library(phytools)
getMRCA(phy, tip)
blue_ASV_Rhodomicrobium_placement <- getMRCA(blue_ASV_tree, tip=c("b08a057837436b2ee5b89eaa6f8015b7", "d5a21822a3e5dc64b9bd258eb22df501"))

blue_ASV_tree_Rhodomicrobium <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_Rhodomicrobium_placement, collapse.singles=TRUE)



blueberry_denovo_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_pipeline_2.1.0-b/marker_predicted_and_nsti.tsv",
                                      header=T, sep="\t", row.names=1)

blueberry_gg_marker <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_pipeline/picrust2_full_output_picrust1-prep_2.1.0-b/marker_predicted_and_nsti.tsv",
                                  header=T, sep="\t", row.names=1)


Hyphomicrobiaceae_ids <- as.character(read.table("/home/gavin/tmp/Hyphomicrobiaceae_ids.txt", stringsAsFactors = FALSE, header=FALSE)$V1)

non_matched_i <- which(! Hyphomicrobiaceae_ids %in% blue_ASV_tree$tip.label)

Hyphomicrobiaceae_ids[non_matched_i] <- gsub("$", "-cluster", Hyphomicrobiaceae_ids[non_matched_i])




blue_ASV_tree_Rhodomicrobium_related_node <- getMRCA(blue_ASV_tree, c("2728368988", "2576861659"))
blue_ASV_tree_Rhodomicrobium_related <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Rhodomicrobium_related_node, collapse.singles=FALSE)

blue_OTU_tree_Rhodomicrobium_related_node <- getMRCA(blue_OTU_tree, c("2728368988", "2576861659"))
blue_OTU_tree_Rhodomicrobium_related <- extract.clade(phy = blue_OTU_tree, node = blue_OTU_tree_Rhodomicrobium_related_node, collapse.singles=FALSE)


blue_ASV_tree_Rhodomicrobium_node <- getMRCA(blue_ASV_tree, c("649633090", "2576861659"))

blue_ASV_tree_Hyphomicrobiaceae_node <- getMRCA(blue_ASV_tree, as.character(Hyphomicrobiaceae_ids))


blue_ASV_tree_Rhodomicrobium <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Rhodomicrobium_node, collapse.singles=TRUE)

blue_ASV_tree_Hyphomicrobiaceae <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Hyphomicrobiaceae_node, collapse.singles=TRUE)


blue_ASV_tree_test <- extract.clade(phy = blue_ASV_tree, node = 26407, collapse.singles=TRUE)


Desulfovibrionaceae_ids <- as.character(read.table("/home/gavin/tmp/Desulfovibrionaceae_ids.txt", stringsAsFactors = FALSE, header=FALSE)$V1)
non_matched_i <- which(! Desulfovibrionaceae_ids %in% blue_ASV_tree$tip.label)
Desulfovibrionaceae_ids[non_matched_i] <- gsub("$", "-cluster", Desulfovibrionaceae_ids[non_matched_i])
blue_ASV_tree_Desulfovibrionaceae_node <- getMRCA(blue_ASV_tree, as.character(Desulfovibrionaceae_ids))
blue_ASV_tree_Desulfovibrionaceae <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Desulfovibrionaceae_node, collapse.singles=TRUE)


Agrobacterium_ids <- as.character(read.table("/home/gavin/tmp/Agrobacterium_ids.txt", stringsAsFactors = FALSE, header=FALSE)$V1)
non_matched_i <- which(! Agrobacterium_ids %in% blue_ASV_tree$tip.label)
Agrobacterium_ids[non_matched_i] <- gsub("$", "-cluster", Agrobacterium_ids[non_matched_i])
blue_ASV_tree_Agrobacterium_node <- getMRCA(blue_ASV_tree, as.character(Agrobacterium_ids)[c(1,2,3,4)])
blue_ASV_tree_Agrobacterium <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Agrobacterium_node, collapse.singles=TRUE)
blue_ASV_tree_Agrobacterium_node <- getMRCA(blue_ASV_tree, as.character(Agrobacterium_ids)[c(7, 8)])


Bacillus_ASVs <- rownames(blue_ASV_taxa)[grep("Nocardia",  blue_ASV_taxa$Taxon)]
Bacillus_OTUs <- rownames(blue_OTU_taxa)[grep("Nocardia",  blue_OTU_taxa$Taxon)]
blue_ASV_tree_Bacillus_node <- getMRCA(blue_ASV_tree, Bacillus_ASVs[c(1,4)])
blue_ASV_tree_Bacillus <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Bacillus_node, collapse.singles=TRUE)


Rhizobacter_ASVs <- rownames(blue_ASV_taxa)[grep("Rhizobacter",  blue_ASV_taxa$Taxon)]
Rhizobacter_OTUs <- rownames(blue_OTU_taxa)[grep("Rhizobacter",  blue_OTU_taxa$Taxon)]
blue_ASV_tree_Rhizobacter_node <- getMRCA(blue_ASV_tree, Rhizobacter_ASVs)
blue_ASV_tree_Rhizobacter <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_Rhizobacter_node, collapse.singles=TRUE)


clade_ids <- as.character(read.table("/home/gavin/tmp/clade_tips.txt", stringsAsFactors = FALSE, header=FALSE)$V1)
non_matched_i <- which(! clade_ids %in% blue_ASV_tree$tip.label)
clade_ids[non_matched_i] <- gsub("$", "-cluster", clade_ids[non_matched_i])
blue_ASV_tree_clade_node <- getMRCA(blue_ASV_tree, clade_ids)
blue_ASV_tree_clade <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_clade_node, collapse.singles=TRUE)


blue_OTU_tree_clade_node <- getMRCA(blue_OTU_tree, clade_ids)
blue_OTU_tree_clade <- extract.clade(phy = blue_OTU_tree, node = blue_OTU_tree_clade_node, collapse.singles=TRUE)



blue_ASV_tree_clade2_node <- getMRCA(blue_ASV_tree, c("2734481911", "2667527978"))
blue_ASV_tree_clade2 <- extract.clade(phy = blue_ASV_tree, node = blue_ASV_tree_clade2_node, collapse.singles=TRUE)
plot(blue_ASV_tree_clade2)

blue_OTU_tree_clade2_node <- getMRCA(blue_OTU_tree, c("2734481911", "2667527978"))
blue_OTU_tree_clade2 <- extract.clade(phy = blue_OTU_tree, node = blue_OTU_tree_clade2_node, collapse.singles=TRUE)
plot(blue_OTU_tree_clade2)

# TESTING:

get_genus_diff <- function(ASV_file, OTU_file) {

  ASV_taxa <- read.table(ASV_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
  OTU_taxa <- read.table(OTU_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

  ASV_taxa$genus <- get_unclassified_genus(ASV_taxa$Taxon)
  OTU_taxa$genus <- get_unclassified_genus(OTU_taxa$Taxon)

  ASV_taxa_noUnclassified <- ASV_taxa[-grep("Unclassified", ASV_taxa$genus), ]
  OTU_taxa_noUnclassified <- OTU_taxa[-grep("Unclassified", OTU_taxa$genus), ]

  diff_underlying_genera <- c()

  unique_genera <- unique(c(ASV_taxa_noUnclassified$genus, OTU_taxa_noUnclassified$genus))

  for(genus in unique_genera) {

    genus_diff <- length(which(ASV_taxa_noUnclassified$genus == genus)) - length(which(OTU_taxa_noUnclassified$genus == genus))

    diff_underlying_genera <- c(diff_underlying_genera, genus_diff)
  }

  names(diff_underlying_genera) <- unique_genera

  return(diff_underlying_genera)
}


get_genus_ratio <- function(ASV_file, OTU_file) {

  ASV_taxa <- read.table(ASV_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
  OTU_taxa <- read.table(OTU_file, header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

  ASV_taxa$genus <- get_unclassified_genus(ASV_taxa$Taxon)
  OTU_taxa$genus <- get_unclassified_genus(OTU_taxa$Taxon)

  ASV_taxa_noUnclassified <- ASV_taxa[-grep("Unclassified", ASV_taxa$genus), ]
  OTU_taxa_noUnclassified <- OTU_taxa[-grep("Unclassified", OTU_taxa$genus), ]

  ratio_underlying_genera <- c()

  unique_genera <- unique(c(ASV_taxa_noUnclassified$genus, OTU_taxa_noUnclassified$genus))

  for(genus in unique_genera) {

    genus_ratio <- (length(which(ASV_taxa_noUnclassified$genus == genus)) + 1) / (length(which(OTU_taxa_noUnclassified$genus == genus)) + 1)

    ratio_underlying_genera <- c(ratio_underlying_genera, genus_ratio)
  }

  names(ratio_underlying_genera) <- unique_genera

  return(ratio_underlying_genera)
}

hmp_genus_diff <- get_genus_diff("hmp/16S/ASVs_OTUs_taxa_classified/hmp_ASVs_taxa/exported/taxonomy.tsv",
                                 "hmp/16S/ASVs_OTUs_taxa_classified/hmp_OTUs_taxa/exported/taxonomy.tsv")

mammal_genus_diff <- get_genus_diff("iGEM/16S/ASVs_OTUs_taxa_classified/mammal_ASVs_taxa/exported/taxonomy.tsv",
                                    "iGEM/16S/ASVs_OTUs_taxa_classified/mammal_OTUs_taxa/exported/taxonomy.tsv")

ocean_genus_diff <- get_genus_diff("ocean/16S/ASVs_OTUs_taxa_classified/ocean_ASVs_taxa/exported/taxonomy.tsv",
                                   "ocean/16S/ASVs_OTUs_taxa_classified/ocean_OTUs_taxa/exported/taxonomy.tsv")


blueberry_genus_diff <- get_genus_diff("blueberry/16S/ASVs_OTUs_taxa_classified/blue_ASVs_taxa/exported/taxonomy.tsv",
                                       "blueberry/16S/ASVs_OTUs_taxa_classified/blue_OTUs_taxa/exported/taxonomy.tsv")




hmp_genus_ratio <- get_genus_ratio("hmp/16S/ASVs_OTUs_taxa_classified/hmp_ASVs_taxa/exported/taxonomy.tsv",
                                   "hmp/16S/ASVs_OTUs_taxa_classified/hmp_OTUs_taxa/exported/taxonomy.tsv")

mammal_genus_ratio <- get_genus_ratio("iGEM/16S/ASVs_OTUs_taxa_classified/mammal_ASVs_taxa/exported/taxonomy.tsv",
                                      "iGEM/16S/ASVs_OTUs_taxa_classified/mammal_OTUs_taxa/exported/taxonomy.tsv")

ocean_genus_ratio <- get_genus_ratio("ocean/16S/ASVs_OTUs_taxa_classified/ocean_ASVs_taxa/exported/taxonomy.tsv",
                                     "ocean/16S/ASVs_OTUs_taxa_classified/ocean_OTUs_taxa/exported/taxonomy.tsv")


blueberry_genus_ratio <- get_genus_ratio("blueberry/16S/ASVs_OTUs_taxa_classified/blue_ASVs_taxa/exported/taxonomy.tsv",
                                         "blueberry/16S/ASVs_OTUs_taxa_classified/blue_OTUs_taxa/exported/taxonomy.tsv")

boxplot(hmp_genus_diff + 1, mammal_genus_diff + 1, ocean_genus_diff + 1, blueberry_genus_diff + 1, log="y")

boxplot(hmp_genus_ratio, mammal_genus_ratio, ocean_genus_ratio, blueberry_genus_ratio, log="y")
abline(h=1)

