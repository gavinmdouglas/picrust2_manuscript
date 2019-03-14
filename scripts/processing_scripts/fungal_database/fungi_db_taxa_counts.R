rm(list=ls())

library(ape)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/")

tree_18S <- read.tree("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/fungi/fungi_ref_18S/fungi_18S.tre")

tree_labels <- gsub("-cluster", "", tree_18S$tip.label)

taxa_18S <- read.table("18S_fungi_taxa_levels.tsv", header = TRUE, sep="\t", quote = "", comment.char="", stringsAsFactors = FALSE)

taxa_18S <- taxa_18S[which(taxa_18S$assembly %in% tree_labels), ]

sapply(taxa_18S, function(x) { no_na = na.omit(x); length(unique(no_na)) })



tree_ITS <- read.tree("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/fungi/fungi_ref_ITS/fungi_ITS.tre")

tree_labels <- gsub("-cluster", "", tree_ITS$tip.label)

taxa_ITS <- read.table("ITS_fungi_taxa_levels.tsv", header = TRUE, sep="\t", quote = "", comment.char="", stringsAsFactors = FALSE)

taxa_ITS <- taxa_ITS[which(taxa_ITS$assembly %in% tree_labels), ]

sapply(taxa_ITS, function(x) { no_na = na.omit(x); length(unique(no_na)) })
