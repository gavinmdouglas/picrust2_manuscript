### Read in functional databases and get basic counts of functions and measures of table sparsity.
### Calculate trait depth using castor to get an idea of how conserved each trait is on tree.
### Read in old database files (from early PICRUSt2 beta, which was based on PyNAST alignment) for basic sanity checks.

setwd("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic")
source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

library(ape)
library(castor)
library(ggplot2)

# First read in latest functional databases:
rrna <- read.table(gzfile("16S.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
ec <- read.table(gzfile("ec.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
ko <- read.table(gzfile("ko.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
cog <- read.table(gzfile("cog.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
pfam <- read.table(gzfile("pfam.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
tigrfam <- read.table(gzfile("tigrfam.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

# Get basic summary stats per database.
characterize_db(rrna)
characterize_db(cog)
characterize_db(ec)
characterize_db(ko)
characterize_db(pfam)
characterize_db(tigrfam)

# Read in reference tree.
ref_tree <- read.tree("reference.tre")

# Get trait depth (after re-coding to binary presence absence for all traits).
rrna_td <- calc_trait_depth(ref_tree, rrna, 1)
cog_td <- calc_trait_depth(ref_tree, cog, 40)
ec_td <- calc_trait_depth(ref_tree, ec, 40)
ko_td <- calc_trait_depth(ref_tree, ko, 40)
pfam_td <- calc_trait_depth(ref_tree, pfam, 40)
tigrfam_td <- calc_trait_depth(ref_tree, tigrfam, 40)

# Parse out mean and sd in trait depth for each dataset.
rrna_td_summary <- data.frame(t(data.frame(lapply(rrna_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
cog_td_summary <- data.frame(t(data.frame(lapply(cog_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
ec_td_summary <- data.frame(t(data.frame(lapply(ec_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
ko_td_summary <- data.frame(t(data.frame(lapply(ko_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
pfam_td_summary <- data.frame(t(data.frame(lapply(pfam_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
tigrfam_td_summary <- data.frame(t(data.frame(lapply(tigrfam_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))

# Add new column for what category the functions are in.
rrna_td_summary$category <- "16S"
cog_td_summary$category <- "COG"
ec_td_summary$category <- "E.C. Number"
ko_td_summary$category <- "KO"
pfam_td_summary$category <- "Pfam"
tigrfam_td_summary$category <- "TIGRFAM"

# Concatenate the rows into 1 final table.
td_summaries <- rbind(rrna_td_summary, cog_td_summary, ec_td_summary, ko_td_summary, pfam_td_summary, tigrfam_td_summary)

colnames(td_summaries) <- c("mean", "sd", "category")

# Set factors to be in order from most conserved to least.
td_summaries$category <- factor(td_summaries$category, levels=c("16S", "COG", "TIGRFAM", "E.C. Number", "Pfam", "KO"))

dodge <- position_dodge(width = 1)

ggplot(td_summaries, aes(category, mean)) + geom_violin(draw_quantiles=TRUE, fill="cornflowerblue") + theme_bw() + 
  geom_boxplot(width=.1, outlier.colour=NA, alpha = .75) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) + ylab("Mean trait conservation") + 
  xlab("Category")

saveRDS(object=td_summaries, file="/home/gavin/projects/picrust2_manuscript/saved_RDS/trait_depth_summaries.rds")

setwd("/home/gavin/github_repos/picrust_repos/picrust2/testing_datasets/old_defaults")

# Read in original tables to compare with:
rrna_orig <- read.table(gzfile("16S_counts_mean_round_var.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
ec_orig <- read.table(gzfile("ec_level4_counts_mean_round_var.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
ko_orig <- read.table(gzfile("ko_counts_mean_round_var_subset.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
cog_orig <- read.table(gzfile("cog_counts_mean_round_var_subset.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
pfam_orig <- read.table(gzfile("pfam_counts_mean_round_var_subset.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
tigrfam_orig <- read.table(gzfile("tigrfam_counts_mean_round_var.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

characterize_db(rrna_orig)
characterize_db(cog_orig)
characterize_db(ec_orig)
characterize_db(ko_orig)
characterize_db(pfam_orig)
characterize_db(tigrfam_orig)

# Get overlapping genome ids between new and old sets.
overlapping_genomes <- rownames(rrna_orig)[which(rownames(rrna_orig) %in% rownames(rrna))]

# Remove all "cluster" genomes from this set.
overlapping_genomes <- overlapping_genomes[-grep("cluster", overlapping_genomes)]

# Double-check that rows and columns that overlap exactly between the new and old tables are identical
# (i.e. when the tables are subsetted to the same rows and columns).

rrna_subset <- rrna[overlapping_genomes,, drop=FALSE]
rrna_orig_subset <- rrna_orig[overlapping_genomes,, drop=FALSE]
identical(rrna_orig_subset, rrna_subset)

ec_subset <- ec[overlapping_genomes, ]
ec_orig_subset <- ec_orig[overlapping_genomes, colnames(ec)]
identical(ec_orig_subset, ec_subset)

cog_subset <- cog[overlapping_genomes, ]
cog_orig_subset <- cog_orig[overlapping_genomes, colnames(cog)]
identical(cog_orig_subset, cog_subset)

ko_subset <- ko[overlapping_genomes, ]
ko_orig_subset <- ko_orig[overlapping_genomes, colnames(ko)]
identical(ko_orig_subset, ko_subset)

pfam_subset <- pfam[overlapping_genomes, ]
pfam_orig_subset <- pfam_orig[overlapping_genomes, colnames(pfam)]
identical(pfam_orig_subset, pfam_subset)

tigrfam_subset <- tigrfam[overlapping_genomes, ]
tigrfam_orig_subset <- tigrfam_orig[overlapping_genomes, colnames(tigrfam)]
identical(tigrfam_orig_subset, tigrfam_subset)

# Calculate trait depth for original tables to compare with new.
old_ref_tree <- read.tree("img_centroid_16S_aligned.tree")

rrna_orig_td <- calc_trait_depth(old_ref_tree, rrna_orig, 1)
cog_orig_td <- calc_trait_depth(old_ref_tree, cog_orig, 40)
ec_orig_td <- calc_trait_depth(old_ref_tree, ec_orig, 40)
ko_orig_td <- calc_trait_depth(old_ref_tree, ko_orig, 40)
pfam_orig_td <- calc_trait_depth(old_ref_tree, pfam_orig, 40)
tigrfam_orig_td <- calc_trait_depth(old_ref_tree, tigrfam_orig, 40)

# Parse out mean and sd in trait depth for each dataset.
rrna_orig_td_summary <- data.frame(t(data.frame(lapply(rrna_orig_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
cog_orig_td_summary <- data.frame(t(data.frame(lapply(cog_orig_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
ec_orig_td_summary <- data.frame(t(data.frame(lapply(ec_orig_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
ko_orig_td_summary <- data.frame(t(data.frame(lapply(ko_orig_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
pfam_orig_td_summary <- data.frame(t(data.frame(lapply(pfam_orig_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))
tigrfam_orig_td_summary <- data.frame(t(data.frame(lapply(tigrfam_orig_td, function(x){ return(c(x$mean_depth, sqrt(x$var_depth))) }))))

# Add new column for what category the functions are in.
rrna_orig_td_summary$category <- "16S"
cog_orig_td_summary$category <- "COG"
ec_orig_td_summary$category <- "E.C. Number"
ko_orig_td_summary$category <- "KO"
pfam_orig_td_summary$category <- "Pfam"
tigrfam_orig_td_summary$category <- "TIGRFAM"

# Concatenate the rows into 1 final table.
td_orig_summaries <- rbind(rrna_orig_td_summary, cog_orig_td_summary, ec_orig_td_summary,
                      ko_orig_td_summary, pfam_orig_td_summary, tigrfam_orig_td_summary)

colnames(td_orig_summaries) <- c("mean", "sd", "category")

# Set factors to be in order from most conserved to least.
td_orig_summaries$category <- factor(td_orig_summaries$category, levels=c("16S", "COG", "TIGRFAM", "E.C. Number", "Pfam", "KO"))

dodge <- position_dodge(width = 1)

ggplot(td_orig_summaries, aes(category, mean)) + geom_violin(draw_quantiles=TRUE, fill="cornflowerblue") + theme_bw() + 
  geom_boxplot(width=.1, outlier.colour=NA, alpha = 0.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.5)) + ylab("Mean trait conservation") + 
  xlab("Category")
