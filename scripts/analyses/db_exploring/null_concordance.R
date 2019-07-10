### Determine baseline level of concordance between random subsets of each function database.

rm(list=ls(all=TRUE))

library(reshape2)

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")
setwd("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic")


# First read in latest functional databases:
rrna <- read.table(gzfile("16S.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
ec <- read.table(gzfile("ec.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
ko <- read.table(gzfile("ko.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
cog <- read.table(gzfile("cog.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
pfam <- read.table(gzfile("pfam.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
tigrfam <- read.table(gzfile("tigrfam.txt.gz"), header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)

ec_reps <- plot_func_compare_reps(func=ec, method="spearman")
ko_reps <- plot_func_compare_reps(func=ko, method="spearman")
cog_reps <- plot_func_compare_reps(func=cog, method="spearman")
pfam_reps <- plot_func_compare_reps(func=pfam, method="spearman")
tigrfam_reps <- plot_func_compare_reps(func=tigrfam, method="spearman")

cog_reps_melt <- melt(cog_reps)
cog_reps_melt$num_genomes <- gsub("n.", "", cog_reps_melt$variable)
cog_reps_melt$db <- "COG"

tigrfam_reps_melt <- melt(tigrfam_reps)
tigrfam_reps_melt$num_genomes <- gsub("n.", "", tigrfam_reps_melt$variable)
tigrfam_reps_melt$db <- "TIGRFAM"

ec_reps_melt <- melt(ec_reps)
ec_reps_melt$num_genomes <- gsub("n.", "", ec_reps_melt$variable)
ec_reps_melt$db <- "EC"

pfam_reps_melt <- melt(pfam_reps)
pfam_reps_melt$num_genomes <- gsub("n.", "", pfam_reps_melt$variable)
pfam_reps_melt$db <- "Pfam"

ko_reps_melt <- melt(ko_reps)
ko_reps_melt$num_genomes <- gsub("n.", "", ko_reps_melt$variable)
ko_reps_melt$db <- "KO"


combined_reps_melt <- rbind(cog_reps_melt, tigrfam_reps_melt, ec_reps_melt, pfam_reps_melt, ko_reps_melt)

combined_reps_melt$num_genomes <- factor(combined_reps_melt$num_genomes, levels=c(as.character(sort(unique(as.numeric(combined_reps_melt$num_genomes))))))

saveRDS(object = combined_reps_melt,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/null_concordance.rds")
