### Commands to analyze IBD Multi'omics database data as an example of PICRUSt2 ###

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/ibdmdb_data/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

# To install ALDEx2:
# source("https://bioconductor.org/biocLite.R")
# biocLite("ALDEx2")

# Load required packages.
library(ALDEx2)
library(ggplot2)
library(reshape2)
library(cowplot)

# Read in unstratified and stratified pathway abundances.
pathabun_in <- read.table(gzfile("picrust2_pathways/pathways_out/path_abun_unstrat.tsv.gz"),
                          header=T, sep="\t", stringsAsFactors = FALSE,
                          comment.char="", quote="", row.names=1)

pathabun_in_strat_per_seq <- read.table(gzfile("picrust2_pathways/pathways_out_per_seq/path_abun_strat.tsv.gz"),
                                        header=T, sep="\t", stringsAsFactors = FALSE,
                                        comment.char="", quote="")

# Read in table of pathway descriptions (part of PICRUSt2 github repo).
path_descrip <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/description_mapfiles/metacyc_pathways_info_prokaryotes.txt.gz"),
                           header=F, sep="\t", stringsAsFactors = FALSE, quote="", comment.char="")
rownames(path_descrip) <- path_descrip$V1

# Remove unnecessary columns from unstratified pathway table.
col2remove <- c("function", "pathway", "description")
pathabun_in <- pathabun_in[, -which(colnames(pathabun_in) %in% col2remove)]

# Read in mapfile.
hmp2_map <- read.table("hmp2_metadata_16S_reformat.tsv", sep="\t", header=T, comment.char="",
                       quote="", stringsAsFactors = FALSE)

# Add "HMP.2." to the start of all external ids in metadata table.
hmp2_map$External.ID <- paste(rep_len("HMP.2.", length(hmp2_map$External.ID)), hmp2_map$External.ID, sep="")

# Remove 1 row of duplicated ids.
hmp2_map <- hmp2_map[-which(duplicated(hmp2_map$External.ID)),]

# Set rownames to be external ids.
rownames(hmp2_map) <- hmp2_map$External.ID

# Read in ASV BIOM table and remove samples with fewer than 4000 reads.
asv_in <- read.table("hmp2_ibd_16S_unfilt.biom.tsv", header=T, skip = 1, comment.char="", sep="\t", row.names=1)
asv_in <- asv_in[, -which(colSums(asv_in) < 4000)]

# Read in assigned taxonomy of each ASV.
asv_tax_in <- read.table("hmp2_ibd_16S_taxonomy_unfilt.tsv", header=T, sep="\t", stringsAsFactors = FALSE)

# Make new df with columns added for each taxonomic level.
asv_tax_in_levels <- add_tax_cols(asv_tax_in)
rownames(asv_tax_in_levels) <- asv_tax_in_levels$Feature.ID

# Get overlapping samples between BIOM, MAP, and predicted pathways output.
overlapping_samples <- colnames(asv_in)[which(colnames(asv_in) %in% rownames(hmp2_map))]
overlapping_samples <- overlapping_samples[which(overlapping_samples %in% colnames(pathabun_in))]

# Only keep overlapping samples.
asv_in <- asv_in[, overlapping_samples]
pathabun_in <- pathabun_in[,overlapping_samples]
pathabun_in_strat_per_seq <- pathabun_in_strat_per_seq[,c("pathway", "sequence", overlapping_samples)]
hmp2_map <- hmp2_map[overlapping_samples,]

# Susbet to ileum samples and CD vs control patients only.
nonibd_cd_ileum <- hmp2_map[which(hmp2_map$biopsy_location == "Ileum" & hmp2_map$diagnosis != "UC"), "External.ID"]
map_nonibd_cd_ileum <- hmp2_map[nonibd_cd_ileum,]

# Remove duplicated patient ids:
map_nonibd_cd_ileum <- map_nonibd_cd_ileum[! duplicated(map_nonibd_cd_ileum$Participant.ID),]

# Vectors of sample names only.
ileum_nonibd <- rownames(map_nonibd_cd_ileum)[which(map_nonibd_cd_ileum$diagnosis=="nonIBD")]
ileum_cd <- rownames(map_nonibd_cd_ileum)[which(map_nonibd_cd_ileum$diagnosis=="CD")]

### Identify at pathways that differ between control and CD individuals (at lenient cut-off).
cd_ileum_pathabun_in <- pathabun_in[, rownames(map_nonibd_cd_ileum)]

# Remove pathways not found at a depth of 10 in at least 25% of samples.
cd_ileum_pathabun_in_filt <- cd_ileum_pathabun_in[which(rowSums(cd_ileum_pathabun_in > 10) >= 0.25*ncol(cd_ileum_pathabun_in)),]

# Round pathway abundances, run ALDEX2 between CD and controls.
cd_ileum_pathabun_in_filt_aldex <- aldex(round(cd_ileum_pathabun_in_filt), map_nonibd_cd_ileum$diagnosis, effect=TRUE)

# Identify significant pathways.
rownames(cd_ileum_pathabun_in_filt_aldex)[which(cd_ileum_pathabun_in_filt_aldex$wi.eBH < 0.25)]
#"POLYAMINSYN3-PWY

### Identify ASVs and taxa that differ between control and CD individuals (again at lenient cut-off).
cd_ileum_asv_in <- asv_in[, rownames(map_nonibd_cd_ileum)]

# Get tables of abundances broken down by each taxonomic level.
cd_ileum_asv_in_taxa <- abun_by_taxa_level(asv_tax_in_levels, cd_ileum_asv_in)

# Combine all of these abundance tables into a single one.
cd_ileum_all_levels_in <- rbind(cd_ileum_asv_in, cd_ileum_asv_in_taxa$Species, cd_ileum_asv_in_taxa$Genus, cd_ileum_asv_in_taxa$Family,
                                cd_ileum_asv_in_taxa$Order, cd_ileum_asv_in_taxa$Class, cd_ileum_asv_in_taxa$Phylum)

# Remove all features that do not have an abundance of 10 in at least 25% of samples.
cd_ileum_all_levels_in_filt <-  cd_ileum_all_levels_in[which(rowSums(cd_ileum_all_levels_in > 10) >= 0.25*ncol(cd_ileum_all_levels_in)),]

# Run ALDEX2 as above.
cd_ileum_all_levels_in_filt_aldex <- aldex(cd_ileum_all_levels_in_filt, map_nonibd_cd_ileum$diagnosis, effect=TRUE)
sig_taxa <- rownames(cd_ileum_all_levels_in_filt_aldex)[which(cd_ileum_all_levels_in_filt_aldex$wi.eBH < 0.25)]
# [1] "2031d34eae50b711bfb7c1a7b9a22f39" --> k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__                                           
# [2] "336454bed7f3f817495886a809d3b775"  --> k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__                                   
# [3] "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__"
# [4] "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira"     
# [5] "k__Bacteria; p__Proteobacteria

# Plot the significant pathway and taxa abundances (after converting to relative abundance):

# Commands to convert each table to percentages.
cd_ileum_pathabun_in_relab <- data.frame(sweep(cd_ileum_pathabun_in, 2, colSums(cd_ileum_pathabun_in), FUN="/")) * 100
cd_ileum_asv_in_relab <- data.frame(sweep(cd_ileum_asv_in, 2, colSums(cd_ileum_asv_in), FUN="/")) * 100
cd_ileum_asv_in_Species_relab <- data.frame(sweep(cd_ileum_asv_in_taxa$Species, 2, colSums(cd_ileum_asv_in_taxa$Species), FUN="/")) * 100
cd_ileum_asv_in_Genus_relab <- data.frame(sweep(cd_ileum_asv_in_taxa$Genus, 2, colSums(cd_ileum_asv_in_taxa$Genus), FUN="/")) * 100
cd_ileum_asv_in_Phylum_relab <- data.frame(sweep(cd_ileum_asv_in_taxa$Phylum, 2, colSums(cd_ileum_asv_in_taxa$Phylum), FUN="/")) * 100

# Plot boxplots:
par(mfrow=c(2,3))
boxplot(as.numeric(cd_ileum_pathabun_in_relab["POLYAMINSYN3-PWY", ileum_nonibd]),
        as.numeric(cd_ileum_pathabun_in_relab["POLYAMINSYN3-PWY", ileum_cd]),
        main="POLYAMINSYN3-PWY\nsuperpathway of polyamine biosynthesis II", names=c("Non-IBD", "CD"), col="grey", ylab="Relative abundance (%)")

boxplot(as.numeric(cd_ileum_asv_in_relab["2031d34eae50b711bfb7c1a7b9a22f39", ileum_nonibd]),
        as.numeric(cd_ileum_asv_in_relab["2031d34eae50b711bfb7c1a7b9a22f39", ileum_cd]),
        main="ASV for unclassified o__Clostridiales", names=c("Non-IBD", "CD"), col="grey", ylab="Relative abundance (%)")

boxplot(as.numeric(cd_ileum_asv_in_relab["336454bed7f3f817495886a809d3b775", ileum_nonibd]),
        as.numeric(cd_ileum_asv_in_relab["336454bed7f3f817495886a809d3b775", ileum_cd]),
        main="ASV for unclassified g__Lachnospira", names=c("Non-IBD", "CD"), col="grey", ylab="Relative abundance (%)")

boxplot(as.numeric(cd_ileum_asv_in_Species_relab["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__", ileum_nonibd]),
        as.numeric(cd_ileum_asv_in_Species_relab["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__", ileum_cd]),
        main="g__Lachnospira unclassified species", names=c("Non-IBD", "CD"), col="grey", ylab="Relative abundance (%)")

boxplot(as.numeric(cd_ileum_asv_in_Genus_relab["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira", ileum_nonibd]),
        as.numeric(cd_ileum_asv_in_Genus_relab["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira", ileum_cd]),
        main="g__Lachnospira", names=c("Non-IBD", "CD"), col="grey", ylab="Relative abundance (%)")

boxplot(as.numeric(cd_ileum_asv_in_Phylum_relab["k__Bacteria; p__Proteobacteria", ileum_nonibd]),
        as.numeric(cd_ileum_asv_in_Phylum_relab["k__Bacteria; p__Proteobacteria", ileum_cd]),
        main="p__Proteobacteria", names=c("Non-IBD", "CD"), col="grey", ylab="Relative abundance (%)")

# Based on significant taxa instead limit testing to functions contributed by taxa in the 2 groups:
# (1) Those that are at higher RA in nonIBD and (2) those are higher RA in CD

# Determine ASVs which fall within the significant taxa categories (e.g. all ASVs within g__Lachnospira).
nonibd_sig_higher_ASVs <- c("2031d34eae50b711bfb7c1a7b9a22f39", "336454bed7f3f817495886a809d3b775")
nonibd_sig_higher_ASVs <- c(nonibd_sig_higher_ASVs, rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Species == "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__")])
nonibd_sig_higher_ASVs <- c(nonibd_sig_higher_ASVs, rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Genus == "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira")])
nonibd_sig_higher_ASVs <- unique(nonibd_sig_higher_ASVs)
# 17 ASVs

cd_sig_higher_ASVs <- rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Phylum == "k__Bacteria; p__Proteobacteria")]
# 218 ASVs

# Identify pathways that are differentially contributed between the ASVs of interest between the two sample groupings.
nonibd_sig_higher_out <- compare_func_ratios_by_asv_groups(strat_table = pathabun_in_strat_per_seq, ASVs_of_interest = nonibd_sig_higher_ASVs, sample_group1 = ileum_cd, sample_group2 = ileum_nonibd)

cd_sig_higher_out <- compare_func_ratios_by_asv_groups(strat_table = pathabun_in_strat_per_seq, ASVs_of_interest = cd_sig_higher_ASVs, sample_group1 = ileum_cd, sample_group2 = ileum_nonibd)


nonibd_sig_higher_out_wil <- nonibd_sig_higher_out$wilcox
nonibd_sig_higher_out_wil$feature <- paste("nonIBD-higher", nonibd_sig_higher_out_wil$feature, sep="_")

cd_sig_higher_out_wil <- cd_sig_higher_out$wilcox
cd_sig_higher_out_wil$feature <- paste("cd-higher", cd_sig_higher_out_wil$feature, sep="_")

combined_out_wil <- rbind(nonibd_sig_higher_out_wil, cd_sig_higher_out_wil)
combined_out_wil$fdr <- p.adjust(combined_out_wil$wilcox_p, "fdr")

combined_out_wil_sig <- combined_out_wil[which(combined_out_wil$fdr < 0.05),]

combined_out_wil_sig_CD_higher <- gsub("cd-higher_", "", combined_out_wil_sig[grep("cd-higher", combined_out_wil_sig$feature), "feature"])
combined_out_wil_sig_nonIBD_higher <- gsub("nonIBD-higher_" , "", combined_out_wil_sig[grep("nonIBD-higher", combined_out_wil_sig$feature), "feature"])

# Make plots of Proteobacteria/other sig pathways and Clostridia/other sig pathways
cd_sig_higher_ratio_prep <- data.frame(t(cd_sig_higher_out$ratio[combined_out_wil_sig_CD_higher,]), check.names=FALSE)
cd_sig_higher_ratio_prep$sample <- rownames(cd_sig_higher_ratio_prep)
cd_sig_higher_ratio_prep$Diagnosis <- map_nonibd_cd_ileum[rownames(cd_sig_higher_ratio_prep), "diagnosis"]
cd_sig_higher_ratio_prep$Diagnosis[which(cd_sig_higher_ratio_prep$Diagnosis== "nonIBD")] <- "Control"
cd_sig_higher_ratio_prep_melt <- melt(cd_sig_higher_ratio_prep)
cd_sig_higher_ratio_prep_melt$descrip <- factor(paste(cd_sig_higher_ratio_prep_melt$variable, path_descrip[as.character(cd_sig_higher_ratio_prep_melt$variable), "V2"], sep=": "))
cd_sig_higher_ratio_prep_melt$log2ratio <- log2(cd_sig_higher_ratio_prep_melt$value)

ggplot(cd_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=Diagnosis)) + geom_boxplot(width=0.75) +
   coord_flip() + scale_fill_manual(values=c("black", "grey")) + xlab("Pathway") + ylab(expression(log[2]~((Contributed~by~Proteobacteria)/(Contributed~by~Other)))) +
   scale_x_discrete(limits = rev(levels(cd_sig_higher_ratio_prep_melt$descrip)))


# Figure out ranking of significant features based on absolute mean differences for Clostridiales plot.
nonibd_sig_higher_out_wil$renamed_feat <- gsub("nonIBD-higher_", "", nonibd_sig_higher_out_wil$feature)
nonibd_sig_higher_out_wil_sig <- nonibd_sig_higher_out_wil[which(nonibd_sig_higher_out_wil$renamed_feat %in% combined_out_wil_sig_nonIBD_higher),]
nonibd_sig_higher_out_wil_sig_mean_diff_order <- order(nonibd_sig_higher_out_wil_sig$mean_diff)
nonibd_sig_higher_out_wil_sig$full_descrip <- paste(nonibd_sig_higher_out_wil_sig$renamed_feat, path_descrip[as.character(nonibd_sig_higher_out_wil_sig$renamed_feat), "V2"], sep=": ")

nonibd_sig_higher_ratio_prep <- data.frame(t(nonibd_sig_higher_out$ratio[combined_out_wil_sig_nonIBD_higher,]), check.names=FALSE)
nonibd_sig_higher_ratio_prep$sample <- rownames(nonibd_sig_higher_ratio_prep)
nonibd_sig_higher_ratio_prep$Diagnosis <- map_nonibd_cd_ileum[rownames(nonibd_sig_higher_ratio_prep), "diagnosis"]
nonibd_sig_higher_ratio_prep$Diagnosis[which(nonibd_sig_higher_ratio_prep$Diagnosis== "nonIBD")] <- "Control"
nonibd_sig_higher_ratio_prep_melt <- melt(nonibd_sig_higher_ratio_prep)
nonibd_sig_higher_ratio_prep_melt$descrip <- paste(nonibd_sig_higher_ratio_prep_melt$variable, path_descrip[as.character(nonibd_sig_higher_ratio_prep_melt$variable), "V2"], sep=": ")
nonibd_sig_higher_ratio_prep_melt$descrip <- factor(nonibd_sig_higher_ratio_prep_melt$descrip, levels=c(nonibd_sig_higher_out_wil_sig$full_descrip[nonibd_sig_higher_out_wil_sig_mean_diff_order]))
nonibd_sig_higher_ratio_prep_melt$log2ratio <- log2(nonibd_sig_higher_ratio_prep_melt$value)

ggplot(nonibd_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=Diagnosis)) + geom_boxplot(outlier.shape=NA) +
  coord_flip() + scale_fill_manual(values=c("black", "grey")) + xlab("Pathway") + ylab(expression(log[2]~((Contributed~by~Clostridiales)/(Contributed~by~Other))))
