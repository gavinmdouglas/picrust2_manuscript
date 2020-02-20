### Test for features that significantly differ between CD and CN ileum samples with Aldex2.
### Summarize these significant hits and then test for significantly enriched functions in
### pathways of interest.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/count_tables/")

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

library(ALDEx2)
library(ggplot2)
library(stringi)
library(reshape2)
library(sinaplot)
library(cowplot)

# Get mean % contribution by a particular clade. Returns list for each function with % contributed by that clade.
df_clade_mean_contrib <- function(strat_func_by_tax_rel, clade_string, num_cores) {
  
  funcs <- gsub("\\|.+$", "", rownames(strat_func_by_tax_rel))
  
  if(length(which(duplicated(funcs))) > 0) {
    funcs <- funcs[-which(duplicated(funcs))]
  }
  
  func_mean_clade_contrib <- mclapply(funcs, function(x) { 
    single_func_mean_contrib(x, strat_func_by_tax_rel, clade_string)
    
  }, mc.cores=num_cores)
  
  names(func_mean_clade_contrib) <- funcs
  
  return(func_mean_clade_contrib)
  
}

single_func_mean_contrib <- function(func, func_df, clade_string) {
  
  func_rows <- func_df[grep(func, rownames(func_df)),]
  
  if(length(grep(clade_string, rownames(func_rows))) > 0) {
    
    clade_func_abun <- as.numeric(func_rows[grep(clade_string, rownames(func_rows)),])
    func_abun <- colSums(func_rows) + 0.00001
    
    return(mean(clade_func_abun / func_abun) * 100)
    
  } else {
    
    return(0)
    
  }
  
}

# Read in hmp2 metadata table.
hmp2_metadata <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_2018-08-20_col_subset.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE, check.names = TRUE)

# Read in stratified pathways (--per-seq-contrib option used).
pathabun_in_strat <- readRDS("hmp2_pathabun_strat_filt_count_Ileum.rds")

# Transpose stratified table and re-separate pathways and sequences as separate columns.
pathabun_in_strat <- data.frame(t(pathabun_in_strat), check.names = FALSE)
orig_col <- colnames(pathabun_in_strat)
pathabun_in_strat$pathway <- as.character(sapply(rownames(pathabun_in_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][1] } ))
pathabun_in_strat$sequence <- as.character(sapply(rownames(pathabun_in_strat), function(x) { stri_split(str=x, regex="\\|")[[1]][2] } ))
pathabun_in_strat <- pathabun_in_strat[, c("pathway", "sequence", orig_col)]

# Unstrat pathways and phenotypes.
pathabun_in_unstrat <- readRDS("hmp2_pathabun_filt_count_Ileum.rds")
pathabun_in_unstrat <- data.frame(t(pathabun_in_unstrat), check.names = FALSE)
pheno_in_unstrat <- readRDS("hmp2_pheno_count_Ileum.rds")

# Read in ASVs:
asv_abun <- readRDS("hmp2_biom_count_Ileum.rds")

descrip_gzfile <- gzfile('/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 'rt')

path_descrip <- read.table(descrip_gzfile, header=FALSE, sep="\t", row.names=1, comment.char="", quote="", stringsAsFactors = FALSE)
close(descrip_gzfile)

# Read in ASV taxonomy.
asv_tax_in <- read.table("../hmp2_asv_taxa.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)
asv_tax_in_levels <- add_tax_cols(asv_tax_in)
rownames(asv_tax_in_levels) <- asv_tax_in_levels$Feature.ID

# Get metadata subsetted to 16S samples.
biopsy_16S_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "biopsy_16S"), ]
biopsy_16S_meta <- biopsy_16S_meta[which(biopsy_16S_meta$biopsy_location == "Ileum"), ]
biopsy_16S_meta <- biopsy_16S_meta[-which(duplicated(biopsy_16S_meta$Participant.ID)), ]
rownames(biopsy_16S_meta) <- biopsy_16S_meta$Participant.ID

Ileum_CD_samples <- biopsy_16S_meta[which(biopsy_16S_meta$diagnosis == "CD"), "Participant.ID"]
Ileum_nonIBD_samples <- biopsy_16S_meta[which(biopsy_16S_meta$diagnosis == "nonIBD"), "Participant.ID"]
Ileum_CD_nonIBD_samples <- c(Ileum_CD_samples, Ileum_nonIBD_samples)

### Test for pathways and phenotypes that differ between nonIBD and CD individuals (at lenient cut-off).
Ileum_CD_non_IBD_pathabun_in <- pathabun_in_unstrat[, Ileum_CD_nonIBD_samples[which(Ileum_CD_nonIBD_samples %in% colnames(pathabun_in_unstrat))]]
Ileum_CD_non_IBD_pathabun_in_filt <- Ileum_CD_non_IBD_pathabun_in[which(rowSums(Ileum_CD_non_IBD_pathabun_in > 0) >= 0.33 * ncol(Ileum_CD_non_IBD_pathabun_in)), ]
Ileum_CD_non_IBD_pathabun_in_filt_aldex <- aldex(round(Ileum_CD_non_IBD_pathabun_in_filt), biopsy_16S_meta[colnames(Ileum_CD_non_IBD_pathabun_in_filt), "diagnosis"], effect=TRUE)
rownames(Ileum_CD_non_IBD_pathabun_in_filt_aldex)[which(Ileum_CD_non_IBD_pathabun_in_filt_aldex$wi.eBH < 0.2)]
# No significant pathways.

Ileum_CD_non_IBD_pheno_in <- pheno_in_unstrat[, Ileum_CD_nonIBD_samples[which(Ileum_CD_nonIBD_samples %in% colnames(pheno_in_unstrat))]]
Ileum_CD_non_IBD_pheno_in_filt <- Ileum_CD_non_IBD_pheno_in[which(rowSums(Ileum_CD_non_IBD_pheno_in > 0) >= 0.33 * ncol(Ileum_CD_non_IBD_pheno_in)), ]
Ileum_CD_non_IBD_pheno_in_filt_aldex <- aldex(round(Ileum_CD_non_IBD_pheno_in_filt), biopsy_16S_meta[colnames(Ileum_CD_non_IBD_pheno_in_filt), "diagnosis"], effect=TRUE)
rownames(Ileum_CD_non_IBD_pheno_in_filt_aldex)[which(Ileum_CD_non_IBD_pheno_in_filt_aldex$wi.eBH < 0.2)]
# No significant phenotypes.


### Identify ASVs and taxa that differ between control and CD individuals (at lenient cut-off).
Ileum_CD_non_IBD_asv_abun <- asv_abun[, Ileum_CD_nonIBD_samples[which(Ileum_CD_nonIBD_samples %in% colnames(pathabun_in_unstrat))]]
Ileum_CD_non_IBD_asv_abun_taxa <- abun_by_taxa_level(asv_tax_in_levels, Ileum_CD_non_IBD_asv_abun)

Ileum_CD_non_IBD_asv_all_levels <- rbind(Ileum_CD_non_IBD_asv_abun, Ileum_CD_non_IBD_asv_abun_taxa$Species,
                                         Ileum_CD_non_IBD_asv_abun_taxa$Genus, Ileum_CD_non_IBD_asv_abun_taxa$Family,
                                         Ileum_CD_non_IBD_asv_abun_taxa$Order, Ileum_CD_non_IBD_asv_abun_taxa$Class, Ileum_CD_non_IBD_asv_abun_taxa$Phylum)

Ileum_CD_non_IBD_asv_all_levels_filt <-  Ileum_CD_non_IBD_asv_all_levels[which(rowSums(Ileum_CD_non_IBD_asv_all_levels > 0) >= 0.33 * ncol(Ileum_CD_non_IBD_asv_all_levels)), ]

Ileum_CD_non_IBD_asv_all_levels_filt_aldex <- aldex(Ileum_CD_non_IBD_asv_all_levels_filt, biopsy_16S_meta[colnames(Ileum_CD_non_IBD_pheno_in_filt), "diagnosis"], effect=TRUE)
sig_taxa <- rownames(Ileum_CD_non_IBD_asv_all_levels_filt_aldex)[which(Ileum_CD_non_IBD_asv_all_levels_filt_aldex$wi.eBH < 0.2)]

# "336454bed7f3f817495886a809d3b775"  --> k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__
# [2] "2031d34eae50b711bfb7c1a7b9a22f39" --> k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__                                                                                     
# [4] "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__"
# [5] "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira"     
# [6] "k__Bacteria; p__Proteobacteria"  

# Plot these significant hits in separate script.

# Based on significant taxa instead limit testing to functions contributed by taxa in the 2 groups:
# (1) Those that are at higher RA in nonIBD and (2) those are higher RA in CD

# Determine ASVs which fall within the significant taxa categories.
nonibd_sig_higher_ASVs <- c("2031d34eae50b711bfb7c1a7b9a22f39", "336454bed7f3f817495886a809d3b775")
nonibd_sig_higher_ASVs <- c(nonibd_sig_higher_ASVs, rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Species == "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira; s__")])
nonibd_sig_higher_ASVs <- c(nonibd_sig_higher_ASVs, rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Genus == "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Lachnospira")])
# 34 ASVs

cd_sig_higher_ASVs <- rownames(asv_tax_in_levels)[which(asv_tax_in_levels$Phylum == "k__Bacteria; p__Proteobacteria")]
# 192 ASVs

nonibd_sig_higher_out <- compare_func_ratios_by_asv_groups(strat_table = pathabun_in_strat,
                                                           ASVs_of_interest = nonibd_sig_higher_ASVs,
                                                           sample_group1 = Ileum_CD_samples[which(Ileum_CD_samples %in% colnames(pathabun_in_strat))],
                                                           sample_group2 = Ileum_nonIBD_samples[which(Ileum_nonIBD_samples %in% colnames(pathabun_in_strat))])

cd_sig_higher_out <- compare_func_ratios_by_asv_groups(strat_table = pathabun_in_strat,
                                                       ASVs_of_interest = cd_sig_higher_ASVs,
                                                       sample_group1 = Ileum_CD_samples[which(Ileum_CD_samples %in% colnames(pathabun_in_strat))],
                                                       sample_group2 = Ileum_nonIBD_samples[which(Ileum_nonIBD_samples %in% colnames(pathabun_in_strat))])


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
cd_sig_higher_ratio_prep$diagnosis <- biopsy_16S_meta[rownames(cd_sig_higher_ratio_prep), "diagnosis"]
cd_sig_higher_ratio_prep_melt <- melt(cd_sig_higher_ratio_prep)
cd_sig_higher_ratio_prep_melt$descrip <- paste(cd_sig_higher_ratio_prep_melt$variable, path_descrip[as.character(cd_sig_higher_ratio_prep_melt$variable), "V2"], sep=": ")
cd_sig_higher_ratio_prep_melt$log2ratio <- log2(cd_sig_higher_ratio_prep_melt$value)

cd_sig_higher_ratio_prep_melt$descrip <- factor(cd_sig_higher_ratio_prep_melt$descrip,
                                                levels=c("PWY-5189: tetrapyrrole biosynthesis II (from glycine)",
                                                         "PWY-5188: tetrapyrrole biosynthesis I (from glutamate)",
                                                         "PWY1G-0: mycothiol biosynthesis"))

cd_sig_higher_ratio_plot <- ggplot(cd_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=diagnosis)) +
                            geom_boxplot(width=0.75, outlier.shape = NA) +
                            coord_flip() +
                            scale_fill_manual(values=c("black", "grey")) +
                            xlab("Pathway") +
                            ylab("log2((Contributed by Proteobacteria + 1)/(Contributed by other + 1))") +
                            labs(fill="Disease State") +
                            theme(legend.position = c(0.8, 0.2),
                                  legend.background = element_rect(color = "black", 
                                                                   fill = "white", size = 0.2, linetype = "solid"))

# Figure out ranking of significant features based on absolute mean differences.
nonibd_sig_higher_out_wil$renamed_feat <- gsub("nonIBD-higher_", "", nonibd_sig_higher_out_wil$feature)
nonibd_sig_higher_out_wil_sig <- nonibd_sig_higher_out_wil[which(nonibd_sig_higher_out_wil$renamed_feat %in% combined_out_wil_sig_nonIBD_higher),]
nonibd_sig_higher_out_wil_sig_mean_diff_order <- order(nonibd_sig_higher_out_wil_sig$mean_diff)
nonibd_sig_higher_out_wil_sig$full_descrip <- paste(nonibd_sig_higher_out_wil_sig$renamed_feat, path_descrip[as.character(nonibd_sig_higher_out_wil_sig$renamed_feat), "V2"], sep=": ")

nonibd_sig_higher_ratio_prep <- data.frame(t(nonibd_sig_higher_out$ratio[combined_out_wil_sig_nonIBD_higher,]), check.names=FALSE)
nonibd_sig_higher_ratio_prep$sample <- rownames(nonibd_sig_higher_ratio_prep)
nonibd_sig_higher_ratio_prep$diagnosis <- ASV_samples_meta[rownames(nonibd_sig_higher_ratio_prep), "diagnosis"]
nonibd_sig_higher_ratio_prep_melt <- melt(nonibd_sig_higher_ratio_prep)
nonibd_sig_higher_ratio_prep_melt$descrip <- paste(nonibd_sig_higher_ratio_prep_melt$variable, path_descrip[as.character(nonibd_sig_higher_ratio_prep_melt$variable), "V2"], sep=": ")
nonibd_sig_higher_ratio_prep_melt$descrip <- factor(nonibd_sig_higher_ratio_prep_melt$descrip, levels=c(nonibd_sig_higher_out_wil_sig$full_descrip[nonibd_sig_higher_out_wil_sig_mean_diff_order]))
nonibd_sig_higher_ratio_prep_melt$log2ratio <- log2(nonibd_sig_higher_ratio_prep_melt$value)

nonibd_sig_higher_ratio_plot <- ggplot(nonibd_sig_higher_ratio_prep_melt, aes(x=descrip, y=log2ratio, fill=diagnosis)) + geom_boxplot(outlier.shape=NA) +
  coord_flip() + scale_fill_manual(values=c("black", "grey")) + xlab("Pathway") + ylab("log2((Contributed by significant Clostridiales)/(Contributed by other))")


# Write out RDS.
saveRDS(object = nonibd_sig_higher_ratio_prep_melt, file = "../results_out/nonibd_sig_higher_ratio_prep_melt.rds")
saveRDS(object = cd_sig_higher_ratio_prep_melt, file = "../results_out/cd_sig_higher_ratio_prep_melt.rds")
