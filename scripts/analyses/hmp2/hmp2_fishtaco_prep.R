### Prep PICRUSt2 output files for input to FishTaco.
### Subset tables to only those samples in subset of ileum CD vs control samples.

### The 4 input files required for FishTaco:
# Taxa	Sample1	Sample2	Sample3
# Taxon1	0.2	0.5	0.1
# Taxon2	0.3	0.1	0.7
 
# Function	Sample1	Sample2	Sample3
# K00001	0.1	0.1	0.1
# K00002	0.3	0.1	0.1

# Sample	Label
# Sample1	1
# Sample2	0

# Taxa	K00001	K00002	K00003	K00004	K00005	K00006
# Taxon1	1	2	0	1	2	0
# Taxon2	0	3	0	1	0	0

rm(list=ls(all=TRUE))

# Function to write 
write_table_w_rows_as_col <- function(outfile, df2write, row_colname) {
 
  orig_col <- colnames(df2write)
   
  df2write[, row_colname] <- rownames(df2write)
  
  df2write <- df2write[, c(row_colname, orig_col)]

  write.table(x = df2write, file = outfile, quote = FALSE, col.names = TRUE, row.names=FALSE, sep="\t")
    
}

setwd("/home/gavin/projects/hmp2_ibd_working/16S/April2018_redownload/picrust2_full_output_2.1.0-b")

# Read in all input files.

# Mapfile.
hmp2_map <- read.table("../hmp2_ibd_16S_map.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char="")

# ASV table.
hmp2_asv_abun <- read.table("../deblur_output_final/hmp2_ibd_16S.biom.tsv", header=TRUE, skip=1, comment.char="", stringsAsFactors = FALSE, sep="\t", row.names=1)

# Phenotype predicitons.
pheno_metagenome <- read.table("metagenome_out_pheno/pred_metagenome_unstrat.tsv", header=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, row.names=1)
pheno_gc <- read.table("pheno_predicted.tsv", header=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, row.names=1)

# pathway predicitons.
pathway_metagenome <- read.table("pathways_per_seq_by_hand/metagenome_out/pred_metagenome_unstrat.tsv", header=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, row.names=1)
pathway_gc <- read.table("pathways_per_seq_by_hand/asv_pathway_gc.tsv", header=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, row.names=1)

# First subet mapfile to samples of interest.
hmp2_map_ileum <- hmp2_map[which(hmp2_map$biopsy_location == "Ileum"), ]
hmp2_map_ileum <- hmp2_map_ileum[-which(duplicated(hmp2_map_ileum$Participant.ID)), ]
hmp2_map_ileum_noUC <- hmp2_map_ileum[-which(hmp2_map_ileum$diagnosis == "UC"), ]
hmp2_map_ileum_noUC <- hmp2_map_ileum_noUC[which(hmp2_map_ileum_noUC$X.SampleID %in% colnames(pheno_metagenome)), ]

# Subset ASV table to samples of interest, remove any ASVs that are all 0, convert to relabun, and write out.
hmp2_asv_abun_subset <- hmp2_asv_abun[ , hmp2_map_ileum_noUC$X.SampleID]
hmp2_asv_abun_subset <- hmp2_asv_abun_subset[-which(rowSums(hmp2_asv_abun_subset) == 0), ]

hmp2_asv_abun_subset_relab <- data.frame(sweep(hmp2_asv_abun_subset, 2, colSums(hmp2_asv_abun_subset), '/'), check.names=FALSE)
write_table_w_rows_as_col(outfile = "../fishtaco_input/hmp2_asv_ileum_CD_CN.tsv", df2write = hmp2_asv_abun_subset_relab, row_colname = "Taxa")


# Subset metagenome abundance tables, convert to relative abundance, and write.
pheno_metagenome_subset <- pheno_metagenome[, hmp2_map_ileum_noUC$X.SampleID]
pathway_metagenome_subset <- pathway_metagenome[, hmp2_map_ileum_noUC$X.SampleID]
pathway_metagenome_subset <- pathway_metagenome_subset[-which(rowSums(pathway_metagenome_subset) == 0), ]

pheno_metagenome_subset_relab <- data.frame(sweep(pheno_metagenome_subset, 2, colSums(pheno_metagenome_subset), '/'), check.names=FALSE)
pathway_metagenome_subset_relab <- data.frame(sweep(pathway_metagenome_subset, 2, colSums(pathway_metagenome_subset), '/'), check.names=FALSE)

write_table_w_rows_as_col(outfile = "../fishtaco_input/hmp2_pheno_ileum_CD_CN.tsv", df2write = pheno_metagenome_subset_relab, row_colname = "Function")
write_table_w_rows_as_col(outfile = "../fishtaco_input/hmp2_pathway_ileum_CD_CN.tsv", df2write = pathway_metagenome_subset_relab, row_colname = "Function")

# Write out genome content tables with subset of ASVs remaining in ASV abundance table and write out with first column "Taxa".
write_table_w_rows_as_col(outfile = "../fishtaco_input/hmp2_pheno_gc.tsv", df2write = pheno_gc[rownames(hmp2_asv_abun_subset), ], row_colname = "Taxa")
write_table_w_rows_as_col(outfile = "../fishtaco_input/hmp2_pathway_gc.tsv", df2write = pathway_gc[rownames(hmp2_asv_abun_subset), ], row_colname = "Taxa")

# Make sample label file and write out.
ileum_CD_CN_labels <- data.frame(Sample=hmp2_map_ileum_noUC$X.SampleID, Label=NA)
ileum_CD_CN_labels$Label[which(hmp2_map_ileum_noUC$diagnosis == "nonIBD")] <- 0
ileum_CD_CN_labels$Label[which(hmp2_map_ileum_noUC$diagnosis == "CD")] <- 1
write.table(file = "../fishtaco_input/labels_hmp2_ileum_CD_CN.tsv", x = ileum_CD_CN_labels, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
