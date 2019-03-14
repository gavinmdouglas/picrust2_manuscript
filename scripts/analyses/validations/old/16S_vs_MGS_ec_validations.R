### Comparing predicted E.C. numbers based on 16S to MGS "gold standard".
### Also compared to Paprica, which output E.C. number predictions.
### Calculated cosine similarities and spearman correlations between MGS and 16S data and saved as RDS object.
### Plotting also done in this script as well.

library(ggplot2)
library(cowplot)

setwd("/home/gavin/projects/picrust2_manuscript/data/")
source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

### Read in database files (used for calculating null distribution).
ec <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ec.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

### Read in expected metacyc database (based on reference E.C. number database).

### Read in HMP E.C. tables.
hmp_picrust2_ec <- read.table("16S_validation/picrust2_out/hmp_picrust2_ec.tsv",
                               header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
hmp_paprica_ec <- data.frame(t(read.table("16S_validation/paprica_out/hmp_paprica_ec.csv",
                             header=T, sep=",", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="", check.names=FALSE)))
colnames(hmp_paprica_ec) <- gsub("_rerep_seqs.", "", colnames(hmp_paprica_ec))
rownames(hmp_paprica_ec) <- paste("EC:", rownames(hmp_paprica_ec), sep="")

# Remove paprica rows that are all 0s.
hmp_paprica_ec <- hmp_paprica_ec[-which(rowSums(hmp_paprica_ec) == 0),]

hmp_mgs_ec <- read.table("mgs_validation/hmp/humann2_ec_unstrat.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

# Subset samples to just those overlapping between PICRUSt2, PAPRICA, and HUMAnN2 output tables.
hmp_ec_cols2keep <- colnames(hmp_paprica_ec)[which(colnames(hmp_paprica_ec) %in% colnames(hmp_picrust2_ec))]
hmp_ec_cols2keep <- hmp_ec_cols2keep[which(hmp_ec_cols2keep %in% colnames(hmp_mgs_ec))]
hmp_mgs_ec <- hmp_mgs_ec[,hmp_ec_cols2keep]
hmp_paprica_ec <- hmp_paprica_ec[,hmp_ec_cols2keep]
hmp_picrust2_ec <- hmp_picrust2_ec[,hmp_ec_cols2keep]

### Read in mammal E.C. tables.
mammal_picrust2_ec <- read.table("16S_validation/picrust2_out/mammal_picrust2_ec.tsv",
                              header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
mammal_paprica_ec <- data.frame(t(read.table("16S_validation/paprica_out/mammal_paprica_ec.csv",
                                          header=T, sep=",", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="", check.names=FALSE)))
colnames(mammal_paprica_ec) <- gsub("_rerep_seqs.", "", colnames(mammal_paprica_ec))
rownames(mammal_paprica_ec) <- paste("EC:", rownames(mammal_paprica_ec), sep="")

# Remove paprica rows that are all 0s.
mammal_paprica_ec <- mammal_paprica_ec[-which(rowSums(mammal_paprica_ec) == 0),]

mammal_mgs_ec <- read.table("mgs_validation/mammalian_stool/humann2_ec_unstrat.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

# Subset samples to just those overlapping between PICRUSt2, PAPRICA, and HUMAnN2 output tables.
mammal_ec_cols2keep <- colnames(mammal_paprica_ec)[which(colnames(mammal_paprica_ec) %in% colnames(mammal_picrust2_ec))]
mammal_ec_cols2keep <- mammal_ec_cols2keep[which(mammal_ec_cols2keep %in% colnames(mammal_mgs_ec))]
mammal_mgs_ec <- mammal_mgs_ec[,mammal_ec_cols2keep]
mammal_paprica_ec <- mammal_paprica_ec[,mammal_ec_cols2keep]
mammal_picrust2_ec <- mammal_picrust2_ec[,mammal_ec_cols2keep]

### Read in ocean E.C. tables.
ocean_picrust2_ec <- read.table("16S_validation/picrust2_out/ocean_picrust2_ec.tsv",
                                 header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
ocean_paprica_ec <- data.frame(t(read.table("16S_validation/paprica_out/ocean_paprica_ec.csv",
                                             header=T, sep=",", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="", check.names=FALSE)))
colnames(ocean_paprica_ec) <- gsub("_rerep_seqs.", "", colnames(ocean_paprica_ec))
rownames(ocean_paprica_ec) <- paste("EC:", rownames(ocean_paprica_ec), sep="")

# Remove paprica rows that are all 0s.
ocean_paprica_ec <- ocean_paprica_ec[-which(rowSums(ocean_paprica_ec) == 0),]

ocean_mgs_ec <- read.table("mgs_validation/ocean/humann2_ec_unstrat.tsv",
                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

# Subset samples to just those overlapping between PICRUSt2, PAPRICA, and HUMAnN2 output tables.
ocean_ec_cols2keep <- colnames(ocean_paprica_ec)[which(colnames(ocean_paprica_ec) %in% colnames(ocean_picrust2_ec))]
ocean_ec_cols2keep <- ocean_ec_cols2keep[which(ocean_ec_cols2keep %in% colnames(ocean_mgs_ec))]
ocean_mgs_ec <- ocean_mgs_ec[,ocean_ec_cols2keep]
ocean_paprica_ec <- ocean_paprica_ec[,ocean_ec_cols2keep]
ocean_picrust2_ec <- ocean_picrust2_ec[,ocean_ec_cols2keep]

### Read in ocean E.C. tables.
soil_picrust2_ec <- read.table("16S_validation/picrust2_out/soil_picrust2_ec.tsv",
                                header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
soil_paprica_ec <- data.frame(t(read.table("16S_validation/paprica_out/soil_paprica_ec.csv",
                                            header=T, sep=",", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="", check.names=FALSE)))
colnames(soil_paprica_ec) <- gsub("_rerep_seqs.", "", colnames(soil_paprica_ec))
rownames(soil_paprica_ec) <- paste("EC:", rownames(soil_paprica_ec), sep="")

# Remove paprica rows that are all 0s.
soil_paprica_ec <- soil_paprica_ec[-which(rowSums(soil_paprica_ec) == 0),]

soil_mgs_ec <- read.table("mgs_validation/soil_crossbiome/humann2_ec_unstrat.tsv",
                           header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

# Subset samples to just those overlapping between PICRUSt2, PAPRICA, and HUMAnN2 output tables.
soil_ec_cols2keep <- colnames(soil_paprica_ec)[which(colnames(soil_paprica_ec) %in% colnames(soil_picrust2_ec))]
soil_ec_cols2keep <- soil_ec_cols2keep[which(soil_ec_cols2keep %in% colnames(soil_mgs_ec))]
soil_mgs_ec <- soil_mgs_ec[,soil_ec_cols2keep]
soil_paprica_ec <- soil_paprica_ec[,soil_ec_cols2keep]
soil_picrust2_ec <- soil_picrust2_ec[,soil_ec_cols2keep]

# Calculate cosine and spearman correlations.
hmp_ec_mgs_null_cosine <- rand_sample_vs_func_table(db = ec, tab = hmp_mgs_ec, metric="cosine")
hmp_ec_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_picrust2_ec, tab2 = hmp_mgs_ec, cat_string="PICRUSt2", metric="cosine")
hmp_ec_paprica_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_paprica_ec, tab2 = hmp_mgs_ec, cat_string="PAPRICA", metric="cosine")
hmp_ec_cosine_df <- rbind(hmp_ec_mgs_null_cosine, hmp_ec_picrust2_vs_mgs_cosine, hmp_ec_paprica_vs_mgs_cosine)

hmp_ec_mgs_null_spearman <- rand_sample_vs_func_table(db = ec, tab = hmp_mgs_ec, metric="spearman")
hmp_ec_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_picrust2_ec, tab2 = hmp_mgs_ec, cat_string="PICRUSt2", metric="spearman")
hmp_ec_paprica_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_paprica_ec, tab2 = hmp_mgs_ec, cat_string="PAPRICA", metric="spearman")
hmp_ec_spearman_df <- rbind(hmp_ec_mgs_null_spearman, hmp_ec_picrust2_vs_mgs_spearman, hmp_ec_paprica_vs_mgs_spearman)

mammal_ec_mgs_null_cosine <- rand_sample_vs_func_table(db = ec, tab = mammal_mgs_ec, metric="cosine")
mammal_ec_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_picrust2_ec, tab2 = mammal_mgs_ec, cat_string="PICRUSt2", metric="cosine")
mammal_ec_paprica_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_paprica_ec, tab2 = mammal_mgs_ec, cat_string="PAPRICA", metric="cosine")
mammal_ec_cosine_df <- rbind(mammal_ec_mgs_null_cosine, mammal_ec_picrust2_vs_mgs_cosine, mammal_ec_paprica_vs_mgs_cosine)

mammal_ec_mgs_null_spearman <- rand_sample_vs_func_table(db = ec, tab = mammal_mgs_ec, metric="spearman")
mammal_ec_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_picrust2_ec, tab2 = mammal_mgs_ec, cat_string="PICRUSt2", metric="spearman")
mammal_ec_paprica_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_paprica_ec, tab2 = mammal_mgs_ec, cat_string="PAPRICA", metric="spearman")
mammal_ec_spearman_df <- rbind(mammal_ec_mgs_null_spearman, mammal_ec_picrust2_vs_mgs_spearman, mammal_ec_paprica_vs_mgs_spearman)

ocean_ec_mgs_null_cosine <- rand_sample_vs_func_table(db = ec, tab = ocean_mgs_ec, metric="cosine")
ocean_ec_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_picrust2_ec, tab2 = ocean_mgs_ec, cat_string="PICRUSt2", metric="cosine")
ocean_ec_paprica_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_paprica_ec, tab2 = ocean_mgs_ec, cat_string="PAPRICA", metric="cosine")
ocean_ec_cosine_df <- rbind(ocean_ec_mgs_null_cosine, ocean_ec_picrust2_vs_mgs_cosine, ocean_ec_paprica_vs_mgs_cosine)

ocean_ec_mgs_null_spearman <- rand_sample_vs_func_table(db = ec, tab = ocean_mgs_ec, metric="spearman")
ocean_ec_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_picrust2_ec, tab2 = ocean_mgs_ec, cat_string="PICRUSt2", metric="spearman")
ocean_ec_paprica_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_paprica_ec, tab2 = ocean_mgs_ec, cat_string="PAPRICA", metric="spearman")
ocean_ec_spearman_df <- rbind(ocean_ec_mgs_null_spearman, ocean_ec_picrust2_vs_mgs_spearman, ocean_ec_paprica_vs_mgs_spearman)

soil_ec_mgs_null_cosine <- rand_sample_vs_func_table(db = ec, tab = soil_mgs_ec, metric="cosine")
soil_ec_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = soil_picrust2_ec, tab2 = soil_mgs_ec, cat_string="PICRUSt2", metric="cosine")
soil_ec_paprica_vs_mgs_cosine <- cor_all_cols(tab1 = soil_paprica_ec, tab2 = soil_mgs_ec, cat_string="PAPRICA", metric="cosine")
soil_ec_cosine_df <- rbind(soil_ec_mgs_null_cosine, soil_ec_picrust2_vs_mgs_cosine, soil_ec_paprica_vs_mgs_cosine)

soil_ec_mgs_null_spearman <- rand_sample_vs_func_table(db = ec, tab = soil_mgs_ec, metric="spearman")
soil_ec_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = soil_picrust2_ec, tab2 = soil_mgs_ec, cat_string="PICRUSt2", metric="spearman")
soil_ec_paprica_vs_mgs_spearman <- cor_all_cols(tab1 = soil_paprica_ec, tab2 = soil_mgs_ec, cat_string="PAPRICA", metric="spearman")
soil_ec_spearman_df <- rbind(soil_ec_mgs_null_spearman, soil_ec_picrust2_vs_mgs_spearman, soil_ec_paprica_vs_mgs_spearman)

# Add dataset identified to each df.
hmp_ec_cosine_df$dataset <- "HMP"
hmp_ec_spearman_df$dataset <- "HMP"

mammal_ec_cosine_df$dataset <- "Mammal"
mammal_ec_spearman_df$dataset <- "Mammal"

ocean_ec_cosine_df$dataset <- "Ocean"
ocean_ec_spearman_df$dataset <- "Ocean"

soil_ec_cosine_df$dataset <- "Soil"
soil_ec_spearman_df$dataset <- "Soil"

# Combine all of these dataframes into a single one for cosine and spearman measures.
combined_ec_cosine_df <- do.call(rbind, list(hmp_ec_cosine_df, mammal_ec_cosine_df, ocean_ec_cosine_df, soil_ec_cosine_df))
combined_ec_spearman_df <- do.call(rbind, list(hmp_ec_spearman_df, mammal_ec_spearman_df, ocean_ec_spearman_df, soil_ec_spearman_df))

combined_ec_cosine_plot <- ggplot(combined_ec_cosine_df, aes(x=dataset, y=metric, fill=cat)) +
                                  geom_boxplot() +
                                  ylab("Cosine similarity") +
                                  xlab("Dataset") +
                                  ylim(low=0.25, high=1) +
                                  scale_fill_manual(values=c("grey", "#7CAE00", "#619CFF"))

combined_ec_spearman_plot <- ggplot(combined_ec_spearman_df, aes(x=dataset, y=metric, fill=cat)) +
  geom_boxplot() +
  ylab("Spearman correlation") +
  xlab("Dataset") +
  ylim(low=0.25, high=1) +
  scale_fill_manual(values=c("grey", "#7CAE00", "#619CFF"))

plot_grid(combined_ec_spearman_plot, combined_ec_cosine_plot, labels = c("A", "B"))

wilcox.test(hmp_ec_mgs_null_cosine$metric, hmp_ec_picrust2_vs_mgs_cosine$metric)
# W = 1357, p-value < 2.2e-16
wilcox.test(hmp_ec_picrust2_vs_mgs_cosine$metric, hmp_ec_paprica_vs_mgs_cosine$metric)
# W = 18282, p-value < 2.2e-16

wilcox.test(hmp_ec_mgs_null_spearman$metric, hmp_ec_picrust2_vs_mgs_spearman$metric)
# W = 270, p-value < 2.2e-16
wilcox.test(hmp_ec_picrust2_vs_mgs_spearman$metric, hmp_ec_paprica_vs_mgs_spearman$metric)
# W = 18624, p-value < 2.2e-16

wilcox.test(mammal_ec_mgs_null_cosine$metric, mammal_ec_picrust2_vs_mgs_cosine$metric)
# W = 0, p-value = 2.835e-06
wilcox.test(mammal_ec_picrust2_vs_mgs_cosine$metric, mammal_ec_paprica_vs_mgs_cosine$metric)
# W = 112, p-value = 0.000275
wilcox.test(mammal_ec_mgs_null_spearman$metric, mammal_ec_picrust2_vs_mgs_spearman$metric)
# W = 0, p-value = 2.835e-06
wilcox.test(mammal_ec_picrust2_vs_mgs_spearman$metric, mammal_ec_paprica_vs_mgs_spearman$metric)
# W = 121, p-value = 2.835e-06

wilcox.test(ocean_ec_mgs_null_cosine$metric, ocean_ec_picrust2_vs_mgs_cosine$metric)
# W = 0, p-value = 0.002165
wilcox.test(ocean_ec_picrust2_vs_mgs_cosine$metric, ocean_ec_paprica_vs_mgs_cosine$metric)
# W = 35, p-value = 0.004329
wilcox.test(ocean_ec_mgs_null_spearman$metric, ocean_ec_picrust2_vs_mgs_spearman$metric)
# W = 0, p-value = 0.002165
wilcox.test(ocean_ec_picrust2_vs_mgs_spearman$metric, ocean_ec_paprica_vs_mgs_spearman$metric)
# W = 36, p-value = 0.002165

wilcox.test(soil_ec_mgs_null_cosine$metric, soil_ec_picrust2_vs_mgs_cosine$metric)
# W = 78, p-value = 0.3761
wilcox.test(soil_ec_picrust2_vs_mgs_cosine$metric, soil_ec_paprica_vs_mgs_cosine$metric)
# W = 191, p-value = 9.472e-07
wilcox.test(soil_ec_mgs_null_spearman$metric, soil_ec_picrust2_vs_mgs_spearman$metric)
# W = 0, p-value = 4.985e-08
wilcox.test(soil_ec_picrust2_vs_mgs_spearman$metric, soil_ec_paprica_vs_mgs_spearman$metric)
# W = 92, p-value = 0.8036

# Save RDS objects:
saveRDS(object = combined_ec_spearman_df, file = "/home/gavin/projects/picrust2_manuscript/saved_RDS/combined_ec_spearman_df.rds")
saveRDS(object = combined_ec_cosine_df, file = "/home/gavin/projects/picrust2_manuscript/saved_RDS/combined_ec_cosine_df.rds")
