### Comparing predicted functions based on 16S to MGS "gold standard".
### Comparisons made between KEGG orthologs and MetaCyc pathways.
### For KEGG Orthologs the comparison with PICRUSt1 is also made.

setwd("/home/gavin/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

### Read in database files (used for calculating null distribution).
ko <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ko.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)

pathabun <- t(read.table(gzfile("/home/gavin/projects/picrust_pipeline/IMG_db/ref_genome_metacyc/path_abun_unstrat.tsv"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE))

pathcov <- t(read.table(gzfile("/home/gavin/projects/picrust_pipeline/IMG_db/ref_genome_metacyc/path_cov_unstrat.tsv"),
                         row.names=1, header=T, sep="\t", stringsAsFactors = FALSE))

### Read in expected metacyc database (based on reference E.C. number database).

### Read in hmp tables.
hmp_predicted_ko <- read.table("picrust2_out/hmp_picrust2_ko.tsv",
                               header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

hmp_predicted_ko_gg97 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_GG97/pred_metagenome_unstrat.tsv",
                               header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

hmp_predicted_ko_picrust1 <- read.table("picrust1_out/hmp_picrust1_ko.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")
colnames(hmp_predicted_ko_picrust1) <- gsub(".nonchimera.fasta", "", colnames(hmp_predicted_ko_picrust1))

hmp_predicted_ko_panfp <- read.table("panfp_out/hmp_panfp_ko.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

hmp_predicted_ko_piphillin <- read.table("piphillin_out/hmp_piphillin_ko.tsv",
                                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

hmp_predicted_ko_tax4fun <- read.table("tax4fun_out/hmp_tax4fun_ko.tsv",
                                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

hmp_mgs_ko <- read.table("../mgs_validation/hmp/humann2_ko_unstrat.tsv",
                         header=T, sep="\t", row.names=1)

hmp_predicted_pathabun <- read.table("picrust2_out/hmp_picrust2_pathabun.tsv",
                                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

hmp_predicted_pathcov <- read.table("picrust2_out/hmp_picrust2_pathcov.tsv",
                                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

hmp_mgs_pathabun <- read.table("../mgs_validation/hmp/humann2_pathabun_unstrat.tsv",
                               header=T, sep="\t", row.names=1)

hmp_mgs_pathcov <- read.table("../mgs_validation/hmp/humann2_pathcov_unstrat.tsv",
                              header=T, sep="\t", row.names=1)

# Only keep HMP samples that overlapped between both pipelines (due to some samples having low depth after limiting it to reference OTUs and rarefaction to 2000 reads).
hmp_overlapping_samples <- colnames(hmp_predicted_ko_picrust1)[which(colnames(hmp_predicted_ko_picrust1) %in% colnames(hmp_predicted_ko))]
hmp_overlapping_samples <- hmp_overlapping_samples[which(hmp_overlapping_samples %in% colnames(hmp_predicted_ko_tax4fun))]
hmp_overlapping_samples <- hmp_overlapping_samples[which(hmp_overlapping_samples %in% colnames(hmp_mgs_ko))]

# Subset to specific columns.
hmp_predicted_ko <- hmp_predicted_ko[,hmp_overlapping_samples]
hmp_predicted_ko_gg97 <- hmp_predicted_ko_gg97[,hmp_overlapping_samples]
hmp_predicted_ko_picrust1 <- hmp_predicted_ko_picrust1[,hmp_overlapping_samples]
hmp_predicted_ko_panfp <- hmp_predicted_ko_panfp[,hmp_overlapping_samples]
hmp_predicted_ko_piphillin <- hmp_predicted_ko_piphillin[,hmp_overlapping_samples]
hmp_predicted_ko_tax4fun <- hmp_predicted_ko_tax4fun[,hmp_overlapping_samples]
hmp_mgs_ko <- hmp_mgs_ko[,hmp_overlapping_samples]
hmp_predicted_pathabun <- hmp_predicted_pathabun[,hmp_overlapping_samples]
hmp_predicted_pathcov <- hmp_predicted_pathcov[,hmp_overlapping_samples]
hmp_mgs_pathabun <- hmp_mgs_pathabun[,hmp_overlapping_samples]
hmp_mgs_pathcov <- hmp_mgs_pathcov[,hmp_overlapping_samples]

### Read in mammal tables.
mammal_predicted_ko <- read.table("picrust2_out/mammal_picrust2_ko.tsv",
                               header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
mammal_predicted_ko_gg97 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_GG97/pred_metagenome_unstrat.tsv",
                                    header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_picrust1 <- read.table("picrust1_out/mammal_picrust1_ko.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")
colnames(mammal_predicted_ko_picrust1) <- gsub(".nonchimera.fasta", "", colnames(mammal_predicted_ko_picrust1))

mammal_predicted_ko_panfp <- read.table("panfp_out/mammal_panfp_ko.tsv",
                                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

mammal_predicted_ko_piphillin <- read.table("piphillin_out/mammal_piphillin_ko.tsv",
                                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

mammal_predicted_ko_tax4fun <- read.table("tax4fun_out/mammal_tax4fun_ko.tsv",
                                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

mammal_mgs_ko <- read.table("../mgs_validation/mammalian_stool/humann2_ko_unstrat.tsv",
                         header=T, sep="\t", row.names=1)

mammal_predicted_pathabun <- read.table("picrust2_out/mammal_picrust2_pathabun.tsv",
                                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_pathcov <- read.table("picrust2_out/mammal_picrust2_pathcov.tsv",
                                    header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_mgs_pathabun <- read.table("../mgs_validation/mammalian_stool/humann2_pathabun_unstrat.tsv",
                               header=T, sep="\t", row.names=1)

mammal_mgs_pathcov <- read.table("../mgs_validation/mammalian_stool/humann2_pathcov_unstrat.tsv",
                              header=T, sep="\t", row.names=1)

# Only keep mammal samples that overlapped between both pipelines (due to some samples having low depth after limiting it to reference OTUs and rarefaction to 2000 reads).
mammal_overlapping_samples <- colnames(mammal_predicted_ko_picrust1)[which(colnames(mammal_predicted_ko_picrust1) %in% colnames(mammal_predicted_ko))]
mammal_overlapping_samples <- mammal_overlapping_samples[which(mammal_overlapping_samples %in% colnames(mammal_predicted_ko_tax4fun))]
mammal_overlapping_samples <- mammal_overlapping_samples[which(mammal_overlapping_samples %in% colnames(mammal_mgs_ko))]

# Subset to specific columns.
mammal_predicted_ko <- mammal_predicted_ko[,mammal_overlapping_samples]
mammal_predicted_ko_gg97 <- mammal_predicted_ko_gg97[,mammal_overlapping_samples]
mammal_predicted_ko_picrust1 <- mammal_predicted_ko_picrust1[,mammal_overlapping_samples]
mammal_predicted_ko_panfp <- mammal_predicted_ko_panfp[,mammal_overlapping_samples]
mammal_predicted_ko_piphillin <- mammal_predicted_ko_piphillin[,mammal_overlapping_samples]
mammal_predicted_ko_tax4fun <- mammal_predicted_ko_tax4fun[,mammal_overlapping_samples]
mammal_mgs_ko <- mammal_mgs_ko[,mammal_overlapping_samples]
mammal_predicted_pathabun <- mammal_predicted_pathabun[,mammal_overlapping_samples]
mammal_predicted_pathcov <- mammal_predicted_pathcov[,mammal_overlapping_samples]
mammal_mgs_pathabun <- mammal_mgs_pathabun[,mammal_overlapping_samples]
mammal_mgs_pathcov <- mammal_mgs_pathcov[,mammal_overlapping_samples]


### Read in ocean tables.
ocean_predicted_ko <- read.table("picrust2_out/ocean_picrust2_ko.tsv",
                                  header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
ocean_predicted_ko_gg97 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_GG97/pred_metagenome_unstrat.tsv",
                                    header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

ocean_predicted_ko_picrust1 <- read.table("picrust1_out/ocean_picrust1_ko.tsv",
                                           header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")
colnames(ocean_predicted_ko_picrust1) <- gsub(".X.L001.R1.001.nonchimera.fasta", "", colnames(ocean_predicted_ko_picrust1))

ocean_predicted_ko_panfp <- read.table("panfp_out/ocean_panfp_ko.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

ocean_predicted_ko_piphillin <- read.table("piphillin_out/ocean_piphillin_ko.tsv",
                                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

ocean_predicted_ko_tax4fun <- read.table("tax4fun_out/ocean_tax4fun_ko.tsv",
                                          header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")
colnames(ocean_predicted_ko_tax4fun) <- gsub(".X.L001.R1.001", "", colnames(ocean_predicted_ko_tax4fun))

ocean_mgs_ko <- read.table("../mgs_validation/ocean/humann2_ko_unstrat.tsv",
                            header=T, sep="\t", row.names=1)

ocean_predicted_pathabun <- read.table("picrust2_out/ocean_picrust2_pathabun.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

ocean_predicted_pathcov <- read.table("picrust2_out/ocean_picrust2_pathcov.tsv",
                                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

ocean_mgs_pathabun <- read.table("../mgs_validation/ocean/humann2_pathabun_unstrat.tsv",
                                  header=T, sep="\t", row.names=1)

ocean_mgs_pathcov <- read.table("../mgs_validation/ocean/humann2_pathcov_unstrat.tsv",
                                 header=T, sep="\t", row.names=1)

# Only keep ocean samples that overlapped between both pipelines (due to some samples having low depth after limiting it to reference OTUs and rarefaction to 2000 reads).
ocean_overlapping_samples <- colnames(ocean_predicted_ko_picrust1)[which(colnames(ocean_predicted_ko_picrust1) %in% colnames(ocean_predicted_ko))]
ocean_overlapping_samples <- ocean_overlapping_samples[which(ocean_overlapping_samples %in% colnames(ocean_predicted_ko_tax4fun))]
ocean_overlapping_samples <- ocean_overlapping_samples[which(ocean_overlapping_samples %in% colnames(ocean_mgs_ko))]

# Subset to specific columns.
ocean_predicted_ko <- ocean_predicted_ko[,ocean_overlapping_samples]
ocean_predicted_ko_gg97 <- ocean_predicted_ko_gg97[,ocean_overlapping_samples]
ocean_predicted_ko_picrust1 <- ocean_predicted_ko_picrust1[,ocean_overlapping_samples]
ocean_predicted_ko_panfp <- ocean_predicted_ko_panfp[,ocean_overlapping_samples]
ocean_predicted_ko_piphillin <- ocean_predicted_ko_piphillin[,ocean_overlapping_samples]
ocean_predicted_ko_tax4fun <- ocean_predicted_ko_tax4fun[,ocean_overlapping_samples]
ocean_mgs_ko <- ocean_mgs_ko[,ocean_overlapping_samples]
ocean_predicted_pathabun <- ocean_predicted_pathabun[,ocean_overlapping_samples]
ocean_predicted_pathcov <- ocean_predicted_pathcov[,ocean_overlapping_samples]
ocean_mgs_pathabun <- ocean_mgs_pathabun[,ocean_overlapping_samples]
ocean_mgs_pathcov <- ocean_mgs_pathcov[,ocean_overlapping_samples]


### Read in soil tables.
soil_predicted_ko <- read.table("picrust2_out/soil_picrust2_ko.tsv",
                                  header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")
soil_predicted_ko_gg97 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_GG97/pred_metagenome_unstrat.tsv",
                                      header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

soil_predicted_ko_picrust1 <- read.table("picrust1_out/soil_picrust1_ko.tsv",
                                           header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, skip=1, comment.char="")
colnames(soil_predicted_ko_picrust1) <- gsub(".nonchimera.fna", "", colnames(soil_predicted_ko_picrust1))

soil_predicted_ko_panfp <- read.table("panfp_out/soil_panfp_ko.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

soil_predicted_ko_piphillin <- read.table("piphillin_out/soil_piphillin_ko.tsv",
                                            header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

soil_predicted_ko_tax4fun <- read.table("tax4fun_out/soil_tax4fun_ko.tsv",
                                          header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, comment.char="")

soil_mgs_ko <- read.table("../mgs_validation/soil_crossbiome/humann2_ko_unstrat.tsv",
                            header=T, sep="\t", row.names=1)

soil_predicted_pathabun <- read.table("picrust2_out/soil_picrust2_pathabun.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

soil_predicted_pathcov <- read.table("picrust2_out/soil_picrust2_pathcov.tsv",
                                       header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

soil_mgs_pathabun <- read.table("../mgs_validation/soil_crossbiome/humann2_pathabun_unstrat.tsv",
                                  header=T, sep="\t", row.names=1)

soil_mgs_pathcov <- read.table("../mgs_validation/soil_crossbiome/humann2_pathcov_unstrat.tsv",
                                 header=T, sep="\t", row.names=1)

# Only keep soil samples that overlapped between both pipelines (due to some samples having low depth after limiting it to reference OTUs and rarefaction to 2000 reads).
soil_overlapping_samples <- colnames(soil_predicted_ko_picrust1)[which(colnames(soil_predicted_ko_picrust1) %in% colnames(soil_predicted_ko))]
soil_overlapping_samples <- soil_overlapping_samples[which(soil_overlapping_samples %in% colnames(soil_predicted_ko_tax4fun))]
soil_overlapping_samples <- soil_overlapping_samples[which(soil_overlapping_samples %in% colnames(soil_mgs_ko))]

# Subset to specific columns.
soil_predicted_ko <- soil_predicted_ko[,soil_overlapping_samples]
soil_predicted_ko_gg97 <- soil_predicted_ko_gg97[,soil_overlapping_samples]
soil_predicted_ko_picrust1 <- soil_predicted_ko_picrust1[,soil_overlapping_samples]
soil_predicted_ko_panfp <- soil_predicted_ko_panfp[,soil_overlapping_samples]
soil_predicted_ko_piphillin <- soil_predicted_ko_piphillin[,soil_overlapping_samples]
soil_predicted_ko_tax4fun <- soil_predicted_ko_tax4fun[,soil_overlapping_samples]
soil_mgs_ko <- soil_mgs_ko[,soil_overlapping_samples]
soil_predicted_pathabun <- soil_predicted_pathabun[,soil_overlapping_samples]
soil_predicted_pathcov <- soil_predicted_pathcov[,soil_overlapping_samples]
soil_mgs_pathabun <- soil_mgs_pathabun[,soil_overlapping_samples]
soil_mgs_pathcov <- soil_mgs_pathcov[,soil_overlapping_samples]

### Calculate Cosine similarity and Spearman's correlation between 16S and MGS samples.
hmp_ko_mgs_null_cosine <- rand_sample_vs_func_table(db = ko, tab = hmp_mgs_ko, metric="cosine")
hmp_predicted_ko_picrust1_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_ko_picrust1, tab2 = hmp_mgs_ko, cat_string="PICRUSt1", metric="cosine")
hmp_ko_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_ko, tab2 = hmp_mgs_ko, cat_string="PICRUSt2", metric="cosine")
hmp_ko_panfp_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_ko_panfp, tab2 = hmp_mgs_ko, cat_string="PanFP", metric="cosine")
hmp_ko_piphillin_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_ko_piphillin, tab2 = hmp_mgs_ko, cat_string="Piphillin", metric="cosine")
hmp_ko_tax4fun_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_ko_tax4fun, tab2 = hmp_mgs_ko, cat_string="Tax4Fun", metric="cosine")
hmp_ko_cosine_df <- rbind(hmp_ko_mgs_null_cosine, hmp_predicted_ko_picrust1_vs_mgs_cosine, hmp_ko_picrust2_vs_mgs_cosine,
                       hmp_ko_panfp_vs_mgs_cosine, hmp_ko_piphillin_vs_mgs_cosine, hmp_ko_tax4fun_vs_mgs_cosine)

hmp_ko_mgs_null_spearman <- rand_sample_vs_func_table(db = ko, tab = hmp_mgs_ko, metric="spearman")
hmp_predicted_ko_picrust1_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_ko_picrust1, tab2 = hmp_mgs_ko, cat_string="PICRUSt1", metric="spearman")
hmp_ko_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_ko, tab2 = hmp_mgs_ko, cat_string="PICRUSt2", metric="spearman")
hmp_ko_panfp_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_ko_panfp, tab2 = hmp_mgs_ko, cat_string="PanFP", metric="spearman")
hmp_ko_piphillin_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_ko_piphillin, tab2 = hmp_mgs_ko, cat_string="Piphillin", metric="spearman")
hmp_ko_tax4fun_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_ko_tax4fun, tab2 = hmp_mgs_ko, cat_string="Tax4Fun", metric="spearman")
hmp_ko_spearman_df <- rbind(hmp_ko_mgs_null_spearman, hmp_predicted_ko_picrust1_vs_mgs_spearman, hmp_ko_picrust2_vs_mgs_spearman,
                            hmp_ko_panfp_vs_mgs_spearman, hmp_ko_piphillin_vs_mgs_spearman, hmp_ko_tax4fun_vs_mgs_spearman)


hmp_pathabun_mgs_null_cosine <- rand_sample_vs_func_table(db = pathabun, tab = hmp_mgs_pathabun, metric="cosine")
hmp_pathabun_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_pathabun, tab2 = hmp_mgs_pathabun, cat_string="PICRUSt2", metric="cosine")

hmp_pathabun_mgs_null_spearman <- rand_sample_vs_func_table(db = pathabun, tab = hmp_mgs_pathabun, metric="spearman")
hmp_pathabun_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_pathabun, tab2 = hmp_mgs_pathabun, cat_string="PICRUSt2", metric="spearman")

hmp_pathcov_mgs_null_cosine <- rand_sample_vs_func_table(db = pathcov, tab = hmp_mgs_pathcov, metric="cosine")
hmp_pathcov_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_pathcov, tab2 = hmp_mgs_pathcov, cat_string="PICRUSt2", metric="cosine")

hmp_pathcov_mgs_null_spearman <- rand_sample_vs_func_table(db = pathcov, tab = hmp_mgs_pathcov, metric="spearman")
hmp_pathcov_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_pathcov, tab2 = hmp_mgs_pathcov, cat_string="PICRUSt2", metric="spearman")


# mammal samples
mammal_ko_mgs_null_cosine <- rand_sample_vs_func_table(db = ko, tab = mammal_mgs_ko, metric="cosine")
mammal_predicted_ko_picrust1_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_ko_picrust1, tab2 = mammal_mgs_ko, cat_string="PICRUSt1", metric="cosine")
mammal_ko_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_ko, tab2 = mammal_mgs_ko, cat_string="PICRUSt2", metric="cosine")
mammal_ko_panfp_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_ko_panfp, tab2 = mammal_mgs_ko, cat_string="PanFP", metric="cosine")
mammal_ko_piphillin_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_ko_piphillin, tab2 = mammal_mgs_ko, cat_string="Piphillin", metric="cosine")
mammal_ko_tax4fun_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_ko_tax4fun, tab2 = mammal_mgs_ko, cat_string="Tax4Fun", metric="cosine")
mammal_ko_cosine_df <- rbind(mammal_ko_mgs_null_cosine, mammal_predicted_ko_picrust1_vs_mgs_cosine, mammal_ko_picrust2_vs_mgs_cosine,
                          mammal_ko_panfp_vs_mgs_cosine, mammal_ko_piphillin_vs_mgs_cosine, mammal_ko_tax4fun_vs_mgs_cosine)

mammal_ko_mgs_null_spearman <- rand_sample_vs_func_table(db = ko, tab = mammal_mgs_ko, metric="spearman")
mammal_predicted_ko_picrust1_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_ko_picrust1, tab2 = mammal_mgs_ko, cat_string="PICRUSt1", metric="spearman")
mammal_ko_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_ko, tab2 = mammal_mgs_ko, cat_string="PICRUSt2", metric="spearman")
mammal_ko_panfp_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_ko_panfp, tab2 = mammal_mgs_ko, cat_string="PanFP", metric="spearman")
mammal_ko_piphillin_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_ko_piphillin, tab2 = mammal_mgs_ko, cat_string="Piphillin", metric="spearman")
mammal_ko_tax4fun_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_ko_tax4fun, tab2 = mammal_mgs_ko, cat_string="Tax4Fun", metric="spearman")
mammal_ko_spearman_df <- rbind(mammal_ko_mgs_null_spearman, mammal_predicted_ko_picrust1_vs_mgs_spearman, mammal_ko_picrust2_vs_mgs_spearman,
                            mammal_ko_panfp_vs_mgs_spearman, mammal_ko_piphillin_vs_mgs_spearman, mammal_ko_tax4fun_vs_mgs_spearman)


mammal_pathabun_mgs_null_cosine <- rand_sample_vs_func_table(db = pathabun, tab = mammal_mgs_pathabun, metric="cosine")
mammal_pathabun_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_pathabun, tab2 = mammal_mgs_pathabun, cat_string="PICRUSt2", metric="cosine")

mammal_pathabun_mgs_null_spearman <- rand_sample_vs_func_table(db = pathabun, tab = mammal_mgs_pathabun, metric="spearman")
mammal_pathabun_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_pathabun, tab2 = mammal_mgs_pathabun, cat_string="PICRUSt2", metric="spearman")

mammal_pathcov_mgs_null_cosine <- rand_sample_vs_func_table(db = pathcov, tab = mammal_mgs_pathcov, metric="cosine")
mammal_pathcov_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_pathcov, tab2 = mammal_mgs_pathcov, cat_string="PICRUSt2", metric="cosine")

mammal_pathcov_mgs_null_spearman <- rand_sample_vs_func_table(db = pathcov, tab = mammal_mgs_pathcov, metric="spearman")
mammal_pathcov_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_pathcov, tab2 = mammal_mgs_pathcov, cat_string="PICRUSt2", metric="spearman")


# Make comparisons for ocean dataset..
ocean_ko_mgs_null_cosine <- rand_sample_vs_func_table(db = ko, tab = ocean_mgs_ko, metric="cosine")
ocean_predicted_ko_picrust1_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_ko_picrust1, tab2 = ocean_mgs_ko, cat_string="PICRUSt1", metric="cosine")
ocean_ko_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_ko, tab2 = ocean_mgs_ko, cat_string="PICRUSt2", metric="cosine")
ocean_ko_panfp_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_ko_panfp, tab2 = ocean_mgs_ko, cat_string="PanFP", metric="cosine")
ocean_ko_piphillin_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_ko_piphillin, tab2 = ocean_mgs_ko, cat_string="Piphillin", metric="cosine")
ocean_ko_tax4fun_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_ko_tax4fun, tab2 = ocean_mgs_ko, cat_string="Tax4Fun", metric="cosine")
ocean_ko_cosine_df <- rbind(ocean_ko_mgs_null_cosine, ocean_predicted_ko_picrust1_vs_mgs_cosine, ocean_ko_picrust2_vs_mgs_cosine,
                          ocean_ko_panfp_vs_mgs_cosine, ocean_ko_piphillin_vs_mgs_cosine, ocean_ko_tax4fun_vs_mgs_cosine)

ocean_ko_mgs_null_spearman <- rand_sample_vs_func_table(db = ko, tab = ocean_mgs_ko, metric="spearman")
ocean_predicted_ko_picrust1_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_ko_picrust1, tab2 = ocean_mgs_ko, cat_string="PICRUSt1", metric="spearman")
ocean_ko_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_ko, tab2 = ocean_mgs_ko, cat_string="PICRUSt2", metric="spearman")
ocean_ko_panfp_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_ko_panfp, tab2 = ocean_mgs_ko, cat_string="PanFP", metric="spearman")
ocean_ko_piphillin_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_ko_piphillin, tab2 = ocean_mgs_ko, cat_string="Piphillin", metric="spearman")
ocean_ko_tax4fun_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_ko_tax4fun, tab2 = ocean_mgs_ko, cat_string="Tax4Fun", metric="spearman")
ocean_ko_spearman_df <- rbind(ocean_ko_mgs_null_spearman, ocean_predicted_ko_picrust1_vs_mgs_spearman, ocean_ko_picrust2_vs_mgs_spearman,
                            ocean_ko_panfp_vs_mgs_spearman, ocean_ko_piphillin_vs_mgs_spearman, ocean_ko_tax4fun_vs_mgs_spearman)


ocean_pathabun_mgs_null_cosine <- rand_sample_vs_func_table(db = pathabun, tab = ocean_mgs_pathabun, metric="cosine")
ocean_pathabun_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_pathabun, tab2 = ocean_mgs_pathabun, cat_string="PICRUSt2", metric="cosine")

ocean_pathabun_mgs_null_spearman <- rand_sample_vs_func_table(db = pathabun, tab = ocean_mgs_pathabun, metric="spearman")
ocean_pathabun_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_pathabun, tab2 = ocean_mgs_pathabun, cat_string="PICRUSt2", metric="spearman")

ocean_pathcov_mgs_null_cosine <- rand_sample_vs_func_table(db = pathcov, tab = ocean_mgs_pathcov, metric="cosine")
ocean_pathcov_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_pathcov, tab2 = ocean_mgs_pathcov, cat_string="PICRUSt2", metric="cosine")

ocean_pathcov_mgs_null_spearman <- rand_sample_vs_func_table(db = pathcov, tab = ocean_mgs_pathcov, metric="spearman")
ocean_pathcov_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_pathcov, tab2 = ocean_mgs_pathcov, cat_string="PICRUSt2", metric="spearman")


# Make comparisons for soil dataset..
soil_ko_mgs_null_cosine <- rand_sample_vs_func_table(db = ko, tab = soil_mgs_ko, metric="cosine")
soil_predicted_ko_picrust1_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_ko_picrust1, tab2 = soil_mgs_ko, cat_string="PICRUSt1", metric="cosine")
soil_ko_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_ko, tab2 = soil_mgs_ko, cat_string="PICRUSt2", metric="cosine")
soil_ko_panfp_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_ko_panfp, tab2 = soil_mgs_ko, cat_string="PanFP", metric="cosine")
soil_ko_piphillin_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_ko_piphillin, tab2 = soil_mgs_ko, cat_string="Piphillin", metric="cosine")
soil_ko_tax4fun_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_ko_tax4fun, tab2 = soil_mgs_ko, cat_string="Tax4Fun", metric="cosine")
soil_ko_cosine_df <- rbind(soil_ko_mgs_null_cosine, soil_predicted_ko_picrust1_vs_mgs_cosine, soil_ko_picrust2_vs_mgs_cosine,
                          soil_ko_panfp_vs_mgs_cosine, soil_ko_piphillin_vs_mgs_cosine, soil_ko_tax4fun_vs_mgs_cosine)

soil_ko_mgs_null_spearman <- rand_sample_vs_func_table(db = ko, tab = soil_mgs_ko, metric="spearman")
soil_predicted_ko_picrust1_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_ko_picrust1, tab2 = soil_mgs_ko, cat_string="PICRUSt1", metric="spearman")
soil_ko_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_ko, tab2 = soil_mgs_ko, cat_string="PICRUSt2", metric="spearman")
soil_ko_panfp_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_ko_panfp, tab2 = soil_mgs_ko, cat_string="PanFP", metric="spearman")
soil_ko_piphillin_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_ko_piphillin, tab2 = soil_mgs_ko, cat_string="Piphillin", metric="spearman")
soil_ko_tax4fun_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_ko_tax4fun, tab2 = soil_mgs_ko, cat_string="Tax4Fun", metric="spearman")
soil_ko_spearman_df <- rbind(soil_ko_mgs_null_spearman, soil_predicted_ko_picrust1_vs_mgs_spearman, soil_ko_picrust2_vs_mgs_spearman,
                            soil_ko_panfp_vs_mgs_spearman, soil_ko_piphillin_vs_mgs_spearman, soil_ko_tax4fun_vs_mgs_spearman)


soil_pathabun_mgs_null_cosine <- rand_sample_vs_func_table(db = pathabun, tab = soil_mgs_pathabun, metric="cosine")
soil_pathabun_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_pathabun, tab2 = soil_mgs_pathabun, cat_string="PICRUSt2", metric="cosine")

soil_pathabun_mgs_null_spearman <- rand_sample_vs_func_table(db = pathabun, tab = soil_mgs_pathabun, metric="spearman")
soil_pathabun_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_pathabun, tab2 = soil_mgs_pathabun, cat_string="PICRUSt2", metric="spearman")

soil_pathcov_mgs_null_cosine <- rand_sample_vs_func_table(db = pathcov, tab = soil_mgs_pathcov, metric="cosine")
soil_pathcov_picrust2_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_pathcov, tab2 = soil_mgs_pathcov, cat_string="PICRUSt2", metric="cosine")

soil_pathcov_mgs_null_spearman <- rand_sample_vs_func_table(db = pathcov, tab = soil_mgs_pathcov, metric="spearman")
soil_pathcov_picrust2_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_pathcov, tab2 = soil_mgs_pathcov, cat_string="PICRUSt2", metric="spearman")

# Get metrics for KO PICRUSt2 predictions limited to ASVs that match the GG database.
hmp_ko_picrust2_gg97_vs_mgs_cosine <- cor_all_cols(tab1 = hmp_predicted_ko_gg97, tab2 = hmp_mgs_ko, cat_string="Ref-only", metric="cosine")
hmp_ko_picrust2_gg97_vs_mgs_spearman <- cor_all_cols(tab1 = hmp_predicted_ko_gg97, tab2 = hmp_mgs_ko, cat_string="Ref-only", metric="spearman")

mammal_ko_picrust2_gg97_vs_mgs_cosine <- cor_all_cols(tab1 = mammal_predicted_ko_gg97, tab2 = mammal_mgs_ko, cat_string="Ref-only", metric="cosine")
mammal_ko_picrust2_gg97_vs_mgs_spearman <- cor_all_cols(tab1 = mammal_predicted_ko_gg97, tab2 = mammal_mgs_ko, cat_string="Ref-only", metric="spearman")

ocean_ko_picrust2_gg97_vs_mgs_cosine <- cor_all_cols(tab1 = ocean_predicted_ko_gg97, tab2 = ocean_mgs_ko, cat_string="Ref-only", metric="cosine")
ocean_ko_picrust2_gg97_vs_mgs_spearman <- cor_all_cols(tab1 = ocean_predicted_ko_gg97, tab2 = ocean_mgs_ko, cat_string="Ref-only", metric="spearman")

soil_ko_picrust2_gg97_vs_mgs_cosine <- cor_all_cols(tab1 = soil_predicted_ko_gg97, tab2 = soil_mgs_ko, cat_string="Ref-only", metric="cosine")
soil_ko_picrust2_gg97_vs_mgs_spearman <- cor_all_cols(tab1 = soil_predicted_ko_gg97, tab2 = soil_mgs_ko, cat_string="Ref-only", metric="spearman")

# Save RDS objects.
saveRDS(object = hmp_ko_cosine_df, file = "../../saved_RDS/16S_vs_MGS_metrics/hmp_ko_cosine_df.rds")
saveRDS(object = hmp_ko_spearman_df, file = "../../saved_RDS/16S_vs_MGS_metrics/hmp_ko_spearman_df.rds")
saveRDS(object = hmp_ko_picrust2_gg97_vs_mgs_cosine, file = "../../saved_RDS/16S_vs_MGS_metrics/hmp_ko_picrust2_gg97_vs_mgs_cosine.rds")
saveRDS(object = hmp_ko_picrust2_gg97_vs_mgs_spearman, file = "../../saved_RDS/16S_vs_MGS_metrics/hmp_ko_picrust2_gg97_vs_mgs_spearman.rds")

saveRDS(object = mammal_ko_cosine_df, file = "../../saved_RDS/16S_vs_MGS_metrics/mammal_ko_cosine_df.rds")
saveRDS(object = mammal_ko_spearman_df, file = "../../saved_RDS/16S_vs_MGS_metrics/mammal_ko_spearman_df.rds")
saveRDS(object = mammal_ko_picrust2_gg97_vs_mgs_cosine, file = "../../saved_RDS/16S_vs_MGS_metrics/mammal_ko_picrust2_gg97_vs_mgs_cosine.rds")
saveRDS(object = mammal_ko_picrust2_gg97_vs_mgs_spearman, file = "../../saved_RDS/16S_vs_MGS_metrics/mammal_ko_picrust2_gg97_vs_mgs_spearman.rds")

saveRDS(object = ocean_ko_cosine_df, file = "../../saved_RDS/16S_vs_MGS_metrics/ocean_ko_cosine_df.rds")
saveRDS(object = ocean_ko_spearman_df, file = "../../saved_RDS/16S_vs_MGS_metrics/ocean_ko_spearman_df.rds")
saveRDS(object = ocean_ko_picrust2_gg97_vs_mgs_cosine, file = "../../saved_RDS/16S_vs_MGS_metrics/ocean_ko_picrust2_gg97_vs_mgs_cosine.rds")
saveRDS(object = ocean_ko_picrust2_gg97_vs_mgs_spearman, file = "../../saved_RDS/16S_vs_MGS_metrics/ocean_ko_picrust2_gg97_vs_mgs_spearman.rds")

saveRDS(object = soil_ko_cosine_df, file = "../../saved_RDS/16S_vs_MGS_metrics/soil_ko_cosine_df.rds")
saveRDS(object = soil_ko_spearman_df, file = "../../saved_RDS/16S_vs_MGS_metrics/soil_ko_spearman_df.rds")
saveRDS(object = soil_ko_picrust2_gg97_vs_mgs_cosine, file = "../../saved_RDS/16S_vs_MGS_metrics/soil_ko_picrust2_gg97_vs_mgs_cosine.rds")
saveRDS(object = soil_ko_picrust2_gg97_vs_mgs_spearman, file = "../../saved_RDS/16S_vs_MGS_metrics/soil_ko_picrust2_gg97_vs_mgs_spearman.rds")

# Combine path abundance and coverage metrics into dataframes.
hmp_pathabun_cosine_df <- rbind(hmp_pathabun_mgs_null_cosine, hmp_pathabun_picrust2_vs_mgs_cosine)
hmp_pathabun_cosine_df$dataset <- "HMP"

hmp_pathcov_cosine_df <- rbind(hmp_pathcov_mgs_null_cosine, hmp_pathcov_picrust2_vs_mgs_cosine)
hmp_pathcov_cosine_df$dataset <- "HMP"

hmp_pathabun_spearman_df <- rbind(hmp_pathabun_mgs_null_spearman, hmp_pathabun_picrust2_vs_mgs_spearman)
hmp_pathabun_spearman_df$dataset <- "HMP"

hmp_pathcov_spearman_df <- rbind(hmp_pathcov_mgs_null_spearman, hmp_pathcov_picrust2_vs_mgs_spearman)
hmp_pathcov_spearman_df$dataset <- "HMP"

mammal_pathabun_cosine_df <- rbind(mammal_pathabun_mgs_null_cosine, mammal_pathabun_picrust2_vs_mgs_cosine)
mammal_pathabun_cosine_df$dataset <- "Mammal"

mammal_pathcov_cosine_df <- rbind(mammal_pathcov_mgs_null_cosine, mammal_pathcov_picrust2_vs_mgs_cosine)
mammal_pathcov_cosine_df$dataset <- "Mammal"

mammal_pathabun_spearman_df <- rbind(mammal_pathabun_mgs_null_spearman, mammal_pathabun_picrust2_vs_mgs_spearman)
mammal_pathabun_spearman_df$dataset <- "Mammal"

mammal_pathcov_spearman_df <- rbind(mammal_pathcov_mgs_null_spearman, mammal_pathcov_picrust2_vs_mgs_spearman)
mammal_pathcov_spearman_df$dataset <- "Mammal"

ocean_pathabun_cosine_df <- rbind(ocean_pathabun_mgs_null_cosine, ocean_pathabun_picrust2_vs_mgs_cosine)
ocean_pathabun_cosine_df$dataset <- "Ocean"

ocean_pathcov_cosine_df <- rbind(ocean_pathcov_mgs_null_cosine, ocean_pathcov_picrust2_vs_mgs_cosine)
ocean_pathcov_cosine_df$dataset <- "Ocean"

ocean_pathabun_spearman_df <- rbind(ocean_pathabun_mgs_null_spearman, ocean_pathabun_picrust2_vs_mgs_spearman)
ocean_pathabun_spearman_df$dataset <- "Ocean"

ocean_pathcov_spearman_df <- rbind(ocean_pathcov_mgs_null_spearman, ocean_pathcov_picrust2_vs_mgs_spearman)
ocean_pathcov_spearman_df$dataset <- "Ocean"

soil_pathabun_cosine_df <- rbind(soil_pathabun_mgs_null_cosine, soil_pathabun_picrust2_vs_mgs_cosine)
soil_pathabun_cosine_df$dataset <- "Soil"

soil_pathcov_cosine_df <- rbind(soil_pathcov_mgs_null_cosine, soil_pathcov_picrust2_vs_mgs_cosine)
soil_pathcov_cosine_df$dataset <- "Soil"

soil_pathabun_spearman_df <- rbind(soil_pathabun_mgs_null_spearman, soil_pathabun_picrust2_vs_mgs_spearman)
soil_pathabun_spearman_df$dataset <- "Soil"

soil_pathcov_spearman_df <- rbind(soil_pathcov_mgs_null_spearman, soil_pathcov_picrust2_vs_mgs_spearman)
soil_pathcov_spearman_df$dataset <- "Soil"

combined_pathabun_cosine_df <- rbind(hmp_pathabun_cosine_df, mammal_pathabun_cosine_df, ocean_pathabun_cosine_df, soil_pathabun_cosine_df)
combined_pathcov_cosine_df <- rbind(hmp_pathcov_cosine_df, mammal_pathcov_cosine_df, ocean_pathcov_cosine_df, soil_pathcov_cosine_df)
combined_pathabun_spearman_df <- rbind(hmp_pathabun_spearman_df, mammal_pathabun_spearman_df, ocean_pathabun_spearman_df, soil_pathabun_spearman_df)
combined_pathcov_spearman_df <- rbind(hmp_pathcov_spearman_df, mammal_pathcov_spearman_df, ocean_pathcov_spearman_df, soil_pathcov_spearman_df)

saveRDS(object = combined_pathabun_cosine_df, file = "../../saved_RDS/16S_vs_MGS_metrics/combined_pathabun_cosine_df.rds")
saveRDS(object = combined_pathabun_spearman_df, file = "../../saved_RDS/16S_vs_MGS_metrics/combined_pathabun_spearman_df.rds")
saveRDS(object = combined_pathcov_cosine_df, file = "../../saved_RDS/16S_vs_MGS_metrics/combined_pathcov_cosine_df.rds")
saveRDS(object = combined_pathcov_spearman_df, file = "../../saved_RDS/16S_vs_MGS_metrics/combined_pathcov_spearman_df.rds")