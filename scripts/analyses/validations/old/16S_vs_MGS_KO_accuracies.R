library(ggplot2)
library(reshape2)

setwd("/home/gavin/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

### Read in possible KOs outputted by each tool.
humann2_kos <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_KOs/humann2_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
panfp_kos <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_KOs/PanFP_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
picrust1_kos <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_KOs/PICRUSt1_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
picrust2_kos <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_KOs/PICRUSt2_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
piphillin_kos <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_KOs/piphillin_KOs.txt", header=F, stringsAsFactors = FALSE)$V1
tax4fun_kos <- read.table("/home/gavin/projects/picrust2_manuscript/data/16S_validation/possible_KOs/Tax4Fun_KOs.txt", header=F, stringsAsFactors = FALSE)$V1

### Read in hmp tables.
hmp_predicted_ko <- read.table("picrust2_out/hmp_picrust2_ko.tsv",
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
hmp_predicted_ko_picrust1 <- hmp_predicted_ko_picrust1[,hmp_overlapping_samples]
hmp_predicted_ko_panfp <- hmp_predicted_ko_panfp[,hmp_overlapping_samples]
hmp_predicted_ko_piphillin <- hmp_predicted_ko_piphillin[,hmp_overlapping_samples]
hmp_predicted_ko_tax4fun <- hmp_predicted_ko_tax4fun[,hmp_overlapping_samples]
hmp_mgs_ko <- hmp_mgs_ko[,hmp_overlapping_samples]
hmp_predicted_pathabun <- hmp_predicted_pathabun[,hmp_overlapping_samples]
hmp_predicted_pathcov <- hmp_predicted_pathcov[,hmp_overlapping_samples]
hmp_mgs_pathabun <- hmp_mgs_pathabun[,hmp_overlapping_samples]
hmp_mgs_pathcov <- hmp_mgs_pathcov[,hmp_overlapping_samples]



# Get numbers of which KOs overlap with MGS KO set (i.e. the potential KOs that could be called as present).
length(which(rownames(hmp_predicted_ko_all) %in% rownames(hmp_mgs_ko_all)))
length(which(rownames(hmp_predicted_ko_picrust1_all) %in% rownames(hmp_mgs_ko_all)))
length(which(rownames(hmp_predicted_ko_panfp_all) %in% rownames(hmp_mgs_ko_all)))
length(which(rownames(hmp_predicted_ko_piphillin_all) %in% rownames(hmp_mgs_ko_all)))
length(which(rownames(hmp_predicted_ko_tax4fun_all) %in% rownames(hmp_mgs_ko_all)))

# Accuracy metrics based on presence/absence calls of functions.
# Df1 is assumed to be the gold standard.
calc_accuracy_metrics <- function(df1, df2) {

  # Subset only to columns and rows that overlap between both and convert to present (TRUE) and absent (FALSE)
  cols2keep <- colnames(df1)[which(colnames(df1) %in% colnames(df2))]
  rows2keep <- rownames(df1)[which(rownames(df1) %in% rownames(df2))]
  
  df1 <- df1[rows2keep, cols2keep] > 0
  df2 <- df2[rows2keep, cols2keep] > 0
  
  out_df <- data.frame(matrix(NA, nrow=length(cols2keep), ncol=11))
  colnames(out_df) <- c("sample", "acc", "TP", "TN", "FP", "FN", "PPV", "NPV", "precision", "recall", "fpr")
  
  row_i = 1
  for(sample in colnames(df1)) {
   
    total_func <- length(rows2keep)
    
    overall_acc <- sum(df1[,sample] == df2[,sample])/total_func
  
    num_true_pos <- length(which(which(df1[,sample]) %in% which(df2[,sample])))
    num_true_neg <- length(which(which(! df1[,sample]) %in% which(! df2[,sample])))
    
    num_false_pos <- length(which(which(! df1[,sample]) %in% which(df2[,sample])))
    num_false_neg <- length(which(which(df1[,sample]) %in% which(! df2[,sample])))
    
    ppv <- num_true_pos/(num_true_pos + num_false_pos)
    npv <- num_true_neg/(num_true_neg + num_false_neg)
    precision <- num_true_pos/(num_true_pos + num_false_pos)
    recall <- num_true_pos/(num_true_pos + num_false_neg)
    fpr <- num_false_pos/(num_false_pos + num_true_neg)
    
    out_df[row_i, ] <- c(NA, overall_acc, num_true_pos, num_true_neg, num_false_pos, num_false_neg,
                          ppv, npv, precision, recall, fpr)
    
    row_i = row_i + 1
  }
  
  out_df$sample <- cols2keep
  
  return(out_df)
}

hmp_mgs_ko_all <- add_missing_funcs(hmp_mgs_ko, humann2_kos)
hmp_predicted_ko_all <- add_missing_funcs(hmp_predicted_ko, picrust2_kos)
hmp_predicted_ko_picrust1_all <- add_missing_funcs(hmp_predicted_ko_picrust1, picrust1_kos)
hmp_predicted_ko_panfp_all <- add_missing_funcs(hmp_predicted_ko_panfp, panfp_kos)
hmp_predicted_ko_piphillin_all <- add_missing_funcs(hmp_predicted_ko_piphillin, piphillin_kos)
hmp_predicted_ko_tax4fun_all <- add_missing_funcs(hmp_predicted_ko_tax4fun, tax4fun_kos)

hmp_picrust2_metrics <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_all)
hmp_picrust1_metrics <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_picrust1_all)
hmp_panfp_metrics <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_panfp_all)
hmp_piphillin_metrics <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_piphillin_all)
hmp_tax4fun_metrics <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_tax4fun_all)

# Subset to kos overlapping in all tools.
hmp_overlapping_kos <- rownames(hmp_predicted_ko_all)[which(rownames(hmp_predicted_ko_all) %in% rownames(hmp_predicted_ko_picrust1_all))]
hmp_overlapping_kos <- hmp_overlapping_kos[which(hmp_overlapping_kos %in% rownames(hmp_predicted_ko_panfp_all))]
hmp_overlapping_kos <- hmp_overlapping_kos[which(hmp_overlapping_kos %in% rownames(hmp_predicted_ko_piphillin_all))]
hmp_overlapping_kos <- hmp_overlapping_kos[which(hmp_overlapping_kos %in% rownames(hmp_predicted_ko_tax4fun_all))]

hmp_predicted_ko_overlapKOs <- hmp_predicted_ko_all[hmp_overlapping_kos,]
hmp_predicted_ko_picrust1_overlapKOs <- hmp_predicted_ko_picrust1_all[hmp_overlapping_kos,]
hmp_predicted_ko_panfp_overlapKOs <- hmp_predicted_ko_panfp_all[hmp_overlapping_kos,]
hmp_predicted_ko_piphillin_overlapKOs <- hmp_predicted_ko_piphillin_all[hmp_overlapping_kos,]
hmp_predicted_ko_tax4fun_overlapKOs <- hmp_predicted_ko_tax4fun_all[hmp_overlapping_kos,]

hmp_picrust2_metrics_overlapKOs <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_overlapKOs)
hmp_picrust1_metrics_overlapKOs <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_picrust1_overlapKOs)
hmp_panfp_metrics_overlapKOs <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_panfp_overlapKOs)
hmp_piphillin_metrics_overlapKOs <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_piphillin_overlapKOs)
hmp_tax4fun_metrics_overlapKOs <- calc_accuracy_metrics(hmp_mgs_ko_all, hmp_predicted_ko_tax4fun_overlapKOs)

hmp_picrust2_metrics_overlapKOs$tool <- "PICRUSt2"
hmp_picrust1_metrics_overlapKOs$tool <- "PICRUSt1"
hmp_panfp_metrics_overlapKOs$tool <- "PanFP"
hmp_piphillin_metrics_overlapKOs$tool <- "Piphillin"
hmp_tax4fun_metrics_overlapKOs$tool <- "Tax4Fun"


hmp_combined_metrics_overlapKOs <- rbind(hmp_picrust2_metrics_overlapKOs, hmp_picrust1_metrics_overlapKOs,
                                         hmp_panfp_metrics_overlapKOs, hmp_piphillin_metrics_overlapKOs,
                                         hmp_tax4fun_metrics_overlapKOs)

hmp_combined_metrics_overlapKOs <- hmp_combined_metrics_overlapKOs[,c("sample", "acc", "PPV", "NPV", "precision", "recall", "fpr", "tool")]

hmp_combined_metrics_overlapKOs$tool <- factor(hmp_combined_metrics_overlapKOs$tool, levels=c("Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))

hmp_combined_metrics_overlapKOs_melt <- melt(hmp_combined_metrics_overlapKOs)

ggplot(hmp_combined_metrics_overlapKOs_melt, aes(x=tool, y=value)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x")
  
  
boxplot(picrust2_metrics_overlapKOs$acc, picrust1_metrics_overlapKOs$acc,
        piphillin_metrics_overlapKOs$acc, panfp_metrics_overlapKOs$acc, tax4fun_metrics_overlapKOs$acc)


mammal_mgs_ko_all <- add_missing_funcs(mammal_mgs_ko, humann2_kos)
mammal_predicted_ko_all <- add_missing_funcs(mammal_predicted_ko, picrust2_kos)
mammal_predicted_ko_nsti1.5_all <- add_missing_funcs(mammal_predicted_ko_nsti1.5, picrust2_kos)
mammal_predicted_ko_nsti1_all <- add_missing_funcs(mammal_predicted_ko_nsti1, picrust2_kos)
mammal_predicted_ko_nsti0.5_all <- add_missing_funcs(mammal_predicted_ko_nsti0.5, picrust2_kos)
mammal_predicted_ko_nsti0.25_all <- add_missing_funcs(mammal_predicted_ko_nsti0.25, picrust2_kos)
mammal_predicted_ko_nsti0.125_all <- add_missing_funcs(mammal_predicted_ko_nsti0.125, picrust2_kos)
mammal_predicted_ko_nsti0.0625_all <- add_missing_funcs(mammal_predicted_ko_nsti0.0625, picrust2_kos)
mammal_predicted_ko_picrust1_all <- add_missing_funcs(mammal_predicted_ko_picrust1, picrust1_kos)
mammal_predicted_ko_panfp_all <- add_missing_funcs(mammal_predicted_ko_panfp, panfp_kos)
mammal_predicted_ko_piphillin_all <- add_missing_funcs(mammal_predicted_ko_piphillin, piphillin_kos)
mammal_predicted_ko_tax4fun_all <- add_missing_funcs(mammal_predicted_ko_tax4fun, tax4fun_kos)

mammal_picrust2_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_all)
mammal_picrust2_nsti1.5_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_nsti1.5_all)
mammal_picrust2_nsti1_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_nsti1_all)
mammal_picrust2_nsti0.5_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_nsti0.5_all)
mammal_picrust1_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_picrust1_all)
mammal_panfp_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_panfp_all)
mammal_piphillin_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_piphillin_all)
mammal_tax4fun_metrics <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_tax4fun_all)

# Subset to kos overlapping in all tools.
mammal_overlapping_kos <- rownames(mammal_predicted_ko_all)[which(rownames(mammal_predicted_ko_all) %in% rownames(mammal_predicted_ko_picrust1_all))]
mammal_overlapping_kos <- mammal_overlapping_kos[which(mammal_overlapping_kos %in% rownames(mammal_predicted_ko_panfp_all))]
mammal_overlapping_kos <- mammal_overlapping_kos[which(mammal_overlapping_kos %in% rownames(mammal_predicted_ko_piphillin_all))]
mammal_overlapping_kos <- mammal_overlapping_kos[which(mammal_overlapping_kos %in% rownames(mammal_predicted_ko_tax4fun_all))]

mammal_predicted_ko_overlapKOs <- mammal_predicted_ko_all[mammal_overlapping_kos,]
mammal_picrust2_nsti1.5_metrics_overlapKOs <- mammal_predicted_ko_nsti1.5_all[mammal_overlapping_kos,]
mammal_picrust2_nsti1_metrics_overlapKOs <- mammal_predicted_ko_nsti1_all[mammal_overlapping_kos,]
mammal_picrust2_nsti0.5_metrics_overlapKOs <- mammal_predicted_ko_nsti0.5_all[mammal_overlapping_kos,]
mammal_picrust2_nsti0.25_metrics_overlapKOs <- mammal_predicted_ko_nsti0.25_all[mammal_overlapping_kos,]
mammal_picrust2_nsti0.125_metrics_overlapKOs <- mammal_predicted_ko_nsti0.125_all[mammal_overlapping_kos,]
mammal_picrust2_nsti0.0625_metrics_overlapKOs <- mammal_predicted_ko_nsti0.0625_all[mammal_overlapping_kos,]
mammal_predicted_ko_picrust1_overlapKOs <- mammal_predicted_ko_picrust1_all[mammal_overlapping_kos,]
mammal_predicted_ko_panfp_overlapKOs <- mammal_predicted_ko_panfp_all[mammal_overlapping_kos,]
mammal_predicted_ko_piphillin_overlapKOs <- mammal_predicted_ko_piphillin_all[mammal_overlapping_kos,]
mammal_predicted_ko_tax4fun_overlapKOs <- mammal_predicted_ko_tax4fun_all[mammal_overlapping_kos,]

mammal_picrust2_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_overlapKOs)
mammal_picrust2_nsti1.5_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_picrust2_nsti1.5_metrics_overlapKOs)
mammal_picrust2_nsti1_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_picrust2_nsti1_metrics_overlapKOs)
mammal_picrust2_nsti0.5_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_picrust2_nsti0.5_metrics_overlapKOs)
mammal_picrust2_nsti0.25_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_picrust2_nsti0.25_metrics_overlapKOs)
mammal_picrust2_nsti0.125_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_picrust2_nsti0.125_metrics_overlapKOs)
mammal_picrust2_nsti0.0625_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_picrust2_nsti0.0625_metrics_overlapKOs)
mammal_picrust1_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_picrust1_overlapKOs)
mammal_panfp_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_panfp_overlapKOs)
mammal_piphillin_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_piphillin_overlapKOs)
mammal_tax4fun_metrics_overlapKOs <- calc_accuracy_metrics(mammal_mgs_ko_all, mammal_predicted_ko_tax4fun_overlapKOs)

mammal_picrust2_metrics_overlapKOs$tool <- "PICRUSt2"
mammal_picrust2_nsti1.5_metrics_overlapKOs$tool <- "nsti1.5"
mammal_picrust2_nsti1_metrics_overlapKOs$tool <- "nsti1"
mammal_picrust2_nsti0.5_metrics_overlapKOs$tool <- "nsti0.5"
mammal_picrust2_nsti0.25_metrics_overlapKOs$tool <- "nsti0.25"
mammal_picrust2_nsti0.125_metrics_overlapKOs$tool <- "nsti0.125"
mammal_picrust2_nsti0.0625_metrics_overlapKOs$tool <- "nsti0.0625"
mammal_picrust1_metrics_overlapKOs$tool <- "PICRUSt1"
mammal_panfp_metrics_overlapKOs$tool <- "PanFP"
mammal_piphillin_metrics_overlapKOs$tool <- "Piphillin"
mammal_tax4fun_metrics_overlapKOs$tool <- "Tax4Fun"


mammal_combined_metrics_overlapKOs <- rbind(mammal_picrust2_metrics_overlapKOs, mammal_picrust2_nsti1.5_metrics_overlapKOs,
                                            mammal_picrust2_nsti1_metrics_overlapKOs, mammal_picrust2_nsti0.5_metrics_overlapKOs,
                                            mammal_picrust2_nsti0.25_metrics_overlapKOs, mammal_picrust2_nsti0.125_metrics_overlapKOs,
                                            mammal_picrust2_nsti0.0625_metrics_overlapKOs,
                                            mammal_picrust1_metrics_overlapKOs, mammal_panfp_metrics_overlapKOs,
                                            mammal_piphillin_metrics_overlapKOs, mammal_tax4fun_metrics_overlapKOs)

mammal_combined_metrics_overlapKOs <- mammal_combined_metrics_overlapKOs[,c("sample", "acc", "PPV", "NPV", "precision", "recall", "fpr", "tool")]

mammal_combined_metrics_overlapKOs$tool <- factor(mammal_combined_metrics_overlapKOs$tool, levels=c("Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2", "nsti1.5", "nsti1", "nsti0.5", "nsti0.25", "nsti0.125", "nsti0.0625"))

mammal_combined_metrics_overlapKOs_melt <- melt(mammal_combined_metrics_overlapKOs)

ggplot(mammal_combined_metrics_overlapKOs_melt, aes(x=tool, y=value)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x")


ocean_mgs_ko_all <- add_missing_funcs(ocean_mgs_ko, humann2_kos)
ocean_predicted_ko_all <- add_missing_funcs(ocean_predicted_ko, picrust2_kos)
ocean_predicted_ko_picrust1_all <- add_missing_funcs(ocean_predicted_ko_picrust1, picrust1_kos)
ocean_predicted_ko_panfp_all <- add_missing_funcs(ocean_predicted_ko_panfp, panfp_kos)
ocean_predicted_ko_piphillin_all <- add_missing_funcs(ocean_predicted_ko_piphillin, piphillin_kos)
ocean_predicted_ko_tax4fun_all <- add_missing_funcs(ocean_predicted_ko_tax4fun, tax4fun_kos)

ocean_picrust2_metrics <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_all)
ocean_picrust1_metrics <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_picrust1_all)
ocean_panfp_metrics <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_panfp_all)
ocean_piphillin_metrics <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_piphillin_all)
ocean_tax4fun_metrics <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_tax4fun_all)

# Subset to kos overlapping in all tools.
ocean_overlapping_kos <- rownames(ocean_predicted_ko_all)[which(rownames(ocean_predicted_ko_all) %in% rownames(ocean_predicted_ko_picrust1_all))]
ocean_overlapping_kos <- ocean_overlapping_kos[which(ocean_overlapping_kos %in% rownames(ocean_predicted_ko_panfp_all))]
ocean_overlapping_kos <- ocean_overlapping_kos[which(ocean_overlapping_kos %in% rownames(ocean_predicted_ko_piphillin_all))]
ocean_overlapping_kos <- ocean_overlapping_kos[which(ocean_overlapping_kos %in% rownames(ocean_predicted_ko_tax4fun_all))]

ocean_predicted_ko_overlapKOs <- ocean_predicted_ko_all[ocean_overlapping_kos,]
ocean_predicted_ko_picrust1_overlapKOs <- ocean_predicted_ko_picrust1_all[ocean_overlapping_kos,]
ocean_predicted_ko_panfp_overlapKOs <- ocean_predicted_ko_panfp_all[ocean_overlapping_kos,]
ocean_predicted_ko_piphillin_overlapKOs <- ocean_predicted_ko_piphillin_all[ocean_overlapping_kos,]
ocean_predicted_ko_tax4fun_overlapKOs <- ocean_predicted_ko_tax4fun_all[ocean_overlapping_kos,]

ocean_picrust2_metrics_overlapKOs <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_overlapKOs)
ocean_picrust1_metrics_overlapKOs <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_picrust1_overlapKOs)
ocean_panfp_metrics_overlapKOs <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_panfp_overlapKOs)
ocean_piphillin_metrics_overlapKOs <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_piphillin_overlapKOs)
ocean_tax4fun_metrics_overlapKOs <- calc_accuracy_metrics(ocean_mgs_ko_all, ocean_predicted_ko_tax4fun_overlapKOs)

ocean_picrust2_metrics_overlapKOs$tool <- "PICRUSt2"
ocean_picrust1_metrics_overlapKOs$tool <- "PICRUSt1"
ocean_panfp_metrics_overlapKOs$tool <- "PanFP"
ocean_piphillin_metrics_overlapKOs$tool <- "Piphillin"
ocean_tax4fun_metrics_overlapKOs$tool <- "Tax4Fun"


ocean_combined_metrics_overlapKOs <- rbind(ocean_picrust2_metrics_overlapKOs, ocean_picrust1_metrics_overlapKOs,
                                            ocean_panfp_metrics_overlapKOs, ocean_piphillin_metrics_overlapKOs,
                                            ocean_tax4fun_metrics_overlapKOs)

ocean_combined_metrics_overlapKOs <- ocean_combined_metrics_overlapKOs[,c("sample", "acc", "PPV", "NPV", "precision", "recall", "fpr", "tool")]

ocean_combined_metrics_overlapKOs$tool <- factor(ocean_combined_metrics_overlapKOs$tool, levels=c("Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))

ocean_combined_metrics_overlapKOs_melt <- melt(ocean_combined_metrics_overlapKOs)

ggplot(ocean_combined_metrics_overlapKOs_melt, aes(x=tool, y=value)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x")


soil_mgs_ko_all <- add_missing_funcs(soil_mgs_ko, humann2_kos)
soil_predicted_ko_all <- add_missing_funcs(soil_predicted_ko, picrust2_kos)
soil_predicted_ko_picrust1_all <- add_missing_funcs(soil_predicted_ko_picrust1, picrust1_kos)
soil_predicted_ko_panfp_all <- add_missing_funcs(soil_predicted_ko_panfp, panfp_kos)
soil_predicted_ko_piphillin_all <- add_missing_funcs(soil_predicted_ko_piphillin, piphillin_kos)
soil_predicted_ko_tax4fun_all <- add_missing_funcs(soil_predicted_ko_tax4fun, tax4fun_kos)

soil_picrust2_metrics <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_all)
soil_picrust1_metrics <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_picrust1_all)
soil_panfp_metrics <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_panfp_all)
soil_piphillin_metrics <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_piphillin_all)
soil_tax4fun_metrics <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_tax4fun_all)

# Subset to kos overlapping in all tools.
soil_overlapping_kos <- rownames(soil_predicted_ko_all)[which(rownames(soil_predicted_ko_all) %in% rownames(soil_predicted_ko_picrust1_all))]
soil_overlapping_kos <- soil_overlapping_kos[which(soil_overlapping_kos %in% rownames(soil_predicted_ko_panfp_all))]
soil_overlapping_kos <- soil_overlapping_kos[which(soil_overlapping_kos %in% rownames(soil_predicted_ko_piphillin_all))]
soil_overlapping_kos <- soil_overlapping_kos[which(soil_overlapping_kos %in% rownames(soil_predicted_ko_tax4fun_all))]

soil_predicted_ko_overlapKOs <- soil_predicted_ko_all[soil_overlapping_kos,]
soil_predicted_ko_picrust1_overlapKOs <- soil_predicted_ko_picrust1_all[soil_overlapping_kos,]
soil_predicted_ko_panfp_overlapKOs <- soil_predicted_ko_panfp_all[soil_overlapping_kos,]
soil_predicted_ko_piphillin_overlapKOs <- soil_predicted_ko_piphillin_all[soil_overlapping_kos,]
soil_predicted_ko_tax4fun_overlapKOs <- soil_predicted_ko_tax4fun_all[soil_overlapping_kos,]

soil_picrust2_metrics_overlapKOs <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_overlapKOs)
soil_picrust1_metrics_overlapKOs <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_picrust1_overlapKOs)
soil_panfp_metrics_overlapKOs <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_panfp_overlapKOs)
soil_piphillin_metrics_overlapKOs <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_piphillin_overlapKOs)
soil_tax4fun_metrics_overlapKOs <- calc_accuracy_metrics(soil_mgs_ko_all, soil_predicted_ko_tax4fun_overlapKOs)

soil_picrust2_metrics_overlapKOs$tool <- "PICRUSt2"
soil_picrust1_metrics_overlapKOs$tool <- "PICRUSt1"
soil_panfp_metrics_overlapKOs$tool <- "PanFP"
soil_piphillin_metrics_overlapKOs$tool <- "Piphillin"
soil_tax4fun_metrics_overlapKOs$tool <- "Tax4Fun"


soil_combined_metrics_overlapKOs <- rbind(soil_picrust2_metrics_overlapKOs, soil_picrust1_metrics_overlapKOs,
                                           soil_panfp_metrics_overlapKOs, soil_piphillin_metrics_overlapKOs,
                                           soil_tax4fun_metrics_overlapKOs)

soil_combined_metrics_overlapKOs <- soil_combined_metrics_overlapKOs[,c("sample", "acc", "PPV", "NPV", "precision", "recall", "fpr", "tool")]

soil_combined_metrics_overlapKOs$tool <- factor(soil_combined_metrics_overlapKOs$tool, levels=c("Tax4Fun", "PanFP", "Piphillin", "PICRUSt1", "PICRUSt2"))

soil_combined_metrics_overlapKOs_melt <- melt(soil_combined_metrics_overlapKOs)

ggplot(soil_combined_metrics_overlapKOs_melt, aes(x=tool, y=value)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x")







### Read in mammal tables.
mammal_predicted_ko <- read.table("picrust2_out/mammal_picrust2_ko.tsv",
                               header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_nsti1.5 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv",
                                  header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_nsti1 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_nsti_1/pred_metagenome_unstrat.tsv",
                                          header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_nsti0.5 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_nsti0.25 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv",
                                          header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_nsti0.125 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_nsti_0.125/pred_metagenome_unstrat.tsv",
                                          header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, quote="", comment.char="")

mammal_predicted_ko_nsti0.0625 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/picrust2_pipeline/picrust2_full_output/KO_metagenome_out_nsti_0.0625/pred_metagenome_unstrat.tsv",
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
mammal_predicted_ko_nsti1.5 <- mammal_predicted_ko_nsti1.5[,mammal_overlapping_samples]
mammal_predicted_ko_nsti1 <- mammal_predicted_ko_nsti1[,mammal_overlapping_samples]
mammal_predicted_ko_nsti0.5 <- mammal_predicted_ko_nsti0.5[,mammal_overlapping_samples]
mammal_predicted_ko_nsti0.25 <- mammal_predicted_ko_nsti0.25[,mammal_overlapping_samples]
mammal_predicted_ko_nsti0.125 <- mammal_predicted_ko_nsti0.125[,mammal_overlapping_samples]
mammal_predicted_ko_nsti0.0625 <- mammal_predicted_ko_nsti0.0625[,mammal_overlapping_samples]
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