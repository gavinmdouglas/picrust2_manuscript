library(ggplot2)
library(reshape2)

setwd("/home/gavin/projects/picrust_pipeline/IMG_phenotypes/combined_table/LOOCV/picrust2_loocv_out/")

# Function to calc accuracy metrics (and allow for NAs)
calc_accuracy_metrics_NA <- function(df1, df2, category) {
  
  # Subset only to columns and rows that overlap between both and convert to present (TRUE) and absent (FALSE)
  cols2keep <- colnames(df1)[which(colnames(df1) %in% colnames(df2))]
  rows2keep <- rownames(df1)[which(rownames(df1) %in% rownames(df2))]
  
  df1 <- df1[rows2keep, cols2keep] > 0
  df2 <- df2[rows2keep, cols2keep] > 0
  
  out_df <- data.frame(matrix(NA, nrow=length(cols2keep), ncol=12))
  colnames(out_df) <- c("category", "sample", "acc", "TP", "TN", "FP", "FN", "NPV", "precision", "recall", "fpr", "F1")
  
  row_i = 1
  for(sample in colnames(df1)) {
    
    x1 <- df1[,sample]
    x2 <- df2[,sample]
    na_element <- which(is.na(x1))
    
    if(length(na_element) > 0) {
      x1 <- x1[-na_element]
      x2 <- x2[-na_element]
    }
  
    overall_acc <- sum(x1 == x2)/length(x1)
    
    num_true_pos <- length(which(which(x1) %in% which(x2)))
    num_true_neg <- length(which(which(! x1) %in% which(! x2)))
    
    num_false_pos <- length(which(which(! x1) %in% which(x2)))
    num_false_neg <- length(which(which(x1) %in% which(! x2)))
    
    npv <- num_true_neg/(num_true_neg + num_false_neg)
    precision <- num_true_pos/(num_true_pos + num_false_pos)
    recall <- num_true_pos/(num_true_pos + num_false_neg)
    fpr <- num_false_pos/(num_false_pos + num_true_neg)
    F1 <-  2 * ((precision * recall)/(precision + recall))
    
    out_df[row_i, ] <- c(NA, NA, overall_acc, num_true_pos, num_true_neg, num_false_pos, num_false_neg,
                         npv, precision, recall, fpr, F1)
    
    row_i = row_i + 1
  }
  
  out_df$category <- category
  out_df$sample <- cols2keep
  
  return(out_df)
}

pheno_loocv_files <- list.files(".")
pheno_loocv_df <- read.table(pheno_loocv_files[1], header=T, sep="\t", stringsAsFactors = FALSE, comment.char="")
pheno_loocv_files <- pheno_loocv_files[-1]

for(f in pheno_loocv_files) {
  file_in <- read.table(f, header=T, sep="\t", stringsAsFactors = FALSE, comment.char="")
  pheno_loocv_df <- rbind(pheno_loocv_df, file_in)
}



phenotype_db_filt <- phenotype_db[,-which(colnames(phenotype_db) %in% c("Growth_on_cellulose_via_cellobiose", "Chlorate_reducer"))]

write.table(x = phenotype_db_filt, file = "/home/gavin/projects/picrust_pipeline/IMG_phenotypes/combined_table/IMG_phenotypes_mean_filt.txt",
            quote = FALSE, sep = "\t",row.names = TRUE, col.names=TRUE)


phenotype_db <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/IMG_phenotypes_mean_filt.txt",
                           header=TRUE, sep="\t", comment.char="", quote="", row.names=1)

rownames(pheno_loocv_df) <- pheno_loocv_df$sequence
pheno_loocv_df <- pheno_loocv_df[rownames(phenotype_db),]
pheno_loocv_df <- pheno_loocv_df[, -1]

# Inititalize null dataframe by scrambling column.
phenotype_null <- phenotype_db
for(phenotype in colnames(phenotype_db)) {
  non_na_index <- which(! is.na(phenotype_db[,phenotype]))
  phenotype_null[non_na_index, phenotype] <- sample(phenotype_db[non_na_index, phenotype])
}
saveRDS(phenotype_null, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/IMG_phenotypes_null.rds")

phenotype_null <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/IMG_phenotypes_null.rds")

phenotype_db_t <- data.frame(t(phenotype_db))
pheno_loocv_df_t <- data.frame(t(pheno_loocv_df))
phenotype_null_t <- data.frame(t(phenotype_null))

acc_by_phenotype <- calc_accuracy_metrics_NA(phenotype_db, pheno_loocv_df, category="PICRUSt2")
acc_by_phenotype_null <- calc_accuracy_metrics_NA(phenotype_db, phenotype_null, category="Null")

combined_acc_by_phenotype <- rbind(acc_by_phenotype, acc_by_phenotype_null)

write.table(x = combined_acc_by_phenotype, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/IMG_pheno_LOOCV_metrics.tsv",
            quote = FALSE, sep = "\t",row.names = TRUE, col.names=TRUE)


########################################
### TESTING
########################################
# Test that values are correct for random phenotype
test_pheno <- pheno_loocv_df[, "L.histidine_prototroph"]
expected_pheno <- phenotype_db[, "L.histidine_prototroph"]

identical(rownames(pheno_loocv_df), rownames(phenotype_db))

# Output:
# category                 sample       acc   TP    TN  FP  FN       NPV precision    recall       fpr        F1
# 20 PICRUSt2 L.histidine_prototroph 0.9592421 1613 16156 364 391 0.9763703 0.8158827 0.8048902 0.0220339 0.8103492

# Accuracy:
length(which(test_pheno$Galactose_utilizing == expected_pheno$Galactose_utilizing))/length(test_pheno$Galactose_utilizing)

test_pheno_binary <- test_pheno > 0
expected_pheno_binary <- expected_pheno > 0

length(which(which(test_pheno_binary) %in% which(expected_pheno_binary)))
# TP: 1613

length(which(which(! test_pheno_binary) %in% which(! expected_pheno_binary)))
# TN:16156

length(which(which(test_pheno_binary) %in% which(! expected_pheno_binary)))
#FP: 364

length(which(which(! test_pheno_binary) %in% which(expected_pheno_binary)))
#FN: 391

#NPV = TN/(TN+FN) = 16156/(391+16156) = 0.9763703
#precision = TP/(TP + FP) = 1613/(1613+364) = 0.8158827
#recall = TP(TP+FN) = 1613/(1613+391) = 0.8048902
#fpr = FP/(FP+TN) = 364/(364 + 16156) = 0.0220339
#f1 = 2*((precision*recall)/(precision+recall)) = 2*((0.8158827*0.8048902)/(0.8158827+0.8048902)) = 0.8103492
########################################


##### TESTING #####
### Note both these function expect vectors of the same length as input.
# Function to calculate Jaccard similarity.
calc_jaccard_sim <- function(x1, x2) {
  overlapping_1 <- length(which(which(x1 == 1) %in% which(x2 == 1)))
  counts_01 <- length(which(which(x1 == 0) %in% which(x2 == 1)))
  counts_10 <- length(which(which(x1 == 1) %in% which(x2 == 0)))
  return(overlapping_1 / (overlapping_1 + counts_01 + counts_10))
}

# Function to calculate simple percent overlap.
calc_percent_overlap <- function(x1, x2) {
  overlapping <- length(which(x1 == x2))
  return((overlapping / length(x1)) * 100)
}

# Get similarities per phenotype.
sim_by_phenotype <- data.frame(matrix(NA, nrow=43, ncol=2))
colnames(sim_by_phenotype) <- c("jaccard", "per_overlap")
rownames(sim_by_phenotype) <- colnames(phenotype_db)


for(phenotype in colnames(phenotype_db)) {
  sim_by_phenotype[phenotype, "jaccard"] <- calc_jaccard_sim(phenotype_db[, phenotype],
                                                             pheno_loocv_df[, phenotype])
  sim_by_phenotype[phenotype, "per_overlap"] <- calc_percent_overlap(phenotype_db[, phenotype],
                                                                     pheno_loocv_df[, phenotype])
}

sim_by_phenotype$phenotype <- rownames(sim_by_phenotype)

boxplot(per_overlap ~ phenotype, data=sim_by_phenotype)

acc_by_phenotype <- calc_accuracy_metrics_NA(phenotype_db, pheno_loocv_df, category="IMG Phenotypes")

acc_by_phenotype_null <- calc_accuracy_metrics_NA(phenotype_db, phenotype_null, category="Null")





rownames(acc_by_sample) <- gsub("^X", "", acc_by_sample$sample)

plot(pheno_loocv_df[rownames(acc_by_sample), "metadata_NSTI"],
     acc_by_sample$precision)


boxplot(acc_by_sample_null$precision, acc_by_sample$precision)
boxplot(acc_by_sample_null$recall, acc_by_sample$recall)
boxplot(acc_by_sample_null$acc, acc_by_sample$acc)



acc_by_sample <- calc_accuracy_metrics_NA(phenotype_db_t, pheno_loocv_df_t, category="PICRUSt2")
acc_by_sample_null <- calc_accuracy_metrics_NA(phenotype_db_t, phenotype_null_t, category="Null")

combined_acc_by_sample <- rbind(acc_by_sample, acc_by_sample_null)

combined_acc_by_sample_subset <- combined_acc_by_sample[,c("sample", "acc", "precision", "recall", "category")]
colnames(combined_acc_by_sample_subset) <- c("Genome", "Accuracy", "Precision", "Recall", "Category")
combined_acc_by_sample_subset <- combined_acc_by_sample_subset[with(combined_acc_by_sample_subset, order(Category, Genome)),]

combined_acc_by_sample_subset$Category <- factor(combined_acc_by_sample_subset$Category, levels=c("Null", "PICRUSt2"))
combined_acc_by_sample_subset_melt <- melt(combined_acc_by_sample_subset)

ggplot(combined_acc_by_sample_subset_melt, aes(x=Category, y=value, fill=Category)) + geom_boxplot() +
  facet_grid(. ~ variable, scales = "free", space = "free", switch="x") + coord_flip() + ylab("") +
  ylim(c(0, 1.0)) + ylab(c(""))  + guides(fill=FALSE) + scale_fill_manual(values=c("light grey", "#00BFC4"))