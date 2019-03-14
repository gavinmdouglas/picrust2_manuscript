setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

ec_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/mean_func_tables/ec_18S.txt",
                     header=T, sep="\t", row.names=1, check.names = FALSE)

ec_18S_GCF_003184545.1 <- data.frame(t(ec_18S["GCF_003184545.1",]))
ec_18S_GCF_003184545.1 <- ec_18S_GCF_003184545.1[-which(ec_18S_GCF_003184545.1 == 0),, drop=FALSE]

tmp_ec <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/humann2_test/humann2_output/GCF_003184545.1_Asphet1_cds_from_genomic_genefamilies_level4ec_unstratified.tsv",
                     header=F, sep="\t", row.names=1, check.names=FALSE)
tmp_ec <- tmp_ec[-c(1,2),, drop=FALSE]
rownames(tmp_ec) <- gsub("^", "EC:", rownames(tmp_ec))

rownames(ec_18S_GCF_003184545.1)[which(! rownames(ec_18S_GCF_003184545.1) %in% rownames(tmp_ec))]
overlap_ec <- 

ec_blue_16S_nsti2 <- read.table("16S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_18S_nsti2 <- read.table("18S/picrust2_full_output_pipeline/ec_18S_metagenome_out/pred_metagenome_unstrat.tsv",
                                header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_mgs <- read.table("mgs/humann2_final_out/humann2_level4ec_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")

colnames(ec_blue_16S_nsti2) <- gsub("Bact", "Blue", colnames(ec_blue_16S_nsti2))
colnames(ec_blue_18S_nsti2) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti2))

rownames(ec_blue_mgs) <- gsub("^", "EC:", rownames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("_Abundance.RPKs", "", colnames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("BB", "Blue", colnames(ec_blue_mgs))

# Add in ECs called as 0.
possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt")$V1
possible_16S_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/picrust2_ECs.txt")$V1
possible_18S_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/picrust2_ECs_18S.txt")$V1
possible_picrust2_ecs <- factor(unique(c(levels(possible_18S_ecs), levels(possible_16S_ecs))))

ec_blue_mgs <- ec_blue_mgs[-which(rownames(ec_blue_mgs) %in% c("EC:UNMAPPED", "EC:UNGROUPED")),]

# Subset to ECs overlapping across all 3 tables.
ec_blue_mgs_nomiss <- add_missing_funcs(ec_blue_mgs, possible_mgs_ecs)
ec_blue_16S_nsti2_nomiss <- add_missing_funcs(ec_blue_16S_nsti2, possible_picrust2_ecs)
ec_blue_18S_nsti2_nomiss <- add_missing_funcs(ec_blue_18S_nsti2, possible_18S_ecs)

overlapping_ec <- rownames(ec_blue_16S_nsti2_nomiss)[which(rownames(ec_blue_16S_nsti2_nomiss) %in% rownames(ec_blue_mgs_nomiss))]
overlapping_ec <- overlapping_ec[which(overlapping_ec %in% rownames(ec_blue_18S_nsti2_nomiss))]

overlapping_col <- colnames(ec_blue_16S_nsti2_nomiss)[which(colnames(ec_blue_16S_nsti2_nomiss) %in% colnames(ec_blue_mgs_nomiss))]
overlapping_col <- overlapping_col[which(overlapping_col %in% colnames(ec_blue_18S_nsti2_nomiss))]

ec_blue_mgs_nomiss_subset <- ec_blue_mgs_nomiss[overlapping_ec, overlapping_col]
ec_blue_16S_nsti2_nomiss_subset <- ec_blue_16S_nsti2_nomiss[overlapping_ec, overlapping_col]
ec_blue_18S_nsti2_nomiss_subset <- ec_blue_18S_nsti2_nomiss[overlapping_ec, overlapping_col]

# Determine if a combination of 16S and 18S predicitons does better than 16S alone.
ec_blue_16S_relabun <-  data.frame(sweep(ec_blue_16S_nsti2_nomiss_subset, 2, colSums(ec_blue_16S_nsti2_nomiss_subset), '/'))
ec_blue_18S_relabun <-  data.frame(sweep(ec_blue_18S_nsti2_nomiss_subset, 2, colSums(ec_blue_18S_nsti2_nomiss_subset), '/'))
ec_blue_mgs_relabun <-  data.frame(sweep(ec_blue_mgs_nomiss_subset, 2, colSums(ec_blue_mgs_nomiss_subset), '/'))

ec_18S_REF_mean_overlap <- colMeans(ec_18S)[overlapping_ec]

# Function to combine 2 sets of predictions based on different levels
# of subsamplings of each and compare each subsampling to MGS.
combine_predictions <- function(prediction1, prediction2, mgs) {
  
  prediction1 <- ec_blue_16S_relabun$Blue191
  prediction2 <-  ec_blue_18S_relabun$Blue191
  mgs <- ec_blue_mgs_nomiss_subset$Blue191
  ran1 <- rnorm(length(prediction1))
  
  contrib <- c()
  rho <- c()
  precision <- c()
  recall <- c()
  F1 <- c()
  TP <- c()
  TN <- c()
  FP <- c()
  FN <- c()
  
  mgs_binary <- mgs > 0
  
  upper_limit = 10000
  step_size = 50
  
  for(i in seq(from = 0, to = upper_limit, by = 10)) {
    
    for(j in 1:10) {
      
      contrib <- c(contrib, i)
      exp_i <- c(1:length(prediction1))
      
      prediction1_sample <- table(sample(1:length(prediction1), prob=prediction1/sum(prediction1), replace=TRUE, size=i))
      prediction2_sample <- table(sample(1:length(prediction2), prob=prediction2/sum(prediction2), replace=TRUE, size=upper_limit-i))
      
      prediction1_sample_missing_i <- exp_i[which(! exp_i %in% names(prediction1_sample))]
      if(length(prediction1_sample_missing_i) > 0) {
        prediction1_sample_missing_filled <- rep.int(0, length(prediction1_sample_missing_i))
        names(prediction1_sample_missing_filled) <- prediction1_sample_missing_i
        prediction1_sample_all <- c(prediction1_sample, prediction1_sample_missing_filled)
        prediction1_sample_all <- prediction1_sample_all[as.character(exp_i)]
      } else {
        prediction1_sample_all <- prediction1_sample
      }
      
      prediction2_sample_missing_i <- exp_i[which(! exp_i %in% names(prediction2_sample))]
      if(length(prediction2_sample_missing_i) > 0) {
        prediction2_sample_missing_filled <- rep.int(0, length(prediction2_sample_missing_i))
        names(prediction2_sample_missing_filled) <- prediction2_sample_missing_i
        prediction2_sample_all <- c(prediction2_sample, prediction2_sample_missing_filled)
        prediction2_sample_all <- prediction2_sample_all[as.character(exp_i)]
      } else {
        prediction2_sample_all <- prediction2_sample
      }
      
      combined_predict <- prediction1_sample_all + prediction2_sample_all
      
      rho <- c(rho, cor.test(combined_predict, mgs, method="spearman")$estimate)
      
      combined_predict_binary <- combined_predict > 0
      
      num_true_pos <- length(which(which(mgs_binary) %in% which(combined_predict_binary)))
      num_true_neg <- length(which(which(! mgs_binary) %in% which(! combined_predict_binary)))
      
      num_false_pos <- length(which(which(! mgs_binary) %in% which(combined_predict_binary)))
      num_false_neg <- length(which(which(mgs_binary) %in% which(! combined_predict_binary)))
      
      precision_tmp <- num_true_pos/(num_true_pos + num_false_pos)
      precision <- c(precision, precision_tmp)
      
      recall_tmp <- num_true_pos/(num_true_pos + num_false_neg)
      recall <- c(recall, recall_tmp)
      
      F1 <- c(F1, 2 * ((precision_tmp * recall_tmp)/(precision_tmp + recall_tmp)))
      
      TP <- c(TP, num_true_pos)
      TN <- c(TN, num_true_neg)
      FP <- c(FP, num_false_pos)
      FN <- c(FN, num_false_neg)
      
    }
  }
  
  return(data.frame(contrib=contrib,
                    rho=rho,
                    precision=precision,
                    recall=recall,
                    F1=F1,
                    TP=TP,
                    TN=TN,
                    FP=FP,
                    FN=FN))
  
}