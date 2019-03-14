setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

ec_blue_16S_nsti2 <- read.table("16S/picrust2_full_output_pipeline/EC_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_18S_nsti2 <- read.table("18S/picrust2_full_output_1000fungi/ec_18S_counts_metagenome_out/pred_metagenome_unstrat.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_blue_mgs <- read.table("mgs/humann2_final_out/humann2_level4ec_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")


colnames(ec_blue_16S_nsti2) <- gsub("Bact", "Blue", colnames(ec_blue_16S_nsti2))
colnames(ec_blue_18S_nsti2) <- gsub("^", "Blue", colnames(ec_blue_18S_nsti2))

rownames(ec_blue_mgs) <- gsub("^", "EC:", rownames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("_Abundance.RPKs", "", colnames(ec_blue_mgs))
colnames(ec_blue_mgs) <- gsub("BB", "Blue", colnames(ec_blue_mgs))

# Add in ECs called as 0.
possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt", stringsAsFactors = FALSE)$V1
possible_16S_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/picrust2_ECs.txt", stringsAsFactors = FALSE)$V1
possible_18S_ecs <- colnames(read.table("/home/gavin/projects/picrust_pipeline/fungal_genomes/mean_func_tables/ec_18S_counts.txt", header=T, sep="\t", row.names=1, check.names = FALSE))
possible_picrust2_ecs <- unique(c(possible_18S_ecs, possible_16S_ecs))

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

# Function to combine 2 sets of predictions based on different levels
# of subsamplings of each and compare each subsampling to MGS.
combine_predictions <- function(prediction1, prediction2, mgs) {
 
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

combine_predictions_simple_percent <- function(prediction1, prediction2, mgs) {
 
  spearman <- c()
  
  for(i in seq(from=0, to=1, by=0.01)) {
   
    rep_combo <- prediction1*i + prediction2*(1-i)
    
    spearman <- c(spearman, cor.test(rep_combo, mgs, method="spearman")$estimate)
    
  }
  
  return(data.frame(contrib=seq(from=0, to=1, by=0.01), spearman=spearman))
}

Blue198_16S_18S_scc_contrib_simple <- combine_predictions_simple_percent(ec_blue_16S_relabun$Blue198, ec_blue_18S_relabun$Blue198, ec_blue_mgs_nomiss_subset$Blue198)


Blue190_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue190, ec_blue_18S_relabun$Blue190, ec_blue_mgs_nomiss_subset$Blue190)
Blue191_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue191, ec_blue_18S_relabun$Blue191, ec_blue_mgs_nomiss_subset$Blue191)
Blue192_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue192, ec_blue_18S_relabun$Blue192, ec_blue_mgs_nomiss_subset$Blue192)
Blue193_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue193, ec_blue_18S_relabun$Blue193, ec_blue_mgs_nomiss_subset$Blue193)
Blue195_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue195, ec_blue_18S_relabun$Blue195, ec_blue_mgs_nomiss_subset$Blue195)
Blue197_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue197, ec_blue_18S_relabun$Blue197, ec_blue_mgs_nomiss_subset$Blue197)
Blue198_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue198, ec_blue_18S_relabun$Blue198, ec_blue_mgs_nomiss_subset$Blue198)
Blue200_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue200, ec_blue_18S_relabun$Blue200, ec_blue_mgs_nomiss_subset$Blue200)
Blue202_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue202, ec_blue_18S_relabun$Blue202, ec_blue_mgs_nomiss_subset$Blue202)
Blue203_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue203, ec_blue_18S_relabun$Blue203, ec_blue_mgs_nomiss_subset$Blue203)
Blue204_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue204, ec_blue_18S_relabun$Blue204, ec_blue_mgs_nomiss_subset$Blue204)
Blue207_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue207, ec_blue_18S_relabun$Blue207, ec_blue_mgs_nomiss_subset$Blue207)
Blue208_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue208, ec_blue_18S_relabun$Blue208, ec_blue_mgs_nomiss_subset$Blue208)
Blue209_16S_18S_scc_contrib <- combine_predictions(ec_blue_16S_relabun$Blue209, ec_blue_18S_relabun$Blue209, ec_blue_mgs_nomiss_subset$Blue209)

Blue190_16S_18S_scc_contrib$sample <- "Blue190"
Blue191_16S_18S_scc_contrib$sample <- "Blue191"
Blue192_16S_18S_scc_contrib$sample <- "Blue192"
Blue193_16S_18S_scc_contrib$sample <- "Blue193"
Blue195_16S_18S_scc_contrib$sample <- "Blue195"
Blue197_16S_18S_scc_contrib$sample <- "Blue197"
Blue198_16S_18S_scc_contrib$sample <- "Blue198"
Blue200_16S_18S_scc_contrib$sample <- "Blue200"
Blue202_16S_18S_scc_contrib$sample <- "Blue202"
Blue203_16S_18S_scc_contrib$sample <- "Blue203"
Blue204_16S_18S_scc_contrib$sample <- "Blue204"
Blue207_16S_18S_scc_contrib$sample <- "Blue207"
Blue208_16S_18S_scc_contrib$sample <- "Blue208"
Blue209_16S_18S_scc_contrib$sample <- "Blue209"


combined_16S_18S_scc_contrib <- rbind(Blue190_16S_18S_scc_contrib,
                                      Blue191_16S_18S_scc_contrib,
                                      Blue192_16S_18S_scc_contrib,
                                      Blue193_16S_18S_scc_contrib,
                                      Blue195_16S_18S_scc_contrib,
                                      Blue197_16S_18S_scc_contrib,
                                      Blue198_16S_18S_scc_contrib,
                                      Blue200_16S_18S_scc_contrib,
                                      Blue202_16S_18S_scc_contrib,
                                      Blue203_16S_18S_scc_contrib,
                                      Blue204_16S_18S_scc_contrib,
                                      Blue207_16S_18S_scc_contrib,
                                      Blue208_16S_18S_scc_contrib,
                                      Blue209_16S_18S_scc_contrib)

saveRDS(object = combined_16S_18S_scc_contrib,
        file="/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/old/old_blueberry_18S/blueberry_ec_16S_18S_contributions.rds")



######## TESTING ########
combined_16S_18S_scc_contrib_mean <- aggregate(.~contrib + sample, data=combined_16S_18S_scc_contrib, mean)

combined_16S_18S_scc_contrib_mean_max <- data.frame(matrix(NA, nrow=length(unique(combined_16S_18S_scc_contrib_mean$sample)), ncol=4))
rownames(combined_16S_18S_scc_contrib_mean_max) <- unique(combined_16S_18S_scc_contrib_mean$sample)
colnames(combined_16S_18S_scc_contrib_mean_max) <- c("rho_contrib", "precision_contrib", "recall_contrib", "F1_contrib")
for(samp in unique(combined_16S_18S_scc_contrib_mean$sample)) {
  combined_16S_18S_scc_contrib_mean_subset <- combined_16S_18S_scc_contrib_mean[combined_16S_18S_scc_contrib_mean$sample==samp,]
  max_rho_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$rho == max(combined_16S_18S_scc_contrib_mean_subset$rho)), "contrib"][1]
  max_precision_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$precision == max(combined_16S_18S_scc_contrib_mean_subset$precision)), "contrib"][1]
  max_recall_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$recall == max(combined_16S_18S_scc_contrib_mean_subset$recall)), "contrib"][1]
  max_F1_contrib <- combined_16S_18S_scc_contrib_mean_subset[which(combined_16S_18S_scc_contrib_mean_subset$F1 == max(combined_16S_18S_scc_contrib_mean_subset$F1)), "contrib"][1]
  combined_16S_18S_scc_contrib_mean_max[samp, ] <- c(max_rho_contrib, max_precision_contrib, max_recall_contrib, max_F1_contrib)
}


combined_16S_18S_scc_contrib_mean_max_percent <- combined_16S_18S_scc_contrib_mean_max/100

kingdom_percent_tmp <- kingdom_percent
rownames(kingdom_percent_tmp) <- gsub("BB", "Blue", rownames(kingdom_percent_tmp) )
kingdom_percent_tmp <- kingdom_percent_tmp[rownames(combined_16S_18S_scc_contrib_mean_max),]

par(mfrow=c(2,2))
plot(combined_16S_18S_scc_contrib_mean_max_percent$rho_contrib, kingdom_percent_tmp$eukaryota)
plot(combined_16S_18S_scc_contrib_mean_max_percent$recall_contrib, kingdom_percent_tmp$eukaryota)
plot(combined_16S_18S_scc_contrib_mean_max_percent$precision_contrib, kingdom_percent_tmp$eukaryota)
plot(combined_16S_18S_scc_contrib_mean_max_percent$F1_contrib  , kingdom_percent_tmp$eukaryota)

xyplot(rho~contrib|sample, data=combined_16S_18S_scc_contrib_mean)

best_contrib_scc <- c(Blue190_16S_18S_scc_contrib$contrib[which(Blue190_16S_18S_scc_contrib$rho==max(Blue190_16S_18S_scc_contrib$rho))],
                      Blue191_16S_18S_scc_contrib$contrib[which(Blue191_16S_18S_scc_contrib$rho==max(Blue191_16S_18S_scc_contrib$rho))],
                      Blue192_16S_18S_scc_contrib$contrib[which(Blue192_16S_18S_scc_contrib$rho==max(Blue192_16S_18S_scc_contrib$rho))],
                      Blue193_16S_18S_scc_contrib$contrib[which(Blue193_16S_18S_scc_contrib$rho==max(Blue193_16S_18S_scc_contrib$rho))],
                      Blue195_16S_18S_scc_contrib$contrib[which(Blue195_16S_18S_scc_contrib$rho==max(Blue195_16S_18S_scc_contrib$rho))],
                      Blue197_16S_18S_scc_contrib$contrib[which(Blue197_16S_18S_scc_contrib$rho==max(Blue197_16S_18S_scc_contrib$rho))],
                      Blue198_16S_18S_scc_contrib$contrib[which(Blue198_16S_18S_scc_contrib$rho==max(Blue198_16S_18S_scc_contrib$rho))],
                      Blue200_16S_18S_scc_contrib$contrib[which(Blue200_16S_18S_scc_contrib$rho==max(Blue200_16S_18S_scc_contrib$rho))],
                      Blue202_16S_18S_scc_contrib$contrib[which(Blue202_16S_18S_scc_contrib$rho==max(Blue202_16S_18S_scc_contrib$rho))],
                      Blue203_16S_18S_scc_contrib$contrib[which(Blue203_16S_18S_scc_contrib$rho==max(Blue203_16S_18S_scc_contrib$rho))],
                      Blue204_16S_18S_scc_contrib$contrib[which(Blue204_16S_18S_scc_contrib$rho==max(Blue204_16S_18S_scc_contrib$rho))],
                      Blue207_16S_18S_scc_contrib$contrib[which(Blue207_16S_18S_scc_contrib$rho==max(Blue207_16S_18S_scc_contrib$rho))],
                      Blue208_16S_18S_scc_contrib$contrib[which(Blue208_16S_18S_scc_contrib$rho==max(Blue208_16S_18S_scc_contrib$rho))],
                      Blue209_16S_18S_scc_contrib$contrib[which(Blue209_16S_18S_scc_contrib$rho==max(Blue209_16S_18S_scc_contrib$rho))])

kingdom_percent_subset <- kingdom_percent[-which(rownames(kingdom_percent) == "BB3"),]
kingdom_percent_subset$best_contrib_scc <- best_contrib_scc


par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(Blue191_16S_18S_scc_contrib$contrib, Blue191_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue191")
plot(Blue192_16S_18S_scc_contrib$contrib, Blue192_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue192")
plot(Blue193_16S_18S_scc_contrib$contrib, Blue193_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue193")
plot(Blue195_16S_18S_scc_contrib$contrib, Blue195_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue195")
plot(Blue197_16S_18S_scc_contrib$contrib, Blue197_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue197")
plot(Blue190_16S_18S_scc_contrib$contrib, Blue190_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue190")
plot(Blue198_16S_18S_scc_contrib$contrib, Blue198_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue198")
plot(Blue200_16S_18S_scc_contrib$contrib, Blue200_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue200")
plot(Blue202_16S_18S_scc_contrib$contrib, Blue202_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue202")
plot(Blue203_16S_18S_scc_contrib$contrib, Blue203_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue203")
plot(Blue204_16S_18S_scc_contrib$contrib, Blue204_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue204")
plot(Blue207_16S_18S_scc_contrib$contrib, Blue207_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue207")
plot(Blue208_16S_18S_scc_contrib$contrib, Blue208_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue208")
plot(Blue209_16S_18S_scc_contrib$contrib, Blue209_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="Blue209")
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(Blue190_16S_18S_scc_contrib$contrib, Blue190_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue190")
plot(Blue191_16S_18S_scc_contrib$contrib, Blue191_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue191")
plot(Blue192_16S_18S_scc_contrib$contrib, Blue192_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue192")
plot(Blue193_16S_18S_scc_contrib$contrib, Blue193_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue193")
plot(Blue200_16S_18S_scc_contrib$contrib, Blue200_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue200")
plot(Blue202_16S_18S_scc_contrib$contrib, Blue202_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue202")
plot(Blue203_16S_18S_scc_contrib$contrib, Blue203_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue203")
plot(Blue204_16S_18S_scc_contrib$contrib, Blue204_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue204")
plot(Blue195_16S_18S_scc_contrib$contrib, Blue195_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue195")
plot(Blue197_16S_18S_scc_contrib$contrib, Blue197_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue197")
plot(Blue198_16S_18S_scc_contrib$contrib, Blue198_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue198")
plot(Blue207_16S_18S_scc_contrib$contrib, Blue207_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue207")
plot(Blue208_16S_18S_scc_contrib$contrib, Blue208_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue208")
plot(Blue209_16S_18S_scc_contrib$contrib, Blue209_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="Blue209")
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(Blue190_16S_18S_scc_contrib$contrib, Blue190_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue190")
plot(Blue191_16S_18S_scc_contrib$contrib, Blue191_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue191")
plot(Blue192_16S_18S_scc_contrib$contrib, Blue192_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue192")
plot(Blue193_16S_18S_scc_contrib$contrib, Blue193_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue193")
plot(Blue200_16S_18S_scc_contrib$contrib, Blue200_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue200")
plot(Blue202_16S_18S_scc_contrib$contrib, Blue202_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue202")
plot(Blue203_16S_18S_scc_contrib$contrib, Blue203_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue203")
plot(Blue204_16S_18S_scc_contrib$contrib, Blue204_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue204")
plot(Blue195_16S_18S_scc_contrib$contrib, Blue195_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue195")
plot(Blue197_16S_18S_scc_contrib$contrib, Blue197_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue197")
plot(Blue198_16S_18S_scc_contrib$contrib, Blue198_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue198")
plot(Blue207_16S_18S_scc_contrib$contrib, Blue207_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue207")
plot(Blue208_16S_18S_scc_contrib$contrib, Blue208_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue208")
plot(Blue209_16S_18S_scc_contrib$contrib, Blue209_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue209")
par(mar=c(5.1, 4.1, 4.1, 2.1))


par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(Blue190_16S_18S_scc_contrib$contrib, Blue190_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue190")
plot(Blue191_16S_18S_scc_contrib$contrib, Blue191_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue191")
plot(Blue192_16S_18S_scc_contrib$contrib, Blue192_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue192")
plot(Blue193_16S_18S_scc_contrib$contrib, Blue193_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue193")
plot(Blue200_16S_18S_scc_contrib$contrib, Blue200_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue200")
plot(Blue202_16S_18S_scc_contrib$contrib, Blue202_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue202")
plot(Blue203_16S_18S_scc_contrib$contrib, Blue203_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue203")
plot(Blue204_16S_18S_scc_contrib$contrib, Blue204_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue204")
plot(Blue195_16S_18S_scc_contrib$contrib, Blue195_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue195")
plot(Blue197_16S_18S_scc_contrib$contrib, Blue197_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue197")
plot(Blue198_16S_18S_scc_contrib$contrib, Blue198_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue198")
plot(Blue207_16S_18S_scc_contrib$contrib, Blue207_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue207")
plot(Blue208_16S_18S_scc_contrib$contrib, Blue208_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue208")
plot(Blue209_16S_18S_scc_contrib$contrib, Blue209_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="Blue209")
par(mar=c(5.1, 4.1, 4.1, 2.1))



par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(Blue191_16S_18S_scc_contrib$contrib, Blue191_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue191")
plot(Blue192_16S_18S_scc_contrib$contrib, Blue192_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue192")
plot(Blue193_16S_18S_scc_contrib$contrib, Blue193_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue193")
plot(Blue195_16S_18S_scc_contrib$contrib, Blue195_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue195")
plot(Blue197_16S_18S_scc_contrib$contrib, Blue197_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue197")
plot(Blue190_16S_18S_scc_contrib$contrib, Blue190_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue190")
plot(Blue198_16S_18S_scc_contrib$contrib, Blue198_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue198")
plot(Blue200_16S_18S_scc_contrib$contrib, Blue200_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue200")
plot(Blue202_16S_18S_scc_contrib$contrib, Blue202_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue202")
plot(Blue203_16S_18S_scc_contrib$contrib, Blue203_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue203")
plot(Blue204_16S_18S_scc_contrib$contrib, Blue204_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue204")
plot(Blue207_16S_18S_scc_contrib$contrib, Blue207_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue207")
plot(Blue208_16S_18S_scc_contrib$contrib, Blue208_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue208")
plot(Blue209_16S_18S_scc_contrib$contrib, Blue209_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="Blue209")
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Test:

blue_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec, tab = ec_blue_mgs_nomiss_subset)
blue_ec_mgs_null_df_round <- round(blue_ec_mgs_null_df  - 0.00000001)

blue_ec_mgs_null <- cor_all_cols(tab1 = blue_ec_mgs_null_df, tab2 = ec_blue_mgs_nomiss_subset, cat_string="Null", metric="spearman")
blue_ec_mgs_16S <- cor_all_cols(tab1 = ec_blue_16S_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="16S", metric="spearman")
blue_ec_mgs_18S <- cor_all_cols(tab1 = ec_blue_18S_nomiss_subset, tab2 = ec_blue_mgs_nomiss_subset, cat_string="18S", metric="spearman")

blue_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, blue_ec_mgs_null_df, category="Null")
blue_ec_mgs_16S_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_16S_nomiss_subset, category="16S")
blue_ec_mgs_18S_metrics <- calc_accuracy_metrics(ec_blue_mgs_nomiss_subset, ec_blue_18S_nomiss_subset, category="18S")


boxplot(blue_ec_mgs_null$metric, blue_ec_mgs_16S$metric, blue_ec_mgs_18S$metric)

boxplot(blue_ec_mgs_null_metrics$precision, blue_ec_mgs_16S_metrics$precision, blue_ec_mgs_18S_metrics$precision)

boxplot(blue_ec_mgs_null_metrics$recall, blue_ec_mgs_16S_metrics$recall, blue_ec_mgs_18S_metrics$recall)


# Comparing UniRef90 and UniRef50-based EC databases
ec50 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/func_tables/protozoa_fungi_faa_uniref50_hits_level4ec_regrouped.txt",
                   header=T, sep="\t", row.names=1)

ec90 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/func_tables/protozoa_fungi_faa_uniref90_hits_level4ec_regrouped.txt",
                   header=T, sep="\t", row.names=1)

overlapping_ec_uniref <- rownames(ec50)[which(rownames(ec50) %in% rownames(ec90))]
ec50 <- ec50[overlapping_ec_uniref,]
ec90 <- ec90[overlapping_ec_uniref,colnames(ec50)]

cor.test(rowSums(ec50), rowSums(ec90), method="spearman")
