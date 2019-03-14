setwd("/home/gavin/projects/picrust_pipeline/data/validation/eukarya/JGI_SoC/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

ko_soc_16S_nsti2 <- read.table("metaxa2_out/picrust2_full_output_pro/KO_metagenome_out/pred_metagenome_unstrat.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ko_soc_18S_nsti2 <- read.table("metaxa2_out/picrust2_full_output_euk/ko_18S_metagenome_out/pred_metagenome_unstrat.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ko_soc_16S_nsti2 <- ko_soc_16S_nsti2[,-1]


ko_soc_mgs <- read.table("humann2_final_out/humann2_KO_unstratified.tsv",
                          header=T, sep="\t", row.names=1, comment.char="", quote="")
colnames(ko_soc_mgs) <- gsub("_Abundance.RPKs", "", colnames(ko_soc_mgs))

tmp1 <- cor_all_cols(ko_soc_16S_nsti2, ko_soc_mgs, metric = "spearman")
tmp2 <- cor_all_cols(ko_soc_18S_nsti2, ko_soc_mgs, metric = "spearman")

tmp_overlap <- rownames(ko_soc_18S_nsti2)[rownames(ko_soc_18S_nsti2) %in% rownames(ko_soc_mgs)]
ko_soc_16S_nsti2_sub <- ko_soc_16S_nsti2[tmp_overlap,]
ko_soc_mgs_sub <- ko_soc_mgs[tmp_overlap,]
ko_soc_18S_nsti2_sub <- ko_soc_18S_nsti2[tmp_overlap,]
plot(ko_soc_mgs_sub$ARC5Mpool, ko_soc_16S_nsti2_sub$ARC5Mpool)

plot(ko_soc_mgs_sub$ARC5Mpool, ko_soc_18S_nsti2_sub$ARC5Mpool)



ec_soc_16S_nsti2 <- read.table("metaxa2_out/picrust2_full_output_pro/EC_metagenome_out///pred_metagenome_unstrat.tsv",
                               header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti0.05 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_0.05/pred_metagenome_unstrat.tsv",
                               header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti0.1 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_0.1/pred_metagenome_unstrat.tsv",
                                  header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti0.25 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_0.25/pred_metagenome_unstrat.tsv",
                                  header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti0.5 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_0.5/pred_metagenome_unstrat.tsv",
                                  header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti1 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_1/pred_metagenome_unstrat.tsv",
                                 header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti1.5 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_1.5/pred_metagenome_unstrat.tsv",
                                 header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_18S_nsti2 <- read.table("metaxa2_out/picrust2_full_output_euk/EC_metagenome_out_nsti_2/pred_metagenome_unstrat.tsv",
                                 header=T, sep="\t", row.names=1, comment.char="", quote="", check.names=FALSE)

ec_soc_16S_nsti2 <- ec_soc_16S_nsti2[,-1]


ec_soc_mgs <- read.table("humann2_final_out/humann2_level4ec_unstratified.tsv",
                         header=T, sep="\t", row.names=1, comment.char="", quote="")
colnames(ec_soc_mgs) <- gsub("_Abundance.RPKs", "", colnames(ec_soc_mgs))
rownames(ec_soc_mgs) <- gsub("^", "EC:", rownames(ec_soc_mgs))

tmp1 <- cor_all_cols(ec_soc_16S_nsti2, ec_soc_mgs, metric = "spearman")
tmp2 <- cor_all_cols(ec_soc_18S_nsti2, ec_soc_mgs, metric = "spearman")

tmp_overlap <- rownames(ec_soc_18S_nsti2)[rownames(ec_soc_18S_nsti2) %in% rownames(ec_soc_mgs)]
ec_soc_16S_nsti2_sub <- ec_soc_16S_nsti2[tmp_overlap,]
ec_soc_mgs_sub <- ec_soc_mgs[tmp_overlap,]
ec_soc_18S_nsti2_sub <- ec_soc_18S_nsti2[tmp_overlap,]
plot(ec_soc_mgs_sub$ARC5Mpool, ec_soc_16S_nsti2_sub$ARC5Mpool)

plot(ec_soc_mgs_sub$ARC5Mpool, ec_soc_18S_nsti2_sub$ARC5Mpool)

# Add in ECs called as 0.
possible_mgs_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/humann2_ECs.txt")$V1
possible_16S_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/picrust2_ECs.txt")$V1
possible_18S_ecs <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_ECs/picrust2_ECs_18S.txt")$V1
possible_picrust2_ecs <- factor(unique(c(levels(possible_18S_ecs), levels(possible_16S_ecs))))

ec_soc_mgs <- ec_soc_mgs[-which(rownames(ec_soc_mgs) %in% c("EC:UNMAPPED", "EC:UNGROUPED")),]

# Subset to ECs overlapping across all 3 tables.
ec_soc_mgs_nomiss <- add_missing_funcs(ec_soc_mgs, possible_mgs_ecs)
ec_soc_16S_nsti2_nomiss <- add_missing_funcs(ec_soc_16S_nsti2, possible_16S_ecs)
ec_soc_18S_nsti0.05_nomiss <- add_missing_funcs(ec_soc_18S_nsti0.05, possible_18S_ecs)
ec_soc_18S_nsti0.1_nomiss <- add_missing_funcs(ec_soc_18S_nsti0.1, possible_18S_ecs)
ec_soc_18S_nsti0.25_nomiss <- add_missing_funcs(ec_soc_18S_nsti0.25, possible_18S_ecs)
ec_soc_18S_nsti0.5_nomiss <- add_missing_funcs(ec_soc_18S_nsti0.5, possible_18S_ecs)
ec_soc_18S_nsti1_nomiss <- add_missing_funcs(ec_soc_18S_nsti1, possible_18S_ecs)
ec_soc_18S_nsti1.5_nomiss <- add_missing_funcs(ec_soc_18S_nsti1.5, possible_18S_ecs)
ec_soc_18S_nsti2_nomiss <- add_missing_funcs(ec_soc_18S_nsti2, possible_18S_ecs)

overlapping_ec <- rownames(ec_soc_16S_nsti2_nomiss)[which(rownames(ec_soc_16S_nsti2_nomiss) %in% rownames(ec_soc_mgs_nomiss))]
overlapping_ec <- overlapping_ec[which(overlapping_ec %in% rownames(ec_soc_18S_nsti2_nomiss))]

overlapping_col <- colnames(ec_soc_16S_nsti2_nomiss)[which(colnames(ec_soc_16S_nsti2_nomiss) %in% colnames(ec_soc_mgs_nomiss))]
overlapping_col <- overlapping_col[which(overlapping_col %in% colnames(ec_soc_18S_nsti2_nomiss))]

ec_soc_mgs_nomiss_subset <- ec_soc_mgs_nomiss[overlapping_ec, overlapping_col]
ec_soc_16S_nsti2_nomiss_subset <- ec_soc_16S_nsti2_nomiss[overlapping_ec, overlapping_col]


# Comparing UniRef90 and UniRef50-based EC databases
ec_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/mean_func_tables/ec_18S.txt",
                   header=T, sep="\t", row.names=1, check.names=FALSE)

soc_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec_18S, tab = ec_soc_mgs_nomiss_subset)
soc_ec_mgs_null_df_round <- round(soc_ec_mgs_null_df  - 0.00000001)

soc_ec_mgs_null <- cor_all_cols(tab1 = soc_ec_mgs_null_df, tab2 = ec_soc_mgs_nomiss_subset, cat_string="Null", metric="spearman")
soc_ec_mgs_16S <- cor_all_cols(tab1 = ec_soc_16S_nsti2_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="16S", metric="spearman")

soc_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, soc_ec_mgs_null_df_round, category="Null")
soc_ec_mgs_16S_metrics <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_16S_nsti2_nomiss_subset, category="16S")

ec_soc_18S_nsti2_nomiss_subset <- ec_soc_18S_nsti2_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti2 <- cor_all_cols(tab1 = ec_soc_18S_nsti2_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti2 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti2_nomiss_subset, category="18S")

ec_soc_18S_nsti1.5_nomiss_subset <- ec_soc_18S_nsti1.5_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti1.5 <- cor_all_cols(tab1 = ec_soc_18S_nsti1.5_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti1.5 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti1.5_nomiss_subset, category="18S")


ec_soc_18S_nsti1_nomiss_subset <- ec_soc_18S_nsti1_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti1 <- cor_all_cols(tab1 = ec_soc_18S_nsti1_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti1 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti1_nomiss_subset, category="18S")


ec_soc_18S_nsti0.5_nomiss_subset <- ec_soc_18S_nsti0.5_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti0.5 <- cor_all_cols(tab1 = ec_soc_18S_nsti0.5_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti0.5 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti0.5_nomiss_subset, category="18S")


ec_soc_18S_nsti0.25_nomiss_subset <- ec_soc_18S_nsti0.25_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti0.25 <- cor_all_cols(tab1 = ec_soc_18S_nsti0.25_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti0.25 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti0.25_nomiss_subset, category="18S")


ec_soc_18S_nsti0.1_nomiss_subset <- ec_soc_18S_nsti0.1_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti0.1 <- cor_all_cols(tab1 = ec_soc_18S_nsti0.1_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti0.1 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti0.1_nomiss_subset, category="18S")


ec_soc_18S_nsti0.05_nomiss_subset <- ec_soc_18S_nsti0.05_nomiss[, overlapping_col]
soc_ec_mgs_18S_nsti0.05 <- cor_all_cols(tab1 = ec_soc_18S_nsti0.05_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")
soc_ec_mgs_18S_metrics_nsti0.05 <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nsti0.05_nomiss_subset, category="18S")


par(mfrow=c(3,1))
boxplot(soc_ec_mgs_null_metrics$precision, soc_ec_mgs_18S_metrics$precision, soc_ec_mgs_16S_metrics$precision)
boxplot(soc_ec_mgs_null_metrics$recall, soc_ec_mgs_18S_metrics$recall, soc_ec_mgs_16S_metrics$recall)
boxplot(soc_ec_mgs_null$metric, soc_ec_mgs_18S$metric, soc_ec_mgs_16S$metric)

par(mfrow=c(3,1))
boxplot(soc_ec_mgs_null$metric, soc_ec_mgs_18S_nsti2$metric, soc_ec_mgs_18S_nsti1.5$metric, soc_ec_mgs_18S_nsti1$metric,
        soc_ec_mgs_18S_nsti0.5$metric, soc_ec_mgs_18S_nsti0.25$metric, soc_ec_mgs_18S_nsti0.1$metric, soc_ec_mgs_18S_nsti0.05$metric, soc_ec_mgs_16S$metric)

boxplot(soc_ec_mgs_null_metrics$recall, soc_ec_mgs_18S_metrics_nsti2$recall, soc_ec_mgs_18S_metrics_nsti1.5$recall, soc_ec_mgs_18S_metrics_nsti1$recall,
        soc_ec_mgs_18S_metrics_nsti0.5$recall, soc_ec_mgs_18S_metrics_nsti0.25$recall, soc_ec_mgs_18S_metrics_nsti0.1$recall, soc_ec_mgs_18S_metrics_nsti0.05$recall,soc_ec_mgs_16S_metrics$recall)

boxplot(soc_ec_mgs_null_metrics$precision, soc_ec_mgs_18S_metrics_nsti2$precision, soc_ec_mgs_18S_metrics_nsti1.5$precision, soc_ec_mgs_18S_metrics_nsti1$precision,
        soc_ec_mgs_18S_metrics_nsti0.5$precision, soc_ec_mgs_18S_metrics_nsti0.25$precision, soc_ec_mgs_18S_metrics_nsti0.1$precision, soc_ec_mgs_18S_metrics_nsti0.05$precision,soc_ec_mgs_16S_metrics$precision)


# Determine if a combination of 16S and 18S predicitons does better than 16S alone.
ec_soc_16S_relabun <-  data.frame(sweep(ec_soc_16S_nsti2_nomiss_subset, 2, colSums(ec_soc_16S_nsti2_nomiss_subset), '/'))
ec_soc_18S_relabun <-  data.frame(sweep(ec_soc_18S_nsti2_nomiss_subset, 2, colSums(ec_soc_18S_nsti2_nomiss_subset), '/'))

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


SRR5819848_16S_18S_scc_contrib <- combine_predictions(ec_soc_16S_relabun$SRR5819848, ec_soc_18S_relabun$SRR5819848, ec_soc_mgs_nomiss_subset$SRR5819848)
SRR5819969_16S_18S_scc_contrib <- combine_predictions(ec_soc_16S_relabun$SRR5819969, ec_soc_18S_relabun$SRR5819969, ec_soc_mgs_nomiss_subset$SRR5819969)
SRR5819968_16S_18S_scc_contrib <- combine_predictions(ec_soc_16S_relabun$SRR5819968, ec_soc_18S_relabun$SRR5819968, ec_soc_mgs_nomiss_subset$SRR5819968)
ARC7Mpool_16S_18S_scc_contrib <- combine_predictions(ec_soc_16S_relabun$ARC7Mpool, ec_soc_18S_relabun$ARC7Mpool, ec_soc_mgs_nomiss_subset$ARC7Mpool)
SRR5819818_16S_18S_scc_contrib <- combine_predictions(ec_soc_16S_relabun$SRR5819818, ec_soc_18S_relabun$SRR5819818, ec_soc_mgs_nomiss_subset$SRR5819818)
ARC5Mpool_16S_18S_scc_contrib <- combine_predictions(ec_soc_16S_relabun$ARC5Mpool, ec_soc_18S_relabun$ARC5Mpool, ec_soc_mgs_nomiss_subset$ARC5Mpool)


soc190_16S_18S_scc_contrib$sample <- "soc190"
soc191_16S_18S_scc_contrib$sample <- "soc191"
soc192_16S_18S_scc_contrib$sample <- "soc192"
soc193_16S_18S_scc_contrib$sample <- "soc193"
soc195_16S_18S_scc_contrib$sample <- "soc195"
soc197_16S_18S_scc_contrib$sample <- "soc197"
soc198_16S_18S_scc_contrib$sample <- "soc198"
soc200_16S_18S_scc_contrib$sample <- "soc200"
soc202_16S_18S_scc_contrib$sample <- "soc202"
soc203_16S_18S_scc_contrib$sample <- "soc203"
soc204_16S_18S_scc_contrib$sample <- "soc204"
soc207_16S_18S_scc_contrib$sample <- "soc207"
soc208_16S_18S_scc_contrib$sample <- "soc208"
soc209_16S_18S_scc_contrib$sample <- "soc209"


combined_16S_18S_scc_contrib <- rbind(soc190_16S_18S_scc_contrib,
                                      soc191_16S_18S_scc_contrib,
                                      soc192_16S_18S_scc_contrib,
                                      soc193_16S_18S_scc_contrib,
                                      soc195_16S_18S_scc_contrib,
                                      soc197_16S_18S_scc_contrib,
                                      soc198_16S_18S_scc_contrib,
                                      soc200_16S_18S_scc_contrib,
                                      soc202_16S_18S_scc_contrib,
                                      soc203_16S_18S_scc_contrib,
                                      soc204_16S_18S_scc_contrib,
                                      soc207_16S_18S_scc_contrib,
                                      soc208_16S_18S_scc_contrib,
                                      soc209_16S_18S_scc_contrib)

saveRDS(object = combined_16S_18S_scc_contrib,
        file="/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/socberry_ec_16S_18S_contributions.rds")




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
rownames(kingdom_percent_tmp) <- gsub("BB", "soc", rownames(kingdom_percent_tmp) )
kingdom_percent_tmp <- kingdom_percent_tmp[rownames(combined_16S_18S_scc_contrib_mean_max),]

par(mfrow=c(2,2))
plot(combined_16S_18S_scc_contrib_mean_max_percent$rho_contrib, kingdom_percent_tmp$eukaryota)
plot(combined_16S_18S_scc_contrib_mean_max_percent$recall_contrib, kingdom_percent_tmp$eukaryota)
plot(combined_16S_18S_scc_contrib_mean_max_percent$precision_contrib, kingdom_percent_tmp$eukaryota)
plot(combined_16S_18S_scc_contrib_mean_max_percent$F1_contrib  , kingdom_percent_tmp$eukaryota)

xyplot(rho~contrib|sample, data=combined_16S_18S_scc_contrib_mean)

best_contrib_scc <- c(soc190_16S_18S_scc_contrib$contrib[which(soc190_16S_18S_scc_contrib$rho==max(soc190_16S_18S_scc_contrib$rho))],
                      soc191_16S_18S_scc_contrib$contrib[which(soc191_16S_18S_scc_contrib$rho==max(soc191_16S_18S_scc_contrib$rho))],
                      soc192_16S_18S_scc_contrib$contrib[which(soc192_16S_18S_scc_contrib$rho==max(soc192_16S_18S_scc_contrib$rho))],
                      soc193_16S_18S_scc_contrib$contrib[which(soc193_16S_18S_scc_contrib$rho==max(soc193_16S_18S_scc_contrib$rho))],
                      soc195_16S_18S_scc_contrib$contrib[which(soc195_16S_18S_scc_contrib$rho==max(soc195_16S_18S_scc_contrib$rho))],
                      soc197_16S_18S_scc_contrib$contrib[which(soc197_16S_18S_scc_contrib$rho==max(soc197_16S_18S_scc_contrib$rho))],
                      soc198_16S_18S_scc_contrib$contrib[which(soc198_16S_18S_scc_contrib$rho==max(soc198_16S_18S_scc_contrib$rho))],
                      soc200_16S_18S_scc_contrib$contrib[which(soc200_16S_18S_scc_contrib$rho==max(soc200_16S_18S_scc_contrib$rho))],
                      soc202_16S_18S_scc_contrib$contrib[which(soc202_16S_18S_scc_contrib$rho==max(soc202_16S_18S_scc_contrib$rho))],
                      soc203_16S_18S_scc_contrib$contrib[which(soc203_16S_18S_scc_contrib$rho==max(soc203_16S_18S_scc_contrib$rho))],
                      soc204_16S_18S_scc_contrib$contrib[which(soc204_16S_18S_scc_contrib$rho==max(soc204_16S_18S_scc_contrib$rho))],
                      soc207_16S_18S_scc_contrib$contrib[which(soc207_16S_18S_scc_contrib$rho==max(soc207_16S_18S_scc_contrib$rho))],
                      soc208_16S_18S_scc_contrib$contrib[which(soc208_16S_18S_scc_contrib$rho==max(soc208_16S_18S_scc_contrib$rho))],
                      soc209_16S_18S_scc_contrib$contrib[which(soc209_16S_18S_scc_contrib$rho==max(soc209_16S_18S_scc_contrib$rho))])

kingdom_percent_subset <- kingdom_percent[-which(rownames(kingdom_percent) == "BB3"),]
kingdom_percent_subset$best_contrib_scc <- best_contrib_scc


par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(soc191_16S_18S_scc_contrib$contrib, soc191_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc191")
plot(soc192_16S_18S_scc_contrib$contrib, soc192_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc192")
plot(soc193_16S_18S_scc_contrib$contrib, soc193_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc193")
plot(soc195_16S_18S_scc_contrib$contrib, soc195_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc195")
plot(soc197_16S_18S_scc_contrib$contrib, soc197_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc197")
plot(soc190_16S_18S_scc_contrib$contrib, soc190_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc190")
plot(soc198_16S_18S_scc_contrib$contrib, soc198_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc198")
plot(soc200_16S_18S_scc_contrib$contrib, soc200_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc200")
plot(soc202_16S_18S_scc_contrib$contrib, soc202_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc202")
plot(soc203_16S_18S_scc_contrib$contrib, soc203_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc203")
plot(soc204_16S_18S_scc_contrib$contrib, soc204_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc204")
plot(soc207_16S_18S_scc_contrib$contrib, soc207_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc207")
plot(soc208_16S_18S_scc_contrib$contrib, soc208_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc208")
plot(soc209_16S_18S_scc_contrib$contrib, soc209_16S_18S_scc_contrib$rho, xlab="Proportion contributed by 16S", ylab="Spearman's correlation coefficient", main="soc209")
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(soc190_16S_18S_scc_contrib$contrib, soc190_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc190")
plot(soc191_16S_18S_scc_contrib$contrib, soc191_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc191")
plot(soc192_16S_18S_scc_contrib$contrib, soc192_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc192")
plot(soc193_16S_18S_scc_contrib$contrib, soc193_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc193")
plot(soc200_16S_18S_scc_contrib$contrib, soc200_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc200")
plot(soc202_16S_18S_scc_contrib$contrib, soc202_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc202")
plot(soc203_16S_18S_scc_contrib$contrib, soc203_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc203")
plot(soc204_16S_18S_scc_contrib$contrib, soc204_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc204")
plot(soc195_16S_18S_scc_contrib$contrib, soc195_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc195")
plot(soc197_16S_18S_scc_contrib$contrib, soc197_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc197")
plot(soc198_16S_18S_scc_contrib$contrib, soc198_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc198")
plot(soc207_16S_18S_scc_contrib$contrib, soc207_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc207")
plot(soc208_16S_18S_scc_contrib$contrib, soc208_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc208")
plot(soc209_16S_18S_scc_contrib$contrib, soc209_16S_18S_scc_contrib$precision, xlab="Proportion contributed by 16S", ylab="Precision", main="soc209")
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(soc190_16S_18S_scc_contrib$contrib, soc190_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc190")
plot(soc191_16S_18S_scc_contrib$contrib, soc191_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc191")
plot(soc192_16S_18S_scc_contrib$contrib, soc192_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc192")
plot(soc193_16S_18S_scc_contrib$contrib, soc193_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc193")
plot(soc200_16S_18S_scc_contrib$contrib, soc200_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc200")
plot(soc202_16S_18S_scc_contrib$contrib, soc202_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc202")
plot(soc203_16S_18S_scc_contrib$contrib, soc203_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc203")
plot(soc204_16S_18S_scc_contrib$contrib, soc204_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc204")
plot(soc195_16S_18S_scc_contrib$contrib, soc195_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc195")
plot(soc197_16S_18S_scc_contrib$contrib, soc197_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc197")
plot(soc198_16S_18S_scc_contrib$contrib, soc198_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc198")
plot(soc207_16S_18S_scc_contrib$contrib, soc207_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc207")
plot(soc208_16S_18S_scc_contrib$contrib, soc208_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc208")
plot(soc209_16S_18S_scc_contrib$contrib, soc209_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc209")
par(mar=c(5.1, 4.1, 4.1, 2.1))


par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(soc190_16S_18S_scc_contrib$contrib, soc190_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc190")
plot(soc191_16S_18S_scc_contrib$contrib, soc191_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc191")
plot(soc192_16S_18S_scc_contrib$contrib, soc192_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc192")
plot(soc193_16S_18S_scc_contrib$contrib, soc193_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc193")
plot(soc200_16S_18S_scc_contrib$contrib, soc200_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc200")
plot(soc202_16S_18S_scc_contrib$contrib, soc202_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc202")
plot(soc203_16S_18S_scc_contrib$contrib, soc203_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc203")
plot(soc204_16S_18S_scc_contrib$contrib, soc204_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc204")
plot(soc195_16S_18S_scc_contrib$contrib, soc195_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc195")
plot(soc197_16S_18S_scc_contrib$contrib, soc197_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc197")
plot(soc198_16S_18S_scc_contrib$contrib, soc198_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc198")
plot(soc207_16S_18S_scc_contrib$contrib, soc207_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc207")
plot(soc208_16S_18S_scc_contrib$contrib, soc208_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc208")
plot(soc209_16S_18S_scc_contrib$contrib, soc209_16S_18S_scc_contrib$F1, xlab="Proportion contributed by 16S", ylab="F1", main="soc209")
par(mar=c(5.1, 4.1, 4.1, 2.1))



par(mfrow=c(4,4), mar=c(1.5,1.5,1.5,1.5))
plot(soc191_16S_18S_scc_contrib$contrib, soc191_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc191")
plot(soc192_16S_18S_scc_contrib$contrib, soc192_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc192")
plot(soc193_16S_18S_scc_contrib$contrib, soc193_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc193")
plot(soc195_16S_18S_scc_contrib$contrib, soc195_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc195")
plot(soc197_16S_18S_scc_contrib$contrib, soc197_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc197")
plot(soc190_16S_18S_scc_contrib$contrib, soc190_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc190")
plot(soc198_16S_18S_scc_contrib$contrib, soc198_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc198")
plot(soc200_16S_18S_scc_contrib$contrib, soc200_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc200")
plot(soc202_16S_18S_scc_contrib$contrib, soc202_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc202")
plot(soc203_16S_18S_scc_contrib$contrib, soc203_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc203")
plot(soc204_16S_18S_scc_contrib$contrib, soc204_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc204")
plot(soc207_16S_18S_scc_contrib$contrib, soc207_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc207")
plot(soc208_16S_18S_scc_contrib$contrib, soc208_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc208")
plot(soc209_16S_18S_scc_contrib$contrib, soc209_16S_18S_scc_contrib$recall, xlab="Proportion contributed by 16S", ylab="recall", main="soc209")
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Test:

soc_ec_mgs_null_df <- generate_null_mean_db_funcs(db = ec, tab = ec_soc_mgs_nomiss_subset)
soc_ec_mgs_null_df_round <- round(soc_ec_mgs_null_df  - 0.00000001)

soc_ec_mgs_null <- cor_all_cols(tab1 = soc_ec_mgs_null_df, tab2 = ec_soc_mgs_nomiss_subset, cat_string="Null", metric="spearman")
soc_ec_mgs_16S <- cor_all_cols(tab1 = ec_soc_16S_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="16S", metric="spearman")
soc_ec_mgs_18S <- cor_all_cols(tab1 = ec_soc_18S_nomiss_subset, tab2 = ec_soc_mgs_nomiss_subset, cat_string="18S", metric="spearman")

soc_ec_mgs_null_metrics <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, soc_ec_mgs_null_df, category="Null")
soc_ec_mgs_16S_metrics <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_16S_nomiss_subset, category="16S")
soc_ec_mgs_18S_metrics <- calc_accuracy_metrics(ec_soc_mgs_nomiss_subset, ec_soc_18S_nomiss_subset, category="18S")


boxplot(soc_ec_mgs_null$metric, soc_ec_mgs_16S$metric, soc_ec_mgs_18S$metric)

boxplot(soc_ec_mgs_null_metrics$precision, soc_ec_mgs_16S_metrics$precision, soc_ec_mgs_18S_metrics$precision)

boxplot(soc_ec_mgs_null_metrics$recall, soc_ec_mgs_16S_metrics$recall, soc_ec_mgs_18S_metrics$recall)


# Comparing UniRef90 and UniRef50-based EC databases
ec50 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/func_tables/protozoa_fungi_faa_uniref50_hits_level4ec_regrouped.txt",
                   header=T, sep="\t", row.names=1)

ec90 <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/func_tables/protozoa_fungi_faa_uniref90_hits_level4ec_regrouped.txt",
                   header=T, sep="\t", row.names=1)

overlapping_ec_uniref <- rownames(ec50)[which(rownames(ec50) %in% rownames(ec90))]
ec50 <- ec50[overlapping_ec_uniref,]
ec90 <- ec90[overlapping_ec_uniref,colnames(ec50)]

cor.test(rowSums(ec50), rowSums(ec90), method="spearman")
