### Commands to identify significant predicted pathways between bulk, rhizophere, and blueberry root samples in managed fields.
### Write out table of relative abundances of significant pathways for subsequent plotting.

rm(list=ls(all=TRUE))

setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S_all_samples")

library(randomForest)
library(rfPermute)
library(pROC)

blue_map <- read.table("18S_map_1.txt", header=TRUE, sep="\t", comment.char="", stringsAsFactors = FALSE)
rownames(blue_map) <- gsub(".assembled.filtered.nonchimera.fasta", "", blue_map$X.SampleID)
rownames(blue_map) <- gsub("\\.", "-", rownames(blue_map))
blue_root_managed_samples <- rownames(blue_map)[which(blue_map$Description_Managemet == "MngRoot")]
blue_rhizo_managed_samples <- rownames(blue_map)[which(blue_map$Description_Managemet == "MngRhizo")]
blue_bulk_managed_samples <- rownames(blue_map)[which(blue_map$Description_Managemet == "MngBulk")]

managed_samples <- c(blue_root_managed_samples, blue_rhizo_managed_samples, blue_bulk_managed_samples)

# Run RF on ASVs:
blue_asvs <- read.table("deblur_output_exported/blueberry_18S_all.biom.tsv",
                            header=T, sep="\t", row.names=1, check.names=FALSE, comment.char="", skip=1)

blue_asvs_t <- data.frame(t(blue_asvs), check.names=FALSE)

all_asvs <- colnames(blue_asvs_t)

blue_asvs_t_relab <- data.frame(sweep(blue_asvs_t, 1, rowSums(blue_asvs_t), '/'), check.names = FALSE) * 100


blue_asvs_t_relab$group <- blue_map[rownames(blue_asvs_t_relab), "Description_Managemet"]
blue_asvs_t_relab_managed <- blue_asvs_t_relab[which(blue_asvs_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]

blue_asvs_t_relab_managed$group <- factor(blue_asvs_t_relab_managed$group)

# Shouldn't be re-run (see RDS below).
# blue_asvs_t_relab_managed_RF <- rfPermute(x=blue_asvs_t_relab_managed[,1:(ncol(blue_asvs_t_relab_managed)-1)],
#                                               y=blue_asvs_t_relab_managed[ , ncol(blue_asvs_t_relab_managed)],
#                                               ntree=501, importance=TRUE, proximities=TRUE, nrep=1000, num.cores =20)
# 
# saveRDS(object = blue_asvs_t_relab_managed_RF, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/blueberry_RF/blue_asvs_t_relab_managed_RF.rds")
blue_asvs_t_relab_managed_RF <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/blueberry_RF/blue_asvs_t_relab_managed_RF.rds")

# Run RF on pathway abundances:
blue_ec <- read.table("picrust2_full_output/ec_18S_counts_metagenome_out/pred_metagenome_unstrat.tsv",
                      header=T, sep="\t", row.names=1, check.names=FALSE)

blue_ec_t <- data.frame(t(blue_ec), check.names=FALSE)

all_ec <- colnames(blue_ec_t)

blue_ec_t_relab <- data.frame(sweep(blue_ec_t, 1, rowSums(blue_ec_t), '/'), check.names = FALSE) * 100


blue_ec_t_relab$group <- blue_map[rownames(blue_ec_t_relab), "Description_Managemet"]
blue_ec_t_relab_managed <- blue_ec_t_relab[which(blue_ec_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]

blue_ec_t_relab_managed$group <- factor(blue_ec_t_relab_managed$group)

# Shouldn't be re-run (see RDS below).
# blue_ec_t_relab_managed_RF <- rfPermute(x=blue_ec_t_relab_managed[,1:(ncol(blue_ec_t_relab_managed)-1)],
#                                         y=blue_ec_t_relab_managed[ , ncol(blue_ec_t_relab_managed)],
#                                         ntree=501, importance=TRUE, proximities=TRUE, nrep=1000, num.cores =20)
# 
# saveRDS(object = blue_ec_t_relab_managed_RF, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/blueberry_RF/blue_ec_t_relab_managed_RF.rds")
blue_ec_t_relab_managed_RF <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/saved_RDS/blueberry_RF/blue_ec_t_relab_managed_RF.rds")

sig_ecs_i <- as.numeric(which(blue_ec_t_relab_managed_RF$pval[,, "scaled"][,"MeanDecreaseAccuracy"] < 0.001))

sig_ecs <- names(blue_ec_t_relab_managed_RF$pval[,, "scaled"][sig_ecs_i,"MeanDecreaseAccuracy"])

blue_ec_t_relab_managed_sig_subset <- blue_ec_t_relab_managed[, c(sig_ecs, "group")]

blue_ec_t_relab_managed_sig_subset$sample <- rownames(blue_ec_t_relab_managed_sig_subset)

write.table(x = blue_ec_t_relab_managed_sig_subset, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/all_blueberry_sig_ecs.tsv",
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
