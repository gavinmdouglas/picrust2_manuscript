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


blue_pathabun <- read.table("picrust2_full_output/pathways_out/path_abun_unstrat.tsv",
                      header=T, sep="\t", row.names=1, check.names=FALSE)


blue_pathabun_t <- data.frame(t(blue_pathabun), check.names=FALSE)

all_pathabun <- colnames(blue_pathabun_t)

blue_pathabun_t_relab <- data.frame(sweep(blue_pathabun_t, 1, rowSums(blue_pathabun_t), '/'), check.names = FALSE) * 100


blue_pathabun_t_relab$group <- blue_map[rownames(blue_pathabun_t_relab), "Description_Managemet"]
blue_pathabun_t_relab_managed <- blue_pathabun_t_relab[which(blue_pathabun_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]

blue_pathabun_t_relab_managed$group <- factor(blue_pathabun_t_relab_managed$group)

blue_pathabun_t_relab_managed_RF <- rfPermute(x=blue_pathabun_t_relab_managed[,1:(ncol(blue_pathabun_t_relab_managed)-1)],
                                        y=blue_pathabun_t_relab_managed[ , ncol(blue_pathabun_t_relab_managed)],
                                        ntree=501, importance=TRUE, proximities=TRUE, nrep=1000, num.cores =20)

sig_pathabuns_i <- as.numeric(which(blue_pathabun_t_relab_managed_RF$pval[,, "scaled"][,"MeanDecreaseAccuracy"] < 0.001))

sig_pathabuns <- names(blue_pathabun_t_relab_managed_RF$pval[,, "scaled"][sig_pathabuns_i,"MeanDecreaseAccuracy"])

blue_pathabun_t_relab_managed_sig_subset <- blue_pathabun_t_relab_managed[, c(sig_pathabuns, "group")]

blue_pathabun_t_relab_managed_sig_subset$sample <- rownames(blue_pathabun_t_relab_managed_sig_subset)

write.table(x = blue_pathabun_t_relab_managed_sig_subset, file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/all_blueberry_sig_pathways.tsv",
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
