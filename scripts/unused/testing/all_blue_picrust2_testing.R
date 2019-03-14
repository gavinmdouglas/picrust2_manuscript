rm(list=ls(all=TRUE))

setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S_all_samples")

library(vegan)
library(randomForest)
library(pROC)

blue_map <- read.table("18S_map_1.txt", header=TRUE, sep="\t", comment.char="", stringsAsFactors = FALSE)
rownames(blue_map) <- gsub(".assembled.filtered.nonchimera.fasta", "", blue_map$X.SampleID)
rownames(blue_map) <- gsub("\\.", "-", rownames(blue_map))
blue_root_managed_samples <- rownames(blue_map)[which(blue_map$Description_Managemet == "MngRoot")]
blue_rhizo_managed_samples <- rownames(blue_map)[which(blue_map$Description_Managemet == "MngRhizo")]
blue_bulk_managed_samples <- rownames(blue_map)[which(blue_map$Description_Managemet == "MngBulk")]

managed_rhizo_bulk_samples <- c(blue_rhizo_managed_samples, blue_bulk_managed_samples)
managed_samples <- c(blue_root_managed_samples, blue_rhizo_managed_samples, blue_bulk_managed_samples)

blue_asv <- read.table("deblur_output_exported/blueberry_18S_all.biom.tsv",
                      header=T, sep="\t", row.names=1, check.names=FALSE, skip=1, comment.char = "")

fungi_asvs <- read.table("fungi_asvs.txt", header=F, stringsAsFactors = FALSE)$V1

blue_asv <- blue_asv[fungi_asvs, ]

blue_asv_t <- data.frame(t(blue_asv), check.names = FALSE)

blue_asv_t <- blue_asv_t[managed_samples[which(managed_samples %in% rownames(blue_asv_t))], ]

blue_asv_t_relab_managed <- data.frame(sweep(blue_asv_t, 1, rowSums(blue_asv_t), '/'), check.names = FALSE) * 100

blue_asv_t_relab_managed_dist <- vegdist(blue_asv_t_relab_managed, method = "bray")
blue_asv_t_relab_managed_NMDS <- metaMDS(blue_asv_t_relab_managed_dist, k=20, try=100, trymax=1000)

blue_asv_t_relab_managed_NMDS_df <- data.frame(NMDS1=blue_asv_t_relab_managed_NMDS$points[,1],
                                               NMDS2=blue_asv_t_relab_managed_NMDS$points[,2],
                                               Group=as.factor(blue_map[rownames(blue_asv_t_relab_managed_NMDS$points), "Description_Managemet"]))

plot(blue_asv_t_relab_managed_NMDS_df$NMDS1,
     blue_asv_t_relab_managed_NMDS_df$NMDS2,
     col=blue_asv_t_relab_managed_NMDS_df$Group,
     pch=16)

blue_ec <- read.table("picrust2_full_output/ec_18S_counts_metagenome_out/pred_metagenome_unstrat.tsv",
                      header=T, sep="\t", row.names=1, check.names=FALSE)

blue_ec_t <- data.frame(t(blue_ec), check.names=FALSE)

all_ec <- colnames(blue_ec_t)
blue_ec_t <- data.frame(t(blue_ec), check.names = FALSE)

blue_ec_t <- blue_ec_t[managed_samples[which(managed_samples %in% rownames(blue_ec_t))], ]

blue_ec_t_relab_managed <- data.frame(sweep(blue_ec_t, 1, rowSums(blue_ec_t), '/'), check.names = FALSE) * 100

blue_ec_t_relab_managed_dist <- vegdist(blue_ec_t_relab_managed, method = "bray")
blue_ec_t_relab_managed_NMDS <- metaMDS(blue_ec_t_relab_managed_dist, k=20, try=100, trymax=1000)

blue_ec_t_relab_managed_NMDS_df <- data.frame(NMDS1=blue_ec_t_relab_managed_NMDS$points[,1],
                                               NMDS2=blue_ec_t_relab_managed_NMDS$points[,2],
                                               Group=as.factor(blue_map[rownames(blue_ec_t_relab_managed_NMDS$points), "Description_Managemet"]))



plot(blue_ec_t_relab_managed_NMDS_df$NMDS1,
     blue_ec_t_relab_managed_NMDS_df$NMDS2,
     col=blue_ec_t_relab_managed_NMDS_df$Group,
     pch=16)


asv_NMDS_points <- data.frame(blue_asv_t_relab_managed_NMDS$points, check.names = FALSE)
ec_NMDS_points <- data.frame(blue_ec_t_relab_managed_NMDS$points, check.names = FALSE)

colnames(asv_NMDS_points) <- gsub("^", "ASV_", colnames(asv_NMDS_points))
colnames(ec_NMDS_points) <- gsub("^", "EC_", colnames(ec_NMDS_points))

asv_ec_NMDS_points <- cbind(asv_NMDS_points, ec_NMDS_points)

asv_ec_NMDS_points$group <- blue_map[rownames(asv_ec_NMDS_points), "Description_Managemet"]
asv_ec_NMDS_points$group <- factor(asv_ec_NMDS_points$group)

asv_ec_NMDS_points_RF <- randomForest(x=asv_ec_NMDS_points[,1:(ncol(asv_ec_NMDS_points)-1)],
                                               y=asv_ec_NMDS_points[ , ncol(asv_ec_NMDS_points)],
                                               ntree=501, importance=TRUE, proximities=TRUE )

asv_ec_NMDS_points_RF.roc <- multiclass.roc(asv_ec_NMDS_points$group, asv_ec_NMDS_points_RF$votes[,2])
asv_ec_NMDS_points_RF_auc <- auc(asv_ec_NMDS_points_RF.roc)


ec_NMDS_points$group <- blue_map[rownames(ec_NMDS_points), "Description_Managemet"]
ec_NMDS_points$group <- factor(ec_NMDS_points$group)

ec_NMDS_points_RF <- randomForest(x=ec_NMDS_points[,1:(ncol(ec_NMDS_points)-1)],
                                      y=ec_NMDS_points[ , ncol(ec_NMDS_points)],
                                      ntree=501, importance=TRUE, proximities=TRUE )

asv_NMDS_points_RF.roc <- multiclass.roc(asv_NMDS_points$group, asv_NMDS_points_RF$votes[,2])
asv_NMDS_points_RF_auc <- auc(asv_NMDS_points_RF.roc)

asv_NMDS_points$group <- blue_map[rownames(asv_NMDS_points), "Description_Managemet"]
asv_NMDS_points$group <- factor(asv_NMDS_points$group)

asv_NMDS_points_RF <- randomForest(x=asv_NMDS_points[,1:(ncol(asv_NMDS_points)-1)],
                                   y=asv_NMDS_points[ , ncol(asv_NMDS_points)],
                                   ntree=501, importance=TRUE, proximities=TRUE)

asv_NMDS_points_RF.roc <- multiclass.roc(asv_NMDS_points$group, asv_NMDS_points_RF$votes[,2])
asv_NMDS_points_RF_auc <- auc(asv_NMDS_points_RF.roc)


blue_asv_ec_t_relab <- cbind(blue_asv_t_relab, blue_ec_t_relab[, top_ec])
blue_asv_ec_t_relab$group <- blue_map[rownames(blue_asv_ec_t_relab), "Description_Managemet"]
blue_asv_ec_t_relab_managed <- blue_asv_ec_t_relab[which(blue_asv_ec_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]
blue_asv_ec_t_relab_managed$group <- factor(blue_asv_ec_t_relab_managed$group)

blue_asv_ec_t_relab_managed_RF <- randomForest(x=blue_asv_ec_t_relab_managed[,1:(ncol(blue_asv_ec_t_relab_managed)-1)],
                                            y=blue_asv_ec_t_relab_managed[ , ncol(blue_asv_ec_t_relab_managed)],
                                            ntree=501, importance=TRUE, proximities=TRUE )

blue_asv_ec_t_relab_managed_RF_imp <- data.frame(blue_asv_ec_t_relab_managed_RF$importance)

blue_asv_ec_t_relab_managed_RF_imp_ordered <- blue_asv_ec_t_relab_managed_RF_imp[with(blue_asv_ec_t_relab_managed_RF_imp, order(-MeanDecreaseAccuracy)),]

func = "X8cdb5f5b4e1db57ab763c75a69a4eb85"
boxplot(blue_asv_ec_t_relab_managed[blue_bulk_managed_samples, func],
        blue_asv_ec_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_asv_ec_t_relab_managed[blue_root_managed_samples, func])

blue_asv_ec_t_relab_managed_RF.roc <- multiclass.roc(blue_asv_ec_t_relab_managed$group, blue_asv_ec_t_relab_managed_RF$votes[,2])
blue_asv_ec_t_relab_managed_RF_auc <- auc(blue_asv_ec_t_relab_managed_RF.roc)



blue_asv_t_relab$group <- blue_map[rownames(blue_asv_t_relab), "Description_Managemet"]
blue_asv_t_relab_managed <- blue_asv_t_relab[which(blue_asv_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]

blue_asv_t_relab_managed$group <- factor(blue_asv_t_relab_managed$group)

blue_asv_t_relab_managed_RF <- randomForest(x=blue_asv_t_relab_managed[,1:(ncol(blue_asv_t_relab_managed)-1)],
                                           y=blue_asv_t_relab_managed[ , ncol(blue_asv_t_relab_managed)],
                                           ntree=501, importance=TRUE, proximities=TRUE )

blue_asv_t_relab_managed_RF_imp <- data.frame(blue_asv_t_relab_managed_RF$importance)

blue_asv_t_relab_managed_RF_imp_ordered <- blue_asv_t_relab_managed_RF_imp[with(blue_asv_t_relab_managed_RF_imp, order(-MeanDecreaseAccuracy)),]

func = "X8cdb5f5b4e1db57ab763c75a69a4eb85"
boxplot(blue_asv_t_relab_managed[blue_bulk_managed_samples, func],
        blue_asv_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_asv_t_relab_managed[blue_root_managed_samples, func])

blue_asv_t_relab_managed_RF.roc <- multiclass.roc(blue_asv_t_relab_managed$group, blue_asv_t_relab_managed_RF$votes[,2])
blue_asv_t_relab_managed_RF_auc <- auc(blue_asv_t_relab_managed_RF.roc)



blue_ec <- read.table("picrust2_full_output/ec_18S_counts_metagenome_out/pred_metagenome_unstrat.tsv",
                      header=T, sep="\t", row.names=1, check.names=FALSE)


blue_ec_t <- data.frame(t(blue_ec), check.names=FALSE)

all_ec <- colnames(blue_ec_t)

blue_ec_t_relab <- data.frame(sweep(blue_ec_t, 1, rowSums(blue_ec_t), '/'), check.names = FALSE) * 100


blue_ec_t_relab$group <- blue_map[rownames(blue_ec_t_relab), "Description_Managemet"]
blue_ec_t_relab_managed <- blue_ec_t_relab[which(blue_ec_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]

blue_ec_t_relab_managed$group <- factor(blue_ec_t_relab_managed$group)

blue_ec_t_relab_managed_RF <- randomForest(x=blue_ec_t_relab_managed[,1:(ncol(blue_ec_t_relab_managed)-1)],
                                            y=blue_ec_t_relab_managed[ , ncol(blue_ec_t_relab_managed)],
                                           ntree=501, importance=TRUE, proximities=TRUE )

blue_ec_t_relab_managed_RF_imp <- data.frame(blue_ec_t_relab_managed_RF$importance)

blue_ec_t_relab_managed_RF_imp_ordered <- blue_ec_t_relab_managed_RF_imp[with(blue_ec_t_relab_managed_RF_imp, order(-MeanDecreaseAccuracy)),]

func = "EC:1.13.11.33"
boxplot(blue_ec_t_relab_managed[blue_bulk_managed_samples, func],
        blue_ec_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_ec_t_relab_managed[blue_root_managed_samples, func])




length(which(p.adjust(pvalues, "bonferroni") < 0.05))

which(is.na(pvalues))

blue_pathabun_t_relab$group <- blue_map[rownames(blue_pathabun_t_relab), "Description_Managemet"]
blue_pathabun_t_relab_managed <- blue_pathabun_t_relab[which(blue_pathabun_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]

blue_pathabun_t_relab_managed$group <- factor(blue_pathabun_t_relab_managed$group)

blue_pathabun_t_relab_managed_RF <- randomForest(x=blue_pathabun_t_relab_managed[,1:(ncol(blue_pathabun_t_relab_managed)-1)],
                                           y=blue_pathabun_t_relab_managed[ , ncol(blue_pathabun_t_relab_managed)],
                                           ntree=5001, importance=TRUE, proximities=TRUE )

blue_pathabun_t_relab_managed_RF_imp <- data.frame(blue_pathabun_t_relab_managed_RF$importance)

blue_pathabun_t_relab_managed_RF_imp_ordered <- blue_pathabun_t_relab_managed_RF_imp[with(blue_pathabun_t_relab_managed_RF_imp, order(-MeanDecreaseAccuracy)),]

func = "PWY-6609"
boxplot(blue_pathabun_t_relab_managed[blue_bulk_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_root_managed_samples, func])

func = "PWY66-422"
boxplot(blue_pathabun_t_relab_managed[blue_bulk_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_root_managed_samples, func])

func = "PWY-6317"
boxplot(blue_pathabun_t_relab_managed[blue_bulk_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_root_managed_samples, func])

func = "PWY-5083"
boxplot(blue_pathabun_t_relab_managed[blue_bulk_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_root_managed_samples, func])





blue_asv_pathabun_t_relab <- cbind(blue_asv_t_relab, blue_pathabun_t_relab)
blue_asv_pathabun_t_relab$group <- blue_map[rownames(blue_asv_pathabun_t_relab), "Description_Managemet"]
blue_asv_pathabun_t_relab_managed <- blue_asv_pathabun_t_relab[which(blue_asv_pathabun_t_relab$group %in% c("MngRoot", "MngRhizo", "MngBulk")), ]
blue_asv_pathabun_t_relab_managed$group <- factor(blue_asv_pathabun_t_relab_managed$group)

blue_asv_pathabun_t_relab_managed_RF <- randomForest(x=blue_asv_pathabun_t_relab_managed[,1:(ncol(blue_asv_pathabun_t_relab_managed)-1)],
                                               y=blue_asv_pathabun_t_relab_managed[ , ncol(blue_asv_pathabun_t_relab_managed)],
                                               ntree=501, importance=TRUE, proximities=TRUE )

blue_asv_pathabun_t_relab_managed_RF_imp <- data.frame(blue_asv_pathabun_t_relab_managed_RF$importance)

blue_asv_pathabun_t_relab_managed_RF_imp_ordered <- blue_asv_pathabun_t_relab_managed_RF_imp[with(blue_asv_pathabun_t_relab_managed_RF_imp, order(-MeanDecreaseAccuracy)),]

blue_asv_pathabun_t_relab_managed_RF.roc <- multiclass.roc(blue_asv_pathabun_t_relab_managed$group, blue_asv_pathabun_t_relab_managed_RF$votes[,2])
blue_asv_pathabun_t_relab_managed_RF_auc <- auc(blue_asv_pathabun_t_relab_managed_RF.roc)


func = "X8cdb5f5b4e1db57ab763c75a69a4eb85"
boxplot(blue_asv_pathabun_t_relab_managed[blue_bulk_managed_samples, func],
        blue_asv_pathabun_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_asv_pathabun_t_relab_managed[blue_root_managed_samples, func])

blue_asv_pathabun_t_relab_managed_RF.roc <- multiclass.roc(blue_asv_pathabun_t_relab_managed$group, blue_asv_pathabun_t_relab_managed_RF$votes[,2])
blue_asv_pathabun_t_relab_managed_RF_auc <- auc(blue_asv_pathabun_t_relab_managed_RF.roc)








blue_ec <- read.table("picrust2_full_output/ec_18S_counts_metagenome_out/pred_metagenome_unstrat.tsv",
                            header=T, sep="\t", row.names=1, check.names=FALSE)


blue_ec_t <- data.frame(t(blue_ec), check.names=FALSE)

all_ec <- colnames(blue_ec_t)

blue_ec_t_relab <- data.frame(sweep(blue_ec_t, 1, rowSums(blue_ec_t), '/'), check.names = FALSE) * 100
blue_ec_t_relab_managed <- blue_ec_t_relab[managed_samples[which(managed_samples %in% rownames(blue_asv_t))], ]

blue_ec_t_relab_managed <- blue_ec_t_relab_managed[, -which(colSums(blue_ec_t_relab_managed) == 0)]

all_path <- colnames(blue_ec_t_relab_managed)

blue_ec_t_relab_managed$group <- factor(blue_map[rownames(blue_ec_t_relab_managed), "Description_Managemet"])


pvalues <- c()

for(path in all_path) {
  
  subset_df <- blue_ec_t_relab_managed[,c(path, "group")]
  
  colnames(subset_df) <- c("path", "group")
  
  tmp <- kruskal.test(formula=path ~ group, data=subset_df)
  
  pvalues <- c(pvalues, tmp$p.value)
  
}

pvalues_bonf <- p.adjust(pvalues, "bonferroni")

















blue_pathabun <- read.table("picrust2_full_output/pathways_out/path_abun_unstrat.tsv",
                            header=T, sep="\t", row.names=1, check.names=FALSE)


blue_pathabun_t <- data.frame(t(blue_pathabun), check.names=FALSE)

blue_pathabun_t_relab <- data.frame(sweep(blue_pathabun_t, 1, rowSums(blue_pathabun_t), '/'), check.names = FALSE) * 100
blue_pathabun_t_relab_managed <- blue_pathabun_t_relab[managed_samples[which(managed_samples %in% rownames(blue_pathabun_t_relab))], ]

blue_pathabun_t_relab_managed <- blue_pathabun_t_relab_managed[, -which(colSums(blue_pathabun_t_relab_managed) == 0)]

all_path <- colnames(blue_pathabun_t_relab_managed)

blue_pathabun_t_relab_managed$group <- factor(blue_map[rownames(blue_pathabun_t_relab_managed), "Description_Managemet"])


pvalues <- c()

for(path in all_path) {
  
  subset_df <- blue_pathabun_t_relab_managed[,c(path, "group")]
  
  colnames(subset_df) <- c("path", "group")
  
  tmp <- kruskal.test(formula=path ~ group, data=subset_df)
  
  pvalues <- c(pvalues, tmp$p.value)
  
}

pvalues_bonf <- p.adjust(pvalues, "bonf")

# [1] "ANAGLYCOLYSIS-PWY" "LEU-DEG2-PWY"      "NONOXIPENT-PWY"    "PWY-5083"          "PWY-6317"          "PWY-6609"         
# [7] "PWY-7385"          "SO4ASSIM-PWY" 

func <- "SO4ASSIM-PWY"

boxplot(blue_pathabun_t_relab_managed[blue_bulk_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_rhizo_managed_samples, func],
        blue_pathabun_t_relab_managed[blue_root_managed_samples, func])


