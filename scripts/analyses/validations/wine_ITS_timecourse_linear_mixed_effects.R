library(lme4)
library(blme)
library(vegan)
library(ggplot2)

### Prep metadata and subset of samples for input to q2-longitudinal.

wine_picrust2_outdir <- "/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/"

wine_picrust2_ec_unstrat_file <- paste(wine_picrust2_outdir, "EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv", sep="")

wine_picrust2_pathabun_unstrat_file <- paste(wine_picrust2_outdir, "pathways_out_nsti_2.0/path_abun_unstrat.tsv", sep="")

wine_picrust2_ec_unstrat <- read.table(wine_picrust2_ec_unstrat_file, header=TRUE, sep="\t", check.names=FALSE, row.names=1)

wine_picrust2_pathabun_unstrat <- read.table(wine_picrust2_pathabun_unstrat_file, header=TRUE, sep="\t", check.names=FALSE, row.names=1)

wine_ITS_abun <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/2014_swarm.2.otutable.txt",
                            header=T, sep=" ", check.names=FALSE, stringsAsFactors = FALSE, row.names=1)

wine_ITS_tax <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/2014_swarm.2.otutable.taxonomy.txt",
                            header=T, sep="\t", check.names=FALSE, stringsAsFactors = FALSE)

rownames(wine_ITS_tax) <- wine_ITS_tax$otu_id

rownames(wine_ITS_abun) <- gsub(";size.*$", "", rownames(wine_ITS_abun))
rownames(wine_ITS_tax) <- gsub(";size.*$", "", rownames(wine_ITS_tax))

wine_ITS_tax$genus <- gsub("^.*g__", "g__", wine_ITS_tax$`QIIME_taxonomy_augmented(98%)`)
wine_ITS_tax$genus <- gsub(";.*$", "", wine_ITS_tax$genus)

# Remove control samples as well as treasurt and samples starting with "14":
col2keep <- c(grep("^Y._", colnames(wine_picrust2_ec_unstrat), value=TRUE),
              grep("^T._", colnames(wine_picrust2_ec_unstrat), value=TRUE))

# Remove Y1_D4-A, Y1_D0-B, Y1_D0-A, and Y1_D0-C due to low depth. Next lowest sample had 64,139 reads.
col2remove <- c("Y1_D4-A", "Y1_D0-A", "Y1_D0-B", "Y1_D0-C")

col2keep <- col2keep[-which(col2keep %in% col2remove)]

# Read in predicted pathway abundances per sequence.
ec_pred <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/EC_predicted.tsv",
                               header=T, sep="\t", row.names=1, check.names = FALSE)
ec_pred <- data.frame(t(round(ec_pred)), check.names=FALSE)
yeast_asvs <- rownames(wine_ITS_tax)[which(wine_ITS_tax$genus=="g__Saccharomyces")]
ec_pred_yeast <- ec_pred[, yeast_asvs]
ec_pred_other <- ec_pred[, -which(colnames(ec_pred) %in% yeast_asvs)]


fisher_out <- c()

for(i in 1:nrow(ec_pred_yeast)) {
  
  yeast_sum_path <- sum(ec_pred_yeast[i,])
  all_yeast <- sum(rowSums(ec_pred_yeast))
  
  other_sum_path <- sum(ec_pred_other[i,])
  all_other <- sum(rowSums(ec_pred_other))
  
  in_matrix <- matrix(c(yeast_sum_path, all_yeast, other_sum_path, all_other), nrow=2, ncol=2)
  
  fisher_test_out <- fisher.test(in_matrix, alternative="greater")
  
  fisher_out <- c(fisher_out, fisher_test_out$p.value)
}

fisher_out_bonf <- p.adjust(fisher_out, method="bonferroni")

yeast_enriched_ec <- rownames(ec_pred_yeast)[which(fisher_out_bonf < 0.05)]

ec_pred_yeast_means <- rowMeans(ec_pred_yeast[yeast_enriched_ec,])



per_seq_ec <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/wine_fungi/ITS/picrust2_pipeline/picrust2_full_output/ec_ITS_counts_metagenome_out/pred_metagenome_strat.tsv",
                         header=T, sep="\t", check.names = FALSE)
colnames(per_seq_ec)[1] <- "func"

per_seq_ec_yeast <- per_seq_ec[which(per_seq_ec$sequence %in% yeast_asvs) , ]
per_seq_ec_other <- per_seq_ec[-which(per_seq_ec$sequence %in% yeast_asvs) , ]

per_seq_ec_yeast <- per_seq_ec_yeast[, -which(colnames(per_seq_ec_yeast) == "sequence")]
per_seq_ec_other <- per_seq_ec_other[, -which(colnames(per_seq_ec_other) == "sequence")]

per_seq_ec_yeast_sum <- aggregate(. ~ func, data=per_seq_ec_yeast, FUN=sum)
per_seq_ec_other_sum <- aggregate(. ~ func, data=per_seq_ec_other, FUN=sum)

rownames(per_seq_ec_yeast_sum) <- per_seq_ec_yeast_sum$func
rownames(per_seq_ec_other_sum) <- per_seq_ec_other_sum$func

per_seq_ec_yeast_sum <- per_seq_ec_yeast_sum[, col2keep]
per_seq_ec_other_sum <- per_seq_ec_other_sum[, col2keep]

per_seq_ec_yeast_sum_relab <- data.frame(sweep(per_seq_ec_yeast_sum, 2, colSums(per_seq_ec_yeast_sum), `/`), check.names=FALSE)
per_seq_ec_other_sum_relab <- data.frame(sweep(per_seq_ec_other_sum, 2, colSums(per_seq_ec_other_sum), `/`), check.names=FALSE)


ec_pred_yeast_sums <- ec_pred_yeast_sums[which(names(ec_pred_yeast_sums) %in% rownames(per_seq_ec_other_sum))]

spearman_r <- c()
for(sample in colnames(per_seq_ec_other_sum_relab)) {
  spearman_r <- c(spearman_r, cor.test(per_seq_ec_other_sum_relab[, sample], ec_pred_yeast_means, method="spearman")$estimate)
}
names(spearman_r) <- colnames(per_seq_ec_other_sum_relab)

spearman_r_D0 <- spearman_r[grep("D0", names(spearman_r))]

yeast_D1_relab <- wine_ITS_abun_subset_genus_relab[grep("D1", rownames(wine_ITS_abun_subset_genus_relab)), "g__Saccharomyces", drop=FALSE]

per_seq_ec_other_sum_relab_t <- data.frame(t(per_seq_ec_other_sum_relab), check.names=TRUE)
per_seq_ec_other_sum_relab_t <- log(per_seq_ec_other_sum_relab_t + 0.00001)
per_seq_ec_other_sum_relab_t$sample <- gsub("_.*$", "", rownames(per_seq_ec_other_sum_relab_t))
per_seq_ec_other_sum_relab_t$stage <- gsub("^.*D", "", rownames(per_seq_ec_other_sum_relab_t))
per_seq_ec_other_sum_relab_t$stage <- as.numeric(gsub("-.*$", "", per_seq_ec_other_sum_relab_t$stage))
per_seq_ec_other_sum_relab_t$replicate <- factor(gsub("^.*-", "", rownames(per_seq_ec_other_sum_relab_t)))

per_seq_ec_other_sum_relab_t_D0 <- per_seq_ec_other_sum_relab_t[which(per_seq_ec_other_sum_relab_t$stage == 0),]

rownames(yeast_D1_relab) <- gsub("D1", "D0", rownames(yeast_D1_relab))

per_seq_ec_other_sum_relab_t_D0 <- cbind(per_seq_ec_other_sum_relab_t_D0, yeast_D1_relab[rownames(per_seq_ec_other_sum_relab_t_D0),, drop=FALSE])

per_seq_ec_other_sum_relab_t_D0$sample <- factor(per_seq_ec_other_sum_relab_t_D0$sample)

tmp <- lmer(g__Saccharomyces ~ EC.1.1.1.169 + (1|sample), data=per_seq_ec_other_sum_relab_t_D0, REML=FALSE)
tmp2 <- lmer(g__Saccharomyces ~  (1|sample), data=per_seq_ec_other_sum_relab_t_D0, REML=FALSE)


p_values <- c()

for(enriched_ec in grep("EC", colnames(per_seq_ec_other_sum_relab_t_D0), value=TRUE)) {
  
  f1 <- as.formula(paste("g__Saccharomyces", "~", "(1 | sample)", sep = " "))
  f2 <- as.formula(paste("g__Saccharomyces", "~", enriched_ec, "+ (1 | sample)", sep = " "))
  
  tmp_null <- lmer(f1, data=per_seq_ec_other_sum_relab_t_D0, REML=FALSE)
  
  tmp_actual <- lmer(f2, data=per_seq_ec_other_sum_relab_t_D0, REML=FALSE)
  
  tmp3 <- anova(tmp_null, tmp_actual)
  
  p_values <- c(p_values, tmp3$`Pr(>Chisq)`[2])
}









wine_picrust2_ec_unstrat_subset <- wine_picrust2_ec_unstrat[, col2keep]

wine_picrust2_pathabun_unstrat_subset <- wine_picrust2_pathabun_unstrat[, col2keep]

wine_ITS_abun_subset <- wine_ITS_abun[, col2keep]

wine_ITS_abun_subset$genus <- factor(wine_ITS_tax[rownames(wine_ITS_abun_subset), "genus"])

wine_ITS_abun_subset_genus <- aggregate(. ~ genus, data=wine_ITS_abun_subset, FUN=sum)
rownames(wine_ITS_abun_subset_genus) <- wine_ITS_abun_subset_genus$genus
wine_ITS_abun_subset_genus <- wine_ITS_abun_subset_genus[, -which(colnames(wine_ITS_abun_subset_genus) == "genus")]
wine_ITS_abun_subset_genus_relab <- data.frame(t(sweep(wine_ITS_abun_subset_genus, 2, colSums(wine_ITS_abun_subset_genus), `/`) * 100), check.names=FALSE)
wine_ITS_abun_subset_genus_relab_dist <- vegdist(wine_ITS_abun_subset_genus_relab, method = "bray")
wine_ITS_abun_subset_genus_relab_dist_NMDS <- metaMDS(wine_ITS_abun_subset_genus_relab_dist)


wine_ITS_abun_subset_genus_relab$sample <- gsub("_.*$", "", rownames(wine_ITS_abun_subset_genus_relab))
wine_ITS_abun_subset_genus_relab$stage <- gsub("^.*D", "", rownames(wine_ITS_abun_subset_genus_relab))
wine_ITS_abun_subset_genus_relab$stage <- as.numeric(gsub("-.*$", "", wine_ITS_abun_subset_genus_relab$stage))

wine_ITS_abun_subset_genus_relab_dist_NMDS_df <- data.frame(NMDS1=wine_ITS_abun_subset_genus_relab_dist_NMDS$points[,1],
                                                            NMDS2=wine_ITS_abun_subset_genus_relab_dist_NMDS$points[,2],
                                                            Sample=as.factor(wine_ITS_abun_subset_genus_relab[rownames(wine_ITS_abun_subset_genus_relab_dist_NMDS$points), "sample"]),
                                                            Stage=as.factor(wine_ITS_abun_subset_genus_relab[rownames(wine_ITS_abun_subset_genus_relab_dist_NMDS$points), "stage"]))

tol15rainbow_modified=c("white", "black", "#77AADD", "#117755", "#44AA88")

ggplot(data=wine_ITS_abun_subset_genus_relab_dist_NMDS_df, aes(NMDS1, NMDS2, colour=Sample)) +
  geom_point(aes(fill=Sample, size=1.5, shape=Stage)) +
  theme_minimal() +
  scale_fill_manual(values=tol15rainbow_modified) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=6)))


boxplot(wine_ITS_abun_subset_genus_relab$g__Saccharomyces ~ wine_ITS_abun_subset_genus_relab$stage)


# Conver to relab.
wine_picrust2_ec_unstrat_subset_relab <- data.frame(sweep(wine_picrust2_ec_unstrat_subset, 2, colSums(wine_picrust2_ec_unstrat_subset), `/`) * 100, check.names=FALSE)
wine_picrust2_ec_unstrat_subset_relab_t <- data.frame(t(wine_picrust2_ec_unstrat_subset_relab), check.names = FALSE)
wine_picrust2_ec_unstrat_subset_relab_t_dist <- vegdist(wine_picrust2_ec_unstrat_subset_relab_t, method = "bray")
wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS <- metaMDS(wine_picrust2_ec_unstrat_subset_relab_t_dist)


wine_picrust2_ec_unstrat_subset_relab_t$sample <- gsub("_.*$", "", rownames(wine_picrust2_ec_unstrat_subset_relab_t))
wine_picrust2_ec_unstrat_subset_relab_t$stage <- gsub("^.*D", "", rownames(wine_picrust2_ec_unstrat_subset_relab_t))
wine_picrust2_ec_unstrat_subset_relab_t$stage <- as.numeric(gsub("-.*$", "", wine_picrust2_ec_unstrat_subset_relab_t$stage))

wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS_df <- data.frame(NMDS1=wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS$points[,1],
                                                                         NMDS2=wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS$points[,2],
                                                                         Sample=as.factor(wine_picrust2_ec_unstrat_subset_relab_t[rownames(wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS$points), "sample"]),
                                                                         Stage=as.factor(wine_picrust2_ec_unstrat_subset_relab_t[rownames(wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS$points), "stage"]))

tol15rainbow_modified=c("white", "black", "#77AADD", "#117755", "#44AA88")

ggplot(data=wine_picrust2_ec_unstrat_subset_relab_t_dist_NMDS_df, aes(NMDS1, NMDS2, colour=Sample)) +
  geom_point(aes(fill=Sample, size=1.5, shape=Stage)) +
  theme_minimal() +
  scale_fill_manual(values=tol15rainbow_modified) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=6)))



wine_picrust2_pathabun_unstrat_subset_relab <- data.frame(sweep(wine_picrust2_pathabun_unstrat_subset, 2, colSums(wine_picrust2_pathabun_unstrat_subset), `/`) * 100, check.names=FALSE)

wine_picrust2_pathabun_unstrat_subset_relab_t <- data.frame(t(wine_picrust2_pathabun_unstrat_subset_relab), check.names = FALSE)
wine_picrust2_pathabun_unstrat_subset_relab_t_dist <- vegdist(wine_picrust2_pathabun_unstrat_subset_relab_t, method = "bray")
wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS <- metaMDS(wine_picrust2_pathabun_unstrat_subset_relab_t_dist)


wine_picrust2_pathabun_unstrat_subset_relab_t$sample <- gsub("_.*$", "", rownames(wine_picrust2_pathabun_unstrat_subset_relab_t))
wine_picrust2_pathabun_unstrat_subset_relab_t$stage <- gsub("^.*D", "", rownames(wine_picrust2_pathabun_unstrat_subset_relab_t))
wine_picrust2_pathabun_unstrat_subset_relab_t$stage <- as.numeric(gsub("-.*$", "", wine_picrust2_pathabun_unstrat_subset_relab_t$stage))

wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS_df <- data.frame(NMDS1=wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS$points[,1],
                                                            NMDS2=wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS$points[,2],
                                                            Sample=as.factor(wine_picrust2_pathabun_unstrat_subset_relab_t[rownames(wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS$points), "sample"]),
                                                            Stage=as.factor(wine_picrust2_pathabun_unstrat_subset_relab_t[rownames(wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS$points), "stage"]))

tol15rainbow_modified=c("white", "black", "#77AADD", "#117755", "#44AA88")

ggplot(data=wine_picrust2_pathabun_unstrat_subset_relab_t_dist_NMDS_df, aes(NMDS1, NMDS2, colour=Sample)) +
  geom_point(aes(fill=Sample, size=1.5, shape=Stage)) +
  theme_minimal() +
  scale_fill_manual(values=tol15rainbow_modified) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=6)))


# Remove rare features (in fewer than 15% of samples, i.e. 8.4 samples).
wine_picrust2_ec_unstrat_subset_relab_filt <- wine_picrust2_ec_unstrat_subset_relab[-which(rowSums(wine_picrust2_ec_unstrat_subset_relab > 0) < 0.15*56),]

wine_picrust2_pathabun_unstrat_subset_relab_filt <- wine_picrust2_pathabun_unstrat_subset_relab[-which(rowSums(wine_picrust2_pathabun_unstrat_subset_relab > 0) < 0.15*56),]


wine_ec_prepped <- data.frame(t(wine_picrust2_ec_unstrat_subset_relab_filt))
ec2test <- colnames(wine_ec_prepped)
wine_ec_prepped$sample <- gsub("_.*$", "", rownames(wine_ec_prepped))
wine_ec_prepped$stage <- gsub("^.*D", "", rownames(wine_ec_prepped))
wine_ec_prepped$stage <- as.numeric(gsub("-.*$", "", wine_ec_prepped$stage))
# Remove stage 4:
wine_ec_prepped <- wine_ec_prepped[-which(wine_ec_prepped$stage == 4),]
wine_ec_prepped$replicate <- factor(gsub("^.*-", "", rownames(wine_ec_prepped)))

wine_ec_prepped_baseline <- wine_ec_prepped[which(wine_ec_prepped$stage==0),]
wine_ec_prepped_baseline$sample <- as.factor(wine_ec_prepped_baseline$sample)
wine_ec_prepped$sample <- as.factor(wine_ec_prepped$sample)

wine_pathabun_prepped <- data.frame(t(wine_picrust2_pathabun_unstrat_subset_relab_filt))
path2test <- colnames(wine_pathabun_prepped)
wine_pathabun_prepped$sample <- gsub("_.*$", "", rownames(wine_pathabun_prepped))
wine_pathabun_prepped$stage <- gsub("^.*D", "", rownames(wine_pathabun_prepped))
wine_pathabun_prepped$stage <- as.numeric(gsub("-.*$", "", wine_pathabun_prepped$stage))
# Remove stage 4:
wine_pathabun_prepped <- wine_pathabun_prepped[-which(wine_pathabun_prepped$stage == 4),]
wine_pathabun_prepped$replicate <- factor(gsub("^.*-", "", rownames(wine_pathabun_prepped)))


wine_pathabun_prepped_baseline <- wine_pathabun_prepped[which(wine_pathabun_prepped$stage==0),]
wine_pathabun_prepped_baseline$sample <- as.factor(wine_pathabun_prepped_baseline$sample)
wine_pathabun_prepped$sample <- as.factor(wine_pathabun_prepped$sample)

tmp <- lmer(EC.1.1.1.1 ~ sample + (1 | replicate), data=wine_ec_prepped_baseline, REML=FALSE)
boxplot(wine_ec_prepped_baseline$EC.1.1.1.1 ~ wine_ec_prepped_baseline$sample)
tmp <- kruskal.test(EC.1.1.1.101 ~ sample, data=wine_ec_prepped_baseline)
summary(tmp)[[1]][["Pr(>F)"]][1]


tmp2 <- lm(EC.1.1.1.1 ~ stage + sample + sample:replicate, data=wine_ec_prepped)

tmp <- lmer(EC.1.1.1.1 ~ stage + (1 | sample), data=wine_ec_prepped, REML=FALSE)
tmp2 <- lmer(EC.1.1.1.1 ~ (1 | sample), data=wine_ec_prepped, REML=FALSE)
tmp2 <- blmer(EC.1.1.1.1 ~ (1 | sample) + (1 | sample : replicate), data=wine_ec_prepped, REML=FALSE)
tmp <- blmer(EC.1.1.1.1 ~ stage + (1 | sample), data=wine_ec_prepped, REML=FALSE)

anova(tmp2, tmp)

tmp2 <- lm(`EC:1.1.1.1` ~ stage + sample + replicate,  data=wine_ec_prepped)

summary(tmp)

p_values <- c()
for(path in path2test) {
  in_formula <- as.formula(paste(path, "~", "sample", sep = " "))
  kruskal_out <- kruskal.test(in_formula, data=wine_pathabun_prepped_baseline)
  p_values <- c(p_values, kruskal_out$p.value)
}

p_values <- c()
for(path in path2test) {

  f1 <- as.formula(paste(path, "~", "(1 | sample) + (1 | sample : replicate)", sep = " "))
  f2 <- as.formula(paste(path, "~", "stage + (1 | sample) + (1 | sample : replicate)", sep = " "))
  
  tmp_null <- blmer(f1, data=wine_pathabun_prepped, REML=FALSE)
  
  tmp_actual <- blmer(f2, data=wine_pathabun_prepped, REML=FALSE)
  
  tmp3 <- anova(tmp_null, tmp_actual)
  
  p_values <- c(p_values, tmp3$`Pr(>Chisq)`[2])
}



