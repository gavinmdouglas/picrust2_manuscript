### Test for associations between PICRUSt2 predictions and RNA-seq levels
### for certain genes of interest (i.e. those previously known to be biomarkers of CD and/or UC).

rm(list=ls(all=TRUE))

library(blme)
library(factoextra)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")

# Read in HMP2 pathway abundances and transform with arcsine square root transformation.
hmp2_pathabun <- readRDS("prepped_tables/hmp2_pathabun_filt_asin_Ileum.rds")

# Load RNA-seq data.
host_rna_rpm_ileum <- readRDS(file = "prepped_tables/hmp2_rna_rpm_log_Ileum.rds")


# Read in hmp2 metadata table.
hmp2_metadata <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_2018-08-20_col_subset.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

rnaseq_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "host_transcriptomics"), ]
rnaseq_meta <- rnaseq_meta[-which(duplicated(rnaseq_meta$`Participant ID`)), ]
rownames(rnaseq_meta) <- rnaseq_meta$`Participant ID`

diagnosis_map <- data.frame(sample=rnaseq_meta$`Participant ID`,
                            diagnosis=rnaseq_meta$diagnosis,
                            consent_age=rnaseq_meta$consent_age,
                            site_name=rnaseq_meta$site_name,
                            stringsAsFactors = FALSE)

overlapping_samples <- rownames(hmp2_pathabun)[which(rownames(hmp2_pathabun) %in% rownames(host_rna_rpm_ileum))]

diagnosis_map <- diagnosis_map[which(diagnosis_map$sample %in% overlapping_samples), ]

CD_ileum_samples <- diagnosis_map[which(diagnosis_map$diagnosis == "CD"), "sample"]

diagnosis_map_CD <- diagnosis_map[which(diagnosis_map$sample %in% CD_ileum_samples), ]

diagnosis_map_CD <- diagnosis_map_CD[which(diagnosis_map_CD$sample %in% overlapping_samples), ]

diagnosis_map_CD$site_name <- factor(diagnosis_map_CD$site_name)

hmp2_pathabun_ileum_CD <- hmp2_pathabun[CD_ileum_samples, ]

# Remove pathways found in fewer than 33% of samples.
hmp2_pathabun_ileum_CD_filt <- hmp2_pathabun_ileum_CD[, -which(colSums(hmp2_pathabun_ileum_CD > 0) < 9)]

host_rna_rpm_ileum_CD <- host_rna_rpm_ileum[CD_ileum_samples, ]

##### Haberman et al. 2014 --- ileal biopsies at time of diagnosis
### Top ileum CD genes overexpressed vs controls: DUOXA2, MMP3, AQP9, IL8, and DUOX2
### Top ileum CD gene underexpressed vs controls: APOA1, NAT8, AGXT2, CUBN, FAM151A
### Esp. co-expression of DUOX2 and APOA1
### Not described in this paper, but NOD2 is also a major gene of interest since the key genetic variants underlying CD risk are in this gene.

CD_ileal_genes_of_interest <- c("DUOXA2", "MMP3", "AQP9", "IL8", "DUOX2", "APOA1", "NAT8", "AGXT2", "CUBN", "FAM151A", "NOD2")

# Cluster these genes based on complement of spearman correlations to avoid redundant tests.
host_rna_rpm_ileum_CD_goi_spearman_dist <- as.dist(1 - cor(host_rna_rpm_ileum_CD[, CD_ileal_genes_of_interest], method="spearman"))

host_rna_rpm_ileum_CD_goi_hclust <- hcut(host_rna_rpm_ileum_CD_goi_spearman_dist, k = 6)

saveRDS(object = host_rna_rpm_ileum_CD_goi_hclust,
        file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/results_out/host_rna_rpm_ileum_CD_goi_hclust.rds")

# Get representative gene from each cluster.
sort(host_rna_rpm_ileum_CD_goi_hclust$cluster)
cluster_genes_of_interest <- c("DUOX2",  "MMP3", "AQP9", "APOA1", "NAT8", "NOD2")


saveRDS(object=list(pathabun=hmp2_pathabun_ileum_CD_filt, rnaseq=host_rna_rpm_ileum_CD), file = "results_out/prepped_pathabun_rnaseq_tables.rds")

# JUST NEEDED TO RUN THIS ONCE TO GET SIGNIFICANT ASSOCIATIONS.
# WROTE THIS TO RDS FILE SO NO NEED TO RE-RUN UNLESS SOMETHING IS CHANGED.

hmp2_pathabun_vs_rnaseq_CD_ileum <- data.frame(matrix(NA, nrow=length(cluster_genes_of_interest) * ncol(hmp2_pathabun_ileum_CD_filt), ncol=3))
colnames(hmp2_pathabun_vs_rnaseq_CD_ileum) <- c("pathway", "gene", "p")

i = 0
for(pathway in colnames(hmp2_pathabun_ileum_CD_filt)) {
  for(gene in cluster_genes_of_interest) {

    tmp_input <- data.frame(pathway=hmp2_pathabun_ileum_CD_filt[, pathway],
                            gene=host_rna_rpm_ileum_CD[, gene],
                            site=diagnosis_map_CD$site_name,
                            consent_age=diagnosis_map_CD$consent_age)

    hmp2_pathabun_vs_rnaseq_CD_ileum[i, c("pathway", "gene")] <- c(pathway, gene)

    model_out1 <- blmer(formula = gene ~ pathway + (1 | site) + (1 | consent_age), data=tmp_input, REML=FALSE)
    model_out2 <- blmer(formula = gene ~ (1 | site) + (1 | consent_age), data=tmp_input, REML=FALSE)

    anova_out <- anova(model_out1, model_out2)

    hmp2_pathabun_vs_rnaseq_CD_ileum[i, "p"] <- anova_out$`Pr(>Chisq)`[[2]]

    i = i + 1

  }
}

hmp2_pathabun_vs_rnaseq_CD_ileum$fdr <- p.adjust(hmp2_pathabun_vs_rnaseq_CD_ileum$p, "fdr")
hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1 <- hmp2_pathabun_vs_rnaseq_CD_ileum[which(hmp2_pathabun_vs_rnaseq_CD_ileum$fdr < 0.1), ]

saveRDS(object = hmp2_pathabun_vs_rnaseq_CD_ileum, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/results_out/hmp2_pathabun_vs_rnaseq_CD_ileum.rds")
saveRDS(object = hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/results_out/hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.rds")
