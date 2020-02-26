rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/results_out")

rnaseq_sig <- readRDS("hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.rds")

# Add column containing list of other genes in cluster.
host_rna_rpm_ileum_CD_goi_hclust <- readRDS("host_rna_rpm_ileum_CD_goi_hclust.rds")
host_rna_rpm_ileum_CD_goi_hclust$cluster

rnaseq_sig$other_genes_in_cluster <- NA
rnaseq_sig[which(rnaseq_sig$gene == "DUOX2"), "other_genes_in_cluster"] <- "DUOXA2"
rnaseq_sig[which(rnaseq_sig$gene == "MMP3"), "other_genes_in_cluster"] <- "None"

rnaseq_sig <- rnaseq_sig[, c("pathway", "gene", "other_genes_in_cluster", "R", "p", "fdr")]

write.table(x = rnaseq_sig,
            file = "hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.tsv",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE,
            sep="\t")


### Downloaded this table and put into table manually.

