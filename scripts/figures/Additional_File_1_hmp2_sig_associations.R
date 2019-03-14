rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/results_out")

rnaseq_sig <- readRDS("hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.rds")

write.table(x = rnaseq_sig,
            file = "hmp2_pathabun_vs_rnaseq_CD_ileum_fdr0.1.tsv",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE,
            sep="\t")


metabolite_sig <- readRDS("hmp2_pathabun_vs_metabolite_CD_ileum_fdr0.1.rds")

write.table(x = metabolite_sig,
            file = "hmp2_pathabun_vs_metabolite_CD_ileum_fdr0.1.tsv",
            quote = FALSE,
            row.names=FALSE,
            col.names=TRUE,
            sep="\t")

### Downloaded these tables and put into single excel file manually.

