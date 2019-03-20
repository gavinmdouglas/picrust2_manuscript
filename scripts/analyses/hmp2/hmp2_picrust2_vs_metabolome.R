### Test for associations between PICRUSt2 predictions and metabolite levels.
### For certain genes of interest (i.e. those previously known to be biomarkers of CD and/or UC).

rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/analyses/hmp2/hmp2_util_functions.R")

# Read in HMP2 pathway abundances and transform with arcsine square root transformation.
hmp2_pathabun <- readRDS("prepped_tables/hmp2_pathabun_filt_asin_Ileum.rds")

# Load RNA-seq data.
hmp2_metabolite <- readRDS(file = "prepped_tables/hmp2_metabolome_filt_subset_log.rds")

# Read in hmp2 metadata table.
hmp2_metadata <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_2018-08-20_col_subset.csv",
                            header=TRUE, sep=",", comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

metabolite_meta <- hmp2_metadata[which(hmp2_metadata$data_type == "metabolomics"), ]
metabolite_meta <- metabolite_meta[-which(duplicated(metabolite_meta$`Participant ID`)), ]
rownames(metabolite_meta) <- metabolite_meta$`Participant ID`

diagnosis_map <- data.frame(sample=metabolite_meta$`Participant ID`,
                            diagnosis=metabolite_meta$diagnosis,
                            consent_age=metabolite_meta$consent_age,
                            site_name=metabolite_meta$site_name,
                            stringsAsFactors = FALSE)

overlapping_samples <- rownames(hmp2_pathabun)[which(rownames(hmp2_pathabun) %in% rownames(hmp2_metabolite))]

diagnosis_map <- diagnosis_map[which(diagnosis_map$sample %in% overlapping_samples), ]

CD_ileum_samples <- diagnosis_map[which(diagnosis_map$diagnosis == "CD"), "sample"]

diagnosis_map_CD <- diagnosis_map[which(diagnosis_map$sample %in% CD_ileum_samples), ]

diagnosis_map_CD <- diagnosis_map_CD[which(diagnosis_map_CD$sample %in% overlapping_samples), ]

diagnosis_map_CD$site_name <- factor(diagnosis_map_CD$site_name)

hmp2_pathabun_ileum_CD <- hmp2_pathabun[CD_ileum_samples, ]
hmp2_metabolite_ileum_CD <- hmp2_metabolite[CD_ileum_samples, ]

# Remove features found in fewer than 33% of samples.
hmp2_pathabun_ileum_CD_filt <- hmp2_pathabun_ileum_CD[, -which(colSums(hmp2_pathabun_ileum_CD > 0) < 9)]
hmp2_metabolite_ileum_CD_filt <- hmp2_metabolite_ileum_CD[, -which(colSums(hmp2_metabolite_ileum_CD > 0) < 9)]

saveRDS(object=list(pathabun=hmp2_pathabun_ileum_CD_filt, metabolite=hmp2_metabolite_ileum_CD_filt), file = "results_out/prepped_pathabun_metabolite_tables.rds")

hmp2_pathabun_vs_metabolite_CD_ileum <- data.frame(matrix(NA, nrow=ncol(hmp2_metabolite_ileum_CD_filt) * ncol(hmp2_pathabun_ileum_CD_filt), ncol=4))
colnames(hmp2_pathabun_vs_metabolite_CD_ileum) <- c("pathway", "metabolite", "R", "p")

i = 1
for(pathway in colnames(hmp2_pathabun_ileum_CD_filt)) {
  for(metabolite in colnames(hmp2_metabolite_ileum_CD_filt)) {
    
    tmp_input <- data.frame(pathway=hmp2_pathabun_ileum_CD_filt[, pathway],
                            metabolite=hmp2_metabolite_ileum_CD_filt[, metabolite],
                            consent_age=diagnosis_map_CD$consent_age)
    
    # BLMER commands (unused):
    # model_out1 <- blmer(formula = gene ~ pathway + (1 | site), data=tmp_input, REML=FALSE)
    # model_out2 <- blmer(formula = gene ~ (1 | site), data=tmp_input, REML=FALSE)
    # anova_out <- anova(model_out1, model_out2)
    # hmp2_pathabun_vs_metabolite_CD_ileum[i, "p"] <- anova_out$`Pr(>Chisq)`[[2]]
    # hmp2_pathabun_vs_metabolite_CD_ileum[i, c("intercept", "pathway_coef")] <- summary(model_out1)$coefficients[, "Estimate"]
    
    partial_spearman <- pcor.test(tmp_input$pathway, tmp_input$metabolite, tmp_input$consent_age, method = "spearman")
    hmp2_pathabun_vs_metabolite_CD_ileum[i, c("pathway", "metabolite")] <- c(pathway, metabolite)
    hmp2_pathabun_vs_metabolite_CD_ileum[i, c("R", "p")] <- c(partial_spearman$estimate, partial_spearman$p.value)
    
    i = i + 1
    
  }
}

hmp2_pathabun_vs_metabolite_CD_ileum$fdr <- p.adjust(hmp2_pathabun_vs_metabolite_CD_ileum$p, "fdr")

saveRDS(object = hmp2_pathabun_vs_metabolite_CD_ileum, file =
          "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/hmp2_tables/results_out/hmp2_pathabun_vs_metabolite_CD_ileum.rds")
