
rm(list=ls(all=TRUE))

setwd("/home/gavin/projects/picrust_pipeline/data/validation/infant_skin_ITS/")

# Filter OTU table using same criteria as paper.
infant_otus <- read.table("fungal_infant1/data/ninja_otutable_newCut_meta2.biom.tsv",
                          skip=1, header=T, sep="\t", row.names=1, check.names=FALSE, comment.char="")
infant_otus_filt <- infant_otus[, -which(colnames(infant_otus) %in% c("Positive", "Negative"))]
infant_otus_filt <- infant_otus_filt[, colSums(infant_otus_filt) > 50]
infant_otus_filt <- infant_otus_filt[rowSums(infant_otus_filt > 0) > 1, ]

infant_otus_filt_relab <- data.frame(sweep(infant_otus_filt, 2, colSums(infant_otus_filt), '/'), check.names = FALSE)

infant_map <- read.table("fungal_infant1/data/map.txt",
                         header=TRUE, sep="\t", comment.char="", stringsAsFactors = FALSE, row.names=1)

infant_map <- infant_map[colnames(infant_otus_filt),]

infant_map_baby <- infant_map[which(infant_map$motherorbaby_M_B == "B"),]
infant_map_mother <- infant_map[which(infant_map$motherorbaby_M_B == "M"),]

infant_skin_samples <- rownames(infant_map_baby[which(infant_map_baby$bodysite_OralMucosa_Forehead_VolarRight_PalmRight_FootRight_Vaginal_Anal_Feces == "Forehead"),])
infant_oral_samples <- rownames(infant_map_baby[which(infant_map_baby$bodysite_OralMucosa_Forehead_VolarRight_PalmRight_FootRight_Vaginal_Anal_Feces == "Oralmucosa"),])
infant_anal_samples <- rownames(infant_map_baby[which(infant_map_baby$bodysite_OralMucosa_Forehead_VolarRight_PalmRight_FootRight_Vaginal_Anal_Feces == "Anal"),])

infant_map_baby_skin <- infant_map_baby[infant_skin_samples,]
infant_map_baby_oral <- infant_map_baby[infant_oral_samples,]

infant_oral_vaginal_samples <- rownames(infant_map_baby_oral[which(infant_map_baby_oral$Delivery_Vvaginal_Ccs_IcsInoc == "V"),])
infant_oral_csection_samples <- rownames(infant_map_baby_oral[which(infant_map_baby_oral$Delivery_Vvaginal_Ccs_IcsInoc == "C"),])



# Read in stratified pathway abundances, aggregate on OTUs left after filtering OTU table, and convert to relative abundances.
infant_ec_strat <- read.table("picrust2_full_output/ec_ITS_counts_metagenome_out/pred_metagenome_strat.tsv",
                                    header=TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

infant_ec_strat_filt <- infant_ec_strat[which(infant_ec_strat$sequence %in% rownames(infant_otus_filt)), c("sequence", "function", colnames(infant_otus_filt))]

infant_ec_strat_filt_noseq <- infant_ec_strat_filt[, -which(colnames(infant_ec_strat_filt) == "sequence")]

colnames(infant_ec_strat_filt_noseq)[1] <- "func"

infant_ec_strat_summed <- aggregate(. ~ func, data=infant_ec_strat_filt_noseq, FUN=sum)

rownames(infant_ec_strat_summed) <- infant_ec_strat_summed$func

infant_ec_strat_summed <- infant_ec_strat_summed[, -1]

infant_ec_strat_summed_relab <- data.frame(sweep(infant_ec_strat_summed, 2, colSums(infant_ec_strat_summed), '/'), check.names = FALSE) * 100



ecs2test <- infant_ec_strat_filt[which(infant_ec_strat_filt$sequence == "SH180864.07FU_AJ698048_refs"), "function"]

boxplot(as.numeric(infant_otus_filt_relab["SH180864.07FU_AJ698048_refs", infant_oral_csection_samples]),
        as.numeric(infant_otus_filt_relab["SH180864.07FU_AJ698048_refs", infant_oral_vaginal_samples]))
p_val <- c()

for(ec in ecs2test) {
  wilcox_test_out <- wilcox.test(as.numeric(infant_ec_strat_summed_relab[ec, infant_oral_csection_samples]),
                                 as.numeric(infant_ec_strat_summed_relab[ec, infant_oral_vaginal_samples]))  
  p_val <- c(p_val, wilcox_test_out$p.value)
}







# Read in stratified pathway abundances, aggregate on OTUs left after filtering OTU table, and convert to relative abundances.
infant_pathabun_strat <- read.table("picrust2_full_output/pathways_out/path_abun_strat.tsv",
                               header=TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

infant_pathabun_strat_filt <- infant_pathabun_strat[which(infant_pathabun_strat$sequence %in% rownames(infant_otus_filt)), c("sequence", "pathway", colnames(infant_otus_filt))]


infant_pathabun_strat_filt[which(infant_pathabun_strat_filt$sequence == "SH180864.07FU_AJ698048_refs"), "pathway"]


infant_pathabun_strat_filt_noseq <- infant_pathabun_strat_filt[, -which(colnames(infant_pathabun_strat_filt) == "sequence")]

infant_pathabun_strat_summed <- aggregate(. ~ pathway, data=infant_pathabun_strat_filt_noseq, FUN=sum)

rownames(infant_pathabun_strat_summed) <- infant_pathabun_strat_summed$pathway

infant_pathabun_strat_summed <- infant_pathabun_strat_summed[, -1]

infant_pathabun_strat_summed_relab <- data.frame(sweep(infant_pathabun_strat_summed, 2, colSums(infant_pathabun_strat_summed), '/'), check.names = FALSE) * 100

p_val <- c()

paths2test <- infant_pathabun_strat_filt[which(infant_pathabun_strat_filt$sequence == "SH180864.07FU_AJ698048_refs"), "pathway"]

for(pathway in paths2test) {
  wilcox_test_out <- wilcox.test(as.numeric(infant_pathabun_strat_summed_relab[pathway, infant_oral_csection_samples]),
                     as.numeric(infant_pathabun_strat_summed_relab[pathway, infant_oral_vaginal_samples]))  
  p_val <- c(p_val, wilcox_test_out$p.value)
}

paths2test[which(p_val < 0.05)]
# [1] "PWY-6609" "PWY-7118" "PWY-7385"
# Proceeding just with the first one since PWY-7385 is just an engineered pathway.


# Write out table of stratified pathway abundances and write out sample groupings as well.
infant_pathabun_strat_filt[, 3:ncol(infant_pathabun_strat_filt)] <- data.frame(sweep(infant_pathabun_strat_filt[, 3:ncol(infant_pathabun_strat_filt)],
                                                                                     2,
                                                                                     colSums(infant_pathabun_strat_filt[, 3:ncol(infant_pathabun_strat_filt)]),
                                                                                     '/'),
                                                                               check.names=FALSE) * 100

infant_pathabun_strat_filt_sig_path <- infant_pathabun_strat_filt[which(infant_pathabun_strat_filt$pathway %in% c("PWY-6609", "PWY-7118")), ]



write.table(x=infant_pathabun_strat_filt_sig_path,
            file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/infant_mycobiome_sig_path_strat_abun.tsv",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE,
            sep="\t")

oral_birth_sample_groups <- data.frame(sample=c(infant_oral_vaginal_samples, infant_oral_csection_samples),
                                       birth=c(rep_len("V", length(infant_oral_vaginal_samples)),
                                               rep_len("C", length(infant_oral_csection_samples))))

write.table(x=oral_birth_sample_groups,
            file="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/infant_mycobiome_oral_birthtype_samples.tsv",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE,
            sep="\t")






infant_oral_vaginal_ASVs <- names(which(rowSums(infant_otus_filt[, infant_oral_vaginal_samples]) > 0))

infant_ec_predictions <- read.table("picrust2_full_output/ec_ITS_counts_predicted.tsv",
                                    header=TRUE, sep="\t", row.names=1)
infant_ec_predictions_subset_vaginal_oral <- infant_ec_predictions[infant_oral_vaginal_ASVs,]

p_values <- c()
for(ec in colnames(infant_ec_predictions_subset_vaginal_oral)) {
  ec_p <- length(which(infant_ec_predictions_subset_vaginal_oral[, ec] >= infant_ec_predictions["SH180864.07FU_AJ698048_refs", ec]))/nrow(infant_ec_predictions_subset_vaginal_oral)
  p_values <- c(p_values, ec_p)
}

colnames(infant_ec_predictions_subset_vaginal_oral)[which(p_values < 0.05)] 55656555

