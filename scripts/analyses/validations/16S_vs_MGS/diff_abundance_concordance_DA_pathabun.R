rm(list=ls(all.names=TRUE))

library("ALDEx2")
library("DESeq2")

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

# Read in metadata and pathabun tables limited to overlapping pathabuns.
datasets <- c("hmp", "indian", "cameroon", "primate")

pathabun_metrics_out <- list()

for(dataset in datasets) {
  pathabun_metrics_out[[dataset]] <- list()
  pathabun_metrics_out[[dataset]][["infiles"]] <- read_in_pathway_predictions(dataset)
}

validation_groups <- readRDS("../../revisions_unsorted/validation_groupings.rds")

# Change sample names to match tables.
validation_groups$blueberry$group1 <- gsub("-", ".", validation_groups$blueberry$group1)
validation_groups$blueberry$group2 <- gsub("-", ".", validation_groups$blueberry$group2)

validation_groups$primate$group1 <- gsub("^", "X", validation_groups$primate$group1)
validation_groups$primate$group2 <- gsub("^", "X", validation_groups$primate$group2)

# Re-code HMP and cameroon validation categories to make them easier to work with.
validation_groups[["hmp"]] <- validation_groups$hmp_within_oral
validation_groups[["cameroon"]] <- validation_groups$cameroonian

# Comparison for hmp limited to 22 vs 36 samples
# Comparison for indian limited to 51 vs 38 samples
# Comparison for cameroon limited to 19 vs 36 samples
# Comparison for primate limited to 29 vs 29 samples


wilcoxon_out <- list()

for(d in datasets) {
  wilcoxon_out[[d]] <- lapply(pathabun_metrics_out[[d]]$infiles$all_pathabun,
                                  function(x) {
                                    run_wilcoxon_relab_tests(table = x,
                                                              metadata=validation_groups[[d]],
                                                              dataset_name = d)
                                    })
}
names(wilcoxon_out) <- datasets

wilcoxon_out_perf_0.05 <- lapply(datasets, function(d) {
                            summarize_test_performance(table_list = wilcoxon_out[[d]],
                                                       sig_column="wilcox_BH",
                                                       reference_name = "mgs_pathabun",
                                                       sig_cutoff=0.05) })
names(wilcoxon_out_perf_0.05) <- datasets




# Test for significant differences with ALDEx2.
# First get subset of datasets to keep for this analysis (because the next commands take a long time).
for(d in datasets) {
  
  pathabun_metrics_out[[d]]$infiles$all_pathabun_subset <- list()
  
  for(d2keep in c("picrust2_pathabun_nsti2", "mgs_pathabun", "picrust2_pathabun_scrambled_mean")) {
    pathabun_metrics_out[[d]]$infiles$all_pathabun_subset[[d2keep]] <- pathabun_metrics_out[[d]]$infiles$all_pathabun[[d2keep]]
  }
}

# Run ALDEx2
aldex2_out <- list()

for(d in datasets) {
  aldex2_out[[d]] <- lapply(pathabun_metrics_out[[d]]$infiles$all_pathabun_subset,
                              function(x) {
                                run_aldex2_two_groups(table = x, metadata=validation_groups[[d]], dataset_name = d) })
}
names(aldex2_out) <- datasets

aldex2_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = aldex2_out[[d]],
                             sig_column="wi.eBH",
                             reference_name = "mgs_pathabun",
                             sig_cutoff=0.05) })
names(aldex2_out_perf_0.05) <- datasets


# Run DESeq2
deseq2_out_res0.05 <- list()


for(d in datasets) {

  deseq2_out_res0.05[[d]] <- lapply(pathabun_metrics_out[[d]]$infiles$all_pathabun_subset,
                                    function(x) {
                                      run_default_deseq2_two_groups(table = x, metadata=validation_groups[[d]], alpha.set=0.05, dataset_name = d) })
}

names(deseq2_out_res0.05) <- datasets

deseq2_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = deseq2_out_res0.05[[d]],
                             sig_column="padj",
                             reference_name = "mgs_pathabun",
                             sig_cutoff=0.05) })
names(deseq2_out_perf_0.05) <- datasets


# Run DESeq2 with different options specifically recommended for microbiome data
deseq2_GMmeans_out_res0.05 <- list()

for(d in datasets) {
  deseq2_GMmeans_out_res0.05[[d]] <- lapply(pathabun_metrics_out[[d]]$infiles$all_pathabun_subset,
                                            function(x) {
                                              run_GMmean_deseq2_two_groups(table = x, metadata=validation_groups[[d]], alpha.set=0.05, dataset_name = d) })
  
}


deseq2_GMmeans_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = deseq2_GMmeans_out_res0.05[[d]],
                             sig_column="padj",
                             reference_name = "mgs_pathabun",
                             sig_cutoff=0.05) })
names(deseq2_GMmeans_out_perf_0.05) <- datasets


# Run tests for differential prevalence (based on function presence / absence)
fishers_out <- list()

for(d in datasets) {
  fishers_out[[d]] <- mclapply(pathabun_metrics_out[[d]]$infiles$all_pathabun,
                               function(x) {
                                 differential_prevalence(table = x,
                                                         metadata=validation_groups[[d]],
                                                         dataset_name = d) }, mc.cores=10)
}
names(fishers_out) <- datasets

fishers_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = fishers_out[[d]],
                             sig_column="fishers_BH",
                             reference_name = "mgs_pathabun",
                             sig_cutoff=0.05) })
names(fishers_out_perf_0.05) <- datasets

# Save RDS of DA tests.
pathabun_da_out <- list(validation_groups=validation_groups,
                  wilcoxon_out=wilcoxon_out,
                  wilcoxon_out_perf_0.05=wilcoxon_out_perf_0.05,
                  
                  aldex2_out=aldex2_out,
                  aldex2_out_perf_0.05=aldex2_out_perf_0.05,
                  
                  deseq2_out_res0.05=deseq2_out_res0.05,
                  deseq2_out_perf_0.05=deseq2_out_perf_0.05,
                  
                  deseq2_GMmeans_out_res0.05=deseq2_GMmeans_out_res0.05,
                  deseq2_GMmeans_out_perf_0.05=deseq2_GMmeans_out_perf_0.05,
                  
                  fishers_out=fishers_out,
                  fishers_out_perf_0.05=fishers_out_perf_0.05)

saveRDS(object = pathabun_da_out, file = "../../data/saved_RDS/DA_concordance/pathabun_da_out.rds")

### Sanity checks on wilcoxon test on relative abundance.

# tmp <- pathabun_metrics_out$hmp$infiles$all_pathabun$picrust2_pathabun_nsti2
# 
# tmp_normalized <-  data.frame(sweep(tmp, 2, colSums(tmp), `/`)) * 100
# 
# tmp2 <- c()
# 
# test_set1 <- validation_groups$hmp$group1[which(validation_groups$hmp$group1 %in% colnames(tmp_normalized))]
# test_set2 <- validation_groups$hmp$group2[which(validation_groups$hmp$group2 %in% colnames(tmp_normalized))]
# 
# for(pathabun in rownames(tmp_normalized)) {
#   tmp2 <- c(tmp2,
#             wilcox.test(as.numeric(tmp_normalized[pathabun, test_set1]),
#                         as.numeric(tmp_normalized[pathabun, test_set2]))$p.value)
# }
# 
# tmp2_fdr <- p.adjust(tmp2, "fdr")
# 
# identical(tmp2_fdr, wilcoxon_out$hmp$picrust2_pathabun_nsti2$wilcox_BH)
