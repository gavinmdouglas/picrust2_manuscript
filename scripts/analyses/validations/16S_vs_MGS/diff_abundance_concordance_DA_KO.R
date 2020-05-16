rm(list=ls(all.names=TRUE))

library("ALDEx2")
library("DESeq2")

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/")
source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

# Read in metadata and KO tables limited to overlapping KOs.
datasets <- c("hmp", "indian", "cameroon", "primate")

ko_metrics_out <- list()

for(dataset in datasets) {
  ko_metrics_out[[dataset]] <- list()
  ko_metrics_out[[dataset]][["infiles"]] <- read_in_ko_predictions(dataset)
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

######### Add in alternative mgs KOs as well.
cameroon_alt_mgs <- read.table("alt_mgs_output/cameroon_genefamilies_unstratified.tsv",
                               header=TRUE, sep="\t", row.names=1, comment.char="")
colnames(cameroon_alt_mgs) <- gsub("_Abundance.RPKs", "", colnames(cameroon_alt_mgs))
cameroon_metadata <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/cameroon/PRJEB27005.txt",
                                header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")
rownames(cameroon_metadata) <- cameroon_metadata$run_accession
colnames(cameroon_alt_mgs) <- cameroon_metadata[colnames(cameroon_alt_mgs), "sample_alias"]
cameroon_alt_mgs <- cameroon_alt_mgs[which(rownames(cameroon_alt_mgs) %in% rownames(ko_metrics_out$cameroon$infiles$all_kos_overlap$mgs_ko)), ]
cameroon_alt_mgs <- add_missing_funcs(in_df = cameroon_alt_mgs, all_funcs = rownames(ko_metrics_out$cameroon$infiles$all_kos_overlap$mgs_ko))
cameroon_alt_mgs <- cameroon_alt_mgs[rownames(ko_metrics_out$cameroon$infiles$all_kos_overlap$mgs_ko), ]

hmp_alt_mgs <- read.table("alt_mgs_output/hmp_genefamilies_unstratified.tsv",
                          header=TRUE, sep="\t", row.names=1, comment.char="")
colnames(hmp_alt_mgs) <- gsub("_R1_R2_cat_Abundance.RPKs", "", colnames(hmp_alt_mgs))
hmp_alt_mgs <- hmp_alt_mgs[which(rownames(hmp_alt_mgs) %in% rownames(ko_metrics_out$hmp$infiles$all_kos_overlap$mgs_ko)), ]
hmp_alt_mgs <- add_missing_funcs(in_df = hmp_alt_mgs, all_funcs = rownames(ko_metrics_out$hmp$infiles$all_kos_overlap$mgs_ko))
hmp_alt_mgs <- hmp_alt_mgs[rownames(ko_metrics_out$hmp$infiles$all_kos_overlap$mgs_ko), ]


indian_alt_mgs <- read.table("alt_mgs_output/indian_genefamilies_unstratified.tsv",
                             header=TRUE, sep="\t", row.names=1, comment.char="")
colnames(indian_alt_mgs) <- gsub("_Abundance.RPKs", "", colnames(indian_alt_mgs))

indian_metadata <- read.table("alt_mgs_output/PRJNA397112.txt",
                              header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")
rownames(indian_metadata) <- indian_metadata$run_accession
colnames(indian_alt_mgs) <- indian_metadata[colnames(indian_alt_mgs), "sample_alias"]
indian_alt_mgs <- indian_alt_mgs[which(rownames(indian_alt_mgs) %in% rownames(ko_metrics_out$indian$infiles$all_kos_overlap$mgs_ko)), ]
indian_alt_mgs <- add_missing_funcs(in_df = indian_alt_mgs, all_funcs = rownames(ko_metrics_out$indian$infiles$all_kos_overlap$mgs_ko))
indian_alt_mgs <- indian_alt_mgs[rownames(ko_metrics_out$indian$infiles$all_kos_overlap$mgs_ko), ]


primate_alt_mgs <- read.table("alt_mgs_output/primate_genefamilies_unstratified.tsv",
                              header=TRUE, sep="\t", row.names=1, comment.char="")
colnames(primate_alt_mgs) <- gsub("_kneaddata_Abundance.RPKs", "", colnames(primate_alt_mgs))
primate_alt_mgs <- primate_alt_mgs[which(rownames(primate_alt_mgs) %in% rownames(ko_metrics_out$primate$infiles$all_kos_overlap$mgs_ko)), ]
primate_alt_mgs <- add_missing_funcs(in_df = primate_alt_mgs, all_funcs = rownames(ko_metrics_out$primate$infiles$all_kos_overlap$mgs_ko))
primate_alt_mgs <- primate_alt_mgs[rownames(ko_metrics_out$primate$infiles$all_kos_overlap$mgs_ko), ]


ko_metrics_out$cameroon$infiles$all_kos_overlap$mgs_ko_alt <- cameroon_alt_mgs
ko_metrics_out$hmp$infiles$all_kos_overlap$mgs_ko_alt <- hmp_alt_mgs
ko_metrics_out$indian$infiles$all_kos_overlap$mgs_ko_alt <- indian_alt_mgs
ko_metrics_out$primate$infiles$all_kos_overlap$mgs_ko_alt <- primate_alt_mgs
###############

musicc_uscgs <- read.table("../MUSiCC_KEGG_single_copy_genes.txt", stringsAsFactors = FALSE)$V1

wilcoxon_musicc_out <- list()

for(d in datasets) {
  wilcoxon_musicc_out[[d]] <- lapply(ko_metrics_out[[d]]$infiles$all_kos_overlap,
                                  function(x) {
                                    run_wilcoxon_musicc_tests(table = x,
                                                              metadata=validation_groups[[d]],
                                                              uscg_set = musicc_uscgs,
                                                              dataset_name = d)
                                    })
}
names(wilcoxon_musicc_out) <- datasets

wilcoxon_musicc_out_perf_0.05 <- lapply(datasets, function(d) {
                            summarize_test_performance(table_list = wilcoxon_musicc_out[[d]],
                                                       sig_column="wilcox_BH",
                                                       reference_name = "mgs_ko",
                                                       sig_cutoff=0.05) })
names(wilcoxon_musicc_out_perf_0.05) <- datasets




# Test for significant differences with ALDEx2.
# First get subset of datasets to keep for this analysis (corresponding to counts, which are required for ALDEx2 and DESeq2.
for(d in datasets) {
  
  ko_metrics_out[[d]]$infiles$all_kos_overlap_subset <- list()
  
  for(d2keep in c("picrust2_ko_nsti2", "picrust1_ko", "piphillin_ko", "mgs_ko", "mgs_ko_alt", "picrust2_ko_scrambled_mean")) {
    ko_metrics_out[[d]]$infiles$all_kos_overlap_subset[[d2keep]] <- ko_metrics_out[[d]]$infiles$all_kos_overlap[[d2keep]]
  }
}

# Run ALDEx2
aldex2_out <- list()

for(d in datasets) {
  aldex2_out[[d]] <- lapply(ko_metrics_out[[d]]$infiles$all_kos_overlap_subset,
                              function(x) {
                                run_aldex2_two_groups(table = x, metadata=validation_groups[[d]], dataset_name = d) })
}
names(aldex2_out) <- datasets

aldex2_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = aldex2_out[[d]],
                             sig_column="wi.eBH",
                             reference_name = "mgs_ko",
                             sig_cutoff=0.05) })
names(aldex2_out_perf_0.05) <- datasets


# Run DESeq2
deseq2_out_res0.05 <- list()


for(d in datasets) {

  deseq2_out_res0.05[[d]] <- lapply(ko_metrics_out[[d]]$infiles$all_kos_overlap_subset,
                                    function(x) {
                                      run_default_deseq2_two_groups(table = x, metadata=validation_groups[[d]], alpha.set=0.05, dataset_name = d) })
}

names(deseq2_out_res0.05) <- datasets

deseq2_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = deseq2_out_res0.05[[d]],
                             sig_column="padj",
                             reference_name = "mgs_ko",
                             sig_cutoff=0.05) })
names(deseq2_out_perf_0.05) <- datasets


# Run DESeq2 with different options specifically recommended for microbiome data
deseq2_GMmeans_out_res0.05 <- list()

for(d in datasets) {
  deseq2_GMmeans_out_res0.05[[d]] <- lapply(ko_metrics_out[[d]]$infiles$all_kos_overlap_subset,
                                            function(x) {
                                              run_GMmean_deseq2_two_groups(table = x, metadata=validation_groups[[d]], alpha.set=0.05, dataset_name = d) })
  
}


deseq2_GMmeans_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = deseq2_GMmeans_out_res0.05[[d]],
                             sig_column="padj",
                             reference_name = "mgs_ko",
                             sig_cutoff=0.05) })
names(deseq2_GMmeans_out_perf_0.05) <- datasets


# Run tests for differential prevalence (based on function presence / absence)
fishers_out <- list()

for(d in datasets) {
  fishers_out[[d]] <- mclapply(ko_metrics_out[[d]]$infiles$all_kos_overlap,
                               function(x) {
                                 differential_prevalence(table = x,
                                                         metadata=validation_groups[[d]],
                                                         dataset_name = d) }, mc.cores=10)
}
names(fishers_out) <- datasets

fishers_out_perf_0.05 <- lapply(datasets, function(d) {
  summarize_test_performance(table_list = fishers_out[[d]],
                             sig_column="fishers_BH",
                             reference_name = "mgs_ko",
                             sig_cutoff=0.05) })
names(fishers_out_perf_0.05) <- datasets

# Save RDS of DA tests.
ko_da_out <- list(validation_groups=validation_groups,
                  wilcoxon_musicc_out=wilcoxon_musicc_out,
                  wilcoxon_musicc_out_perf_0.05=wilcoxon_musicc_out_perf_0.05,
                  
                  aldex2_out=aldex2_out,
                  aldex2_out_perf_0.05=aldex2_out_perf_0.05,
                  
                  deseq2_out_res0.05=deseq2_out_res0.05,
                  deseq2_out_perf_0.05=deseq2_out_perf_0.05,
                  
                  deseq2_GMmeans_out_res0.05=deseq2_GMmeans_out_res0.05,
                  deseq2_GMmeans_out_perf_0.05=deseq2_GMmeans_out_perf_0.05,
                  
                  fishers_out=fishers_out,
                  fishers_out_perf_0.05=fishers_out_perf_0.05)

saveRDS(object = ko_da_out, file = "../../data/saved_RDS/DA_concordance/ko_da_out.rds")

### Sanity checks on individual steps of this pipeline

# Wilcoxon / MUSICC test
# tmp <- ko_metrics_out$cameroon$infiles$all_kos_overlap$picrust2_ko_nsti2
# 
# tmp_medians <- colMedians(as.matrix(tmp[musicc_uscgs, ]))
# 
# tmp_normalized <-  data.frame(sweep(tmp, 2, colMedians(as.matrix(tmp[musicc_uscgs, ])), `/`))
# tmp_normalized <- tmp_normalized[-which(rownames(tmp_normalized) %in% musicc_uscgs), ]
# 
# tmp2 <- c()
# 
# test_set1 <- validation_groups$cameroon$group1[which(validation_groups$cameroon$group1 %in% colnames(tmp_normalized))]
# test_set2 <- validation_groups$cameroon$group2[which(validation_groups$cameroon$group2 %in% colnames(tmp_normalized))]
# 
# for(ko in rownames(tmp_normalized)) {
#   tmp2 <- c(tmp2,
#             wilcox.test(as.numeric(tmp_normalized[ko, test_set1]),
#                         as.numeric(tmp_normalized[ko, test_set2]))$p.value)
# }
# 
# tmp2_fdr <- p.adjust(tmp2, "fdr")
# 
# identical(tmp2_fdr, wilcoxon_musicc_out$cameroon$picrust2_ko_nsti2$wilcox_BH)



# ALDEx2
# tmp <- round(ko_metrics_out$hmp$infiles$all_kos_overlap$picrust2_ko_nsti2)
# test_set1 <- validation_groups$hmp$group1[which(validation_groups$hmp$group1 %in% colnames(tmp))]
# test_set2 <- validation_groups$hmp$group2[which(validation_groups$hmp$group2 %in% colnames(tmp))]
# 
# 
# aldex2_tmp  <- aldex(reads = tmp[, c(test_set1, test_set2)],
#                       conditions=c(rep("group1", length(test_set1)),
#                                    rep("group2", length(test_set2))),
#                       mc.samples = 128, test = "t", effect = TRUE,
#                       include.sample.summary = FALSE,
#                       verbose = FALSE,
#                       denom = "all",
#                       iterate = FALSE)
# 
# expected_aldex2_out <- aldex2_out$hmp$picrust2_ko_nsti2[rownames(aldex2_tmp), ]
# 
# 
# #weird_outliers_to_remove <- which(is.na(expected_aldex2_out$rab.all))
# #length(weird_outliers_to_remove)
# 
# #expected_aldex2_out <- expected_aldex2_out[-weird_outliers_to_remove, ]
# #aldex2_tmp <- aldex2_tmp[-weird_outliers_to_remove, ]
# 
# plot(expected_aldex2_out$wi.eBH, aldex2_tmp$wi.eBH)

# deseq2 - default and GMeans
# tmp <- floor(ko_metrics_out$indian$infiles$all_kos_overlap$picrust2_ko_nsti2)
# test_set1 <- validation_groups$indian$group1[which(validation_groups$indian$group1 %in% colnames(tmp))]
# test_set2 <- validation_groups$indian$group2[which(validation_groups$indian$group2 %in% colnames(tmp))]
# 
# metadata_tmp <- data.frame(sample=c(test_set1, test_set2),
#                            group=c(rep("group1", length(test_set1)),
#                                    rep("group2", length(test_set2))))
# 
# dds <- DESeqDataSetFromMatrix(countData = tmp[, c(test_set1, test_set2)],
#                               colData = metadata_tmp,
#                               design = ~ group)
# default_deseq2 <- DESeq(dds)
# default_deseq2_results <- results(default_deseq2, alpha=0.05)
# 
# all.equal(default_deseq2_results, deseq2_out_res0.05$indian$picrust2_ko_nsti2)
# 
# 
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# 
# geoMeans <- apply(counts(dds), 1, gm_mean)
# dds_estimated <- estimateSizeFactors(dds, geoMeans = geoMeans)
# dds_estimated = DESeq(dds_estimated, fitType="local" )
# 
# geoMeans_deseq2_results <- results(dds_estimated, alpha=0.05)
# 
# all.equal(geoMeans_deseq2_results, deseq2_GMmeans_out_res0.05$indian$picrust2_ko_nsti2)

# # Fisher's exact tests based on presence / absence
# tmp <- ko_metrics_out$indian$infiles$all_kos_overlap$picrust2_ko_nsti0.5
# 
# tmp2 <- c()
# 
# test_set1 <- validation_groups$indian$group1[which(validation_groups$indian$group1 %in% colnames(tmp))]
# test_set2 <- validation_groups$indian$group2[which(validation_groups$indian$group2 %in% colnames(tmp))]
# 
# for(ko in rownames(tmp)) {
#   tmp2 <- c(tmp2,
#             fisher.test(matrix(c(length(which(as.numeric(tmp[ko, test_set1]) > 0)),
#                                  length(which(as.numeric(tmp[ko, test_set1]) == 0)),
#                                  length(which(as.numeric(tmp[ko, test_set2]) > 0)),
#                                  length(which(as.numeric(tmp[ko, test_set2]) == 0))), nrow=2, ncol=2))$p.value)
# }
# 
# tmp2_fdr <- p.adjust(tmp2, "fdr")
# 
# identical(tmp2_fdr, fishers_out$indian$picrust2_ko_nsti0.5$fishers_BH)
