orig_level4ec <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/func_tables/euk_faa_uniref90_hits_level4ec_regrouped.txt",
                            header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", row.names=1, check.names=FALSE)

new_uniref <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/id_mappings/eukaryotic_microbes_cds_uniref90.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", row.names=1, check.names=FALSE)

new_level4ec <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/id_mappings/eukaryotic_microbes_cds_level4ec_regrouped.txt",
                           header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", row.names=1, check.names=FALSE)

rownames(new_level4ec) <- gsub("^", "EC:", rownames(new_level4ec))

sample_overlap <- colnames(orig_level4ec)[which(colnames(orig_level4ec) %in% colnames(new_level4ec))]
func_overlap <- rownames(orig_level4ec)[which(rownames(orig_level4ec) %in% rownames(new_level4ec))]

new_level4ec_subset <- new_level4ec[func_overlap, sample_overlap]
orig_level4ec_subset <- orig_level4ec[func_overlap, sample_overlap]

tmp_predictions <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline/ec_18S_predicted.tsv",
                              header=T, sep="\t", row.names=1)

tmp_predictions2 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_full_output_pipeline/EC_predicted.tsv",
                              header=T, sep="\t", row.names=1)

tmp_predictions3 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output_pipeline/marker_nsti_predicted.tsv",
                              header=T, sep="\t", row.names=1)

tmp_predictions4 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/picrust2_full_output_pipeline/marker_nsti_predicted.tsv",
                               header=T, sep="\t", row.names=1)