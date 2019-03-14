### Commands to clean-up HUMAnN2 output tables to make them easier to compare with PICRUSt2 tables.

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")

### Run these commands on all expected files given an input and output folder:
wrap_humann2_cleaner <- function(input_dir, output_dir, old_sample=NULL, new_sample=NULL) {
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_pathabundance_unstratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_pathabun_unstrat.tsv", sep="/"),
                        col_str_to_remove="_Abundance", first_col="pathway", strat=FALSE, rm_descrip=TRUE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_pathabundance_stratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_pathabun_strat.tsv", sep="/"),
                        col_str_to_remove="_Abundance", first_col="pathway", strat=TRUE, rm_descrip=TRUE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_pathcoverage_unstratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_pathcov_unstrat.tsv", sep="/"),
                        col_str_to_remove="_Coverage", first_col="pathway", strat=FALSE, rm_descrip=TRUE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_pathcoverage_stratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_pathcov_strat.tsv", sep="/"),
                        col_str_to_remove="_Coverage", first_col="pathway", strat=TRUE, rm_descrip=TRUE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_KO_unstratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_ko_unstrat.tsv", sep="/"),
                        col_str_to_remove="_Abundance.RPKs", first_col="function", strat=FALSE, rm_descrip=FALSE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_KO_stratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_ko_strat.tsv", sep="/"),
                        col_str_to_remove="_Abundance.RPKs", first_col="function", strat=TRUE, rm_descrip=FALSE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_level4ec_unstratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_ec_unstrat.tsv", sep="/"),
                        col_str_to_remove="_Abundance.RPKs", first_col="function", strat=FALSE, rm_descrip=FALSE, str2add="EC:",
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_level4ec_stratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_ec_strat.tsv", sep="/"),
                        col_str_to_remove="_Abundance.RPKs", first_col="function", strat=TRUE, rm_descrip=FALSE, str2add="EC:",
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_pfam_unstratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_pfam_unstrat.tsv", sep="/"),
                        col_str_to_remove="_Abundance.RPKs", first_col="function", strat=FALSE, rm_descrip=FALSE, replace_pf=TRUE,
                        old_sample=old_sample, new_sample=new_sample)
  
  clean_raw_humann2_out(filename = paste(input_dir, "humann2_pfam_stratified.tsv", sep="/"),
                        outfile = paste(output_dir, "humann2_pfam_strat.tsv", sep="/"),
                        col_str_to_remove="_Abundance.RPKs", first_col="function", strat=TRUE, rm_descrip=FALSE, replace_pf=TRUE,
                        old_sample=old_sample, new_sample=new_sample)

}



wrap_humann2_cleaner(input_dir = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/mgs/humann2_final_out",
                     output_dir = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/hmp")

wrap_humann2_cleaner(input_dir = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/mgs/humann2_final_out",
                     output_dir = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/mammalian_stool")

soil_ids_map <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S_to_MGS_ids.txt",
                           header=T, sep="\t", row.names=2, stringsAsFactors = FALSE)

wrap_humann2_cleaner(input_dir = "/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/mgs/humann2_final_out",
                     output_dir = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/soil_crossbiome",
                     old_sample = rownames(soil_ids_map),
                     new_sample = soil_ids_map$sample_id)

ocean_sample_map <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/ocean/ocean_16S_mgs_sample_links.txt",
                               header=T, sep="\t", stringsAsFactors = FALSE)

wrap_humann2_cleaner(input_dir = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/humann2_final_out",
                     output_dir = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/ocean",
                     old_sample = ocean_sample_map$runids_mgs,
                     new_sample = ocean_sample_map$runids_16S)


# Also clean-up blueberry humann2 files.

blueberry_soil_16S_names <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/16S_sample_names.txt",
                                       header=F, stringsAsFactors = FALSE)
blueberry_soil_mgs_names <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/mgs_sample_names.txt",
                                       header=F, stringsAsFactors = FALSE)

# Sanity check that they are in the same order.
#blueberry_sample_map <- data.frame(amplicon16S=blueberry_soil_16S_names$V1,
#                                   mgs=blueberry_soil_mgs_names$V1)
#blueberry_sample_map

wrap_humann2_cleaner(input_dir = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/humann2_final_out",
                     output_dir = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/mgs_validation/blueberry",
                     old_sample = blueberry_soil_mgs_names$V1,
                     new_sample = blueberry_soil_16S_names$V1)
