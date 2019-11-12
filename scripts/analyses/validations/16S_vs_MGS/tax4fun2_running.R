### Commands to run Tax4Fun2

rm(list=ls(all=TRUE))

library(Tax4Fun2)

# Function to format a TSV BIOM file for Tax4Fun2 and to run the pipeline.
# Output files are automatically written to the specified folder.

run_tax4fun2 <- function(in_fasta_path, in_tsv_biom_path, output_dir, num_threads=1, force=FALSE) {
  
  if(file.exists(output_dir)) {
    if(! force) {
      stop("Stopping - output directory already exists and force=FALSE")
    }
  } else {
    dir.create(path = output_dir)
  }
  
  # Format BIOM file for Tax4Fun2.
  in_biom <- read.table(file = in_tsv_biom_path, header=TRUE, skip=1, stringsAsFactors = FALSE, comment.char="", check.names = FALSE, sep="\t")
  colnames(in_biom)[1] <- "ID"
  in_formatted_biom_path <- paste(output_dir, "tax4fun2_input_seqabun.tsv", sep="/")
  write.table(x = in_biom, file = in_formatted_biom_path, col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  blast_out_dir <- paste(output_dir, "blast_out", sep="/")
  
  runRefBlast(path_to_otus = in_fasta_path,
              path_to_reference_data = "/home/gavin/local/Tax4Fun2_ref/Tax4Fun2_ReferenceData_v1.1",
              path_to_temp_folder = blast_out_dir,
              database_mode = "Ref100NR",
              use_force = FALSE,
              num_threads = num_threads)

  makeFunctionalPrediction(path_to_otu_table=in_formatted_biom_path,
                           path_to_reference_data="/home/gavin/local/Tax4Fun2_ref/Tax4Fun2_ReferenceData_v1.1",
                           path_to_temp_folder=blast_out_dir,
                           database_mode = 'Ref100NR',
                           normalize_by_copy_number = TRUE,
                           min_identity_to_reference = 97,
                           include_user_data = FALSE,
                           path_to_user_data = '',
                           use_uproc = TRUE,
                           normalize_pathways = FALSE)
}

# Run pipeline on all datasets.

# HMP
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/hmp_16S_rep_seqs.fasta",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/hmp_16S.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/tax4fun2_working",
             num_threads=40)

# Mammal
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/deblur_output_final/iGEM_16S_rep_seqs.fasta",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/deblur_output_final/iGEM_16S.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/tax4fun2_working",
             num_threads=40)

# Ocean
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/deblur_output_final/ocean_16S_rep_seqs.fasta",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/deblur_output_final/ocean_16S.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/tax4fun2_working",
             num_threads=40)

# Blueberry
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/deblur_output_exported/blueberry_16S_rep_seqs.fna",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/deblur_output_exported/blueberry_16S.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/tax4fun2_working",
             num_threads=40)

# Indian
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/indian/16S_workflow/deblur_output_final/indian_16S_rep_seqs.fasta",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/indian/16S_workflow/deblur_output_final/indian_16S.renamed.with_extra_header.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/indian/16S_workflow/tax4fun2_working",
             num_threads=40)

# Cameroon
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/cameroon/16S_workflow/deblur_output_final/cameroon_16S_rep_seqs.fasta",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/cameroon/16S_workflow/deblur_output_final/cameroon_16S.renamed.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/cameroon/16S_workflow/tax4fun2_working",
             num_threads=40)

# Primate
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/primate/16S/final_files/primate_asvs.fna",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/primate/16S/final_files/primate_16S.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/primate/16S/tax4fun2_working",
             num_threads=40)

# Soil crossbiome
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S/qiime2_artifacts/soil_16S_rep_seqs.fasta",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S/qiime2_artifacts/soil_16S.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/soil_crossbiome/16S/tax4fun2_working",
             num_threads=40)


# Microbial mat
run_tax4fun2(in_fasta_path = "/home/gavin/projects/picrust_pipeline/data/validation/microbial_mat/QIITA_16S/BIOM/47904/gg_rep_seqs.fna",
             in_tsv_biom_path = "/home/gavin/projects/picrust_pipeline/data/validation/microbial_mat/QIITA_16S/BIOM/47904/otu_table.biom.tsv",
             output_dir = "/home/gavin/projects/picrust_pipeline/data/validation/microbial_mat/QIITA_16S/BIOM/tax4fun2_working",
             num_threads=40)


