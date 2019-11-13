### Make PAPRICA pathway abundance output consistent with other tables.

rm(list=ls(all=TRUE))

path_description <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz"),
                               header=FALSE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")

rownames(path_description) <- path_description$V2



clean_write_paprica_pathabun <- function(infile, outfile, path_descrip) {
  paprica_pathabun <- data.frame(t(read.table(infile,header=TRUE, row.names=1, sep=",", stringsAsFactors = FALSE, check.names=FALSE, quote="\"", comment.char="")),
                                 check.names=FALSE)
  colnames(paprica_pathabun) <- gsub("_rerep_seqs.", "", colnames(paprica_pathabun))
  
  paprica_pathabun <- paprica_pathabun[which(rownames(paprica_pathabun) %in% path_descrip$V2), ]
  
  rownames(paprica_pathabun) <- path_descrip[rownames(paprica_pathabun), "V1"]
  
  write.table(x = paprica_pathabun, file = outfile, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
}

datasets_to_paprica_path_out <- list()

datasets_to_paprica_path_out[["cameroon"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/cameroon/16S_workflow/paprica_working/combined.path_tally.csv"
datasets_to_paprica_path_out[["hmp"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/paprica_working/combined.path_tally.csv"
datasets_to_paprica_path_out[["indian"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/indian/16S_workflow/paprica_working/combined.path_tally.csv"
datasets_to_paprica_path_out[["primate"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/primate/16S/paprica_working/combined.path_tally.csv"
datasets_to_paprica_path_out[["mammal"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/paprica_working/combined.path_tally.csv"
datasets_to_paprica_path_out[["ocean"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/paprica_working/combined.path_tally.csv"
datasets_to_paprica_path_out[["blueberry"]] <- "/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/paprica_working/combined.path_tally.csv"

datasets <- c("cameroon", "hmp", "indian", "primate", "mammal", "ocean", "blueberry")

for(d in datasets) {
  
  clean_write_paprica_pathabun(infile=datasets_to_paprica_path_out[[d]],
                               outfile=paste("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/paprica_out/", d, "_paprica_path.tsv", sep=""),
                               path_descrip=path_description)
}


# Also get all possible PAPRICA MetaCyc pathways predicted using similar approach:
paprica_ref_pathabun <- read.table(gzfile("/home/gavin/local/prg/paprica_v0.5.2/ref_genome_database/bacteria/terminal_paths.csv.gz"),
                                   header=TRUE, sep=",", row.names=1, check.names=FALSE, stringsAsFactors = FALSE, comment.char="", quote="\"")

colnames(paprica_ref_pathabun)[which(! colnames(paprica_ref_pathabun) %in% rownames(path_description))]

paprica_ref_pathabun_subset <- paprica_ref_pathabun[, which(colnames(paprica_ref_pathabun) %in% rownames(path_description))]

write.table(x = path_description[colnames(paprica_ref_pathabun_subset), "V1"],
            file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_path/paprica_path.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)



