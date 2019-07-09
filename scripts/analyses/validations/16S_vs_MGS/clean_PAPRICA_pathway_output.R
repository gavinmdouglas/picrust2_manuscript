### Make PAPRICA pathway abundance output consistent with other tables.

path_description <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/description_mapfiles/metacyc_pathways_info_prokaryotes.txt.gz"),
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


clean_write_paprica_pathabun(infile="/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/paprica_out/hmp_paprica_out.path_tally.csv",
                             outfile="/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/paprica_out/hmp_paprica_out.path_tally_clean.csv",
                             path_descrip=path_description)
                             

clean_write_paprica_pathabun(infile="/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/paprica_out/iGEM_paprica_out.path_tally.csv",
                             outfile="/home/gavin/projects/picrust_pipeline/data/validation/iGEM/16S/paprica_out/iGEM_paprica_out.path_tally_clean.csv",
                             path_descrip=path_description)

clean_write_paprica_pathabun(infile="/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/paprica_out/ocean_paprica_out.path_tally.csv",
                             outfile="/home/gavin/projects/picrust_pipeline/data/validation/ocean/16S/paprica_out/ocean_paprica_out.path_tally_clean.csv",
                             path_descrip=path_description)

clean_write_paprica_pathabun(infile="/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/paprica_out/blueberry_paprica_out.path_tally.csv",
                             outfile="/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/paprica_out/blueberry_paprica_out.path_tally_clean.csv",
                             path_descrip=path_description)


# Also get all possible PAPRICA MetaCyc pathways predicted using similar approach:
paprica_ref_pathabun <- read.table("/home/gavin/local/anaconda2_old_backup/envs/paprica/download/paprica/ref_genome_database/bacteria/terminal_paths.csv",
                                   header=TRUE, sep=",", row.names=1, check.names=FALSE, stringsAsFactors = FALSE, comment.char="", quote="\"")

colnames(paprica_ref_pathabun)[which(! colnames(paprica_ref_pathabun) %in% rownames(path_description))]

paprica_ref_pathabun_subset <- paprica_ref_pathabun[, which(colnames(paprica_ref_pathabun) %in% rownames(path_description))]

write.table(x = path_description[colnames(paprica_ref_pathabun_subset), "V1"],
            file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/16S_validation/possible_metacyc_pathways/paprica_bacteria_pathways.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)



