### Subset pathway mapfile to pathways that are found in reference genomes at least once.

pathabun_18S <- data.frame(t(read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/metacyc_output/ref_genome_metacyc_18S/path_abun_unstrat.tsv",
                                        row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

pathabun_ITS <- data.frame(t(read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/metacyc_output/ref_genome_metacyc_ITS/path_abun_unstrat.tsv",
                                        row.names=1, header=T, sep="\t", stringsAsFactors = FALSE)), check.names=FALSE)

# Get union of pathways.
detected_pathways <- unique(c(colnames(pathabun_18S), colnames(pathabun_ITS)))

fungi_mapfile <- read.table("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt",
                            header=FALSE, sep="\t", comment.char="", quote="", stringsAsFactors = FALSE)

fungi_mapfile_subset <- fungi_mapfile[which(fungi_mapfile$V1 %in% detected_pathways), ]

write.table(x = fungi_mapfile_subset, file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi_present.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")