### Commands to create tables of the mean EC number abundances across all reference sequences.
### These tables will then be run with pathway_pipeline.py for a different way of getting the "expected"
### MetaCyc pathway abundances (the alternative was to sum over the MetaCyc pathway abundances predicted within
### each reference genome).

ec <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ec.txt.gz"),
                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

ec_t <- data.frame(t(ec), check.names=FALSE)

ec_t_rowMeans <- rowMeans(ec_t)

ec_means_out <- data.frame(mean_EC=ec_t_rowMeans)
rownames(ec_means_out) <- names(ec_t_rowMeans)

write.table(x = ec_means_out, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_16S.tsv",
            col.names=NA, row.names=TRUE, quote=FALSE, sep="\t")


ec_18S <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/fungi/ec_18S_counts.txt.gz"),
                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

ec_18S_t <- data.frame(t(ec_18S), check.names=FALSE)

ec_18S_t_rowMeans <- rowMeans(ec_18S_t)

ec_18S_means_out <- data.frame(mean_EC=ec_18S_t_rowMeans)
rownames(ec_18S_means_out) <- names(ec_18S_t_rowMeans)

write.table(x = ec_18S_means_out, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_18S.tsv",
            col.names=NA, row.names=TRUE, quote=FALSE, sep="\t")



ec_ITS <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/fungi/ec_ITS_counts.txt.gz"),
                     header=TRUE, sep="\t", row.names=1, check.names=FALSE)

ec_ITS_t <- data.frame(t(ec_ITS), check.names=FALSE)

ec_ITS_t_rowMeans <- rowMeans(ec_ITS_t)

ec_ITS_means_out <- data.frame(mean_EC=ec_ITS_t_rowMeans)
rownames(ec_ITS_means_out) <- names(ec_ITS_t_rowMeans)

write.table(x = ec_ITS_means_out, file = "/home/gavin/gavin_backup/projects/picrust2_manuscript/data/working_tables/ref_wide_mean_ec/mean_ec_ref_ITS.tsv",
            col.names=NA, row.names=TRUE, quote=FALSE, sep="\t")


### "function" was added as the first column name to all of these tables manually.

### Ran these commands on the above tables to get predicted MetaCyc pathway abundances:

#pathway_pipeline.py -i mean_ec_ref_16S.tsv -o mean_ec_ref_16S_pathway --coverage
#pathway_pipeline.py -i mean_ec_ref_ITS.tsv -m /home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt --coverage -o mean_ec_ref_ITS_pathway
#pathway_pipeline.py -i mean_ec_ref_18S.tsv -m /home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/pathway_mapfiles/metacyc_path2rxn_struc_filt_fungi.txt --coverage -o mean_ec_ref_18S_pathway
