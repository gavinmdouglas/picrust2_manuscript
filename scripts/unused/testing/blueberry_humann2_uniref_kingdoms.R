uniref50_tax <- read.table("/home/gavin/databases/uniref50_taxonomy.txt",
                           header=T, sep="\t", row.names=1, stringsAsFactors = FALSE,
                           quote="", comment.char = "")

blueberry_humann2 <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/humann2_final_out/humann2_genefamilies_unstratified.tsv",
                                header=T, sep="\t", row.names=1, stringsAsFactors = FALSE,
                                quote="", comment.char = "")

rownames(blueberry_humann2) <- gsub(": .*$", "", rownames(blueberry_humann2))

overlapping_rows <- rownames(blueberry_humann2)[which(rownames(blueberry_humann2) %in% rownames(uniref50_tax))]

non_overlapping_rows <- rownames(blueberry_humann2)[which(! rownames(blueberry_humann2) %in% rownames(uniref50_tax))]

# Note that ~2/3 of UniRef ids aren't in the latest UniRef database... maybe due to the reduced redundancy?!?

blueberry_humann2_subset <- blueberry_humann2[overlapping_rows,]
uniref50_tax_subset <- uniref50_tax[overlapping_rows,]
rm(uniref50_tax)

blueberry_humann2_subset_percent <- data.frame(sweep(blueberry_humann2_subset, colSums(blueberry_humann2_subset), MARGIN = 2, FUN ="/" ))*100
blueberry_humann2_subset_percent$superkingdom <- uniref50_tax_subset$superkingdom
blueberry_humann2_subset_by_superkingdom <- aggregate(. ~ superkingdom, data=blueberry_humann2_subset_percent, FUN = sum)
rownames(blueberry_humann2_subset_by_superkingdom) <- blueberry_humann2_subset_by_superkingdom$superkingdom
blueberry_humann2_subset_by_superkingdom <- blueberry_humann2_subset_by_superkingdom[, -1]
colnames(blueberry_humann2_subset_by_superkingdom) <- gsub("_Abundance.RPKs", "", colnames(blueberry_humann2_subset_by_superkingdom))



blueberry_metaxa2_summary <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/blueberry_metaxa2_kingdom_summary.txt",
                                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
blueberry_metaxa2_summary_percent <- data.frame(sweep(blueberry_metaxa2_summary, colSums(blueberry_metaxa2_summary), MARGIN = 2, FUN ="/" ))*100

plot(as.numeric(blueberry_humann2_subset_by_superkingdom["Eukaryota",]),
     as.numeric(blueberry_metaxa2_summary_percent["Eukaryota", colnames(blueberry_humann2_subset_by_superkingdom)]))

blueberry_ec_spearman_df <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_ec_spearman_df.rds")
blueberry_ko_spearman_df <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_ko_spearman_df.rds")
blueberry_ec_metrics_df <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_ec_acc_metrics.rds")
blueberry_ko_metrics_df <- readRDS("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/blueberry_acc_metrics.rds")

blueberry_ec_spearman_df_nsti2 <- blueberry_ec_spearman_df[which(blueberry_ec_spearman_df$cat == "NSTI=2"),]
rownames(blueberry_ec_spearman_df_nsti2) <- as.character(blueberry_ec_spearman_df_nsti2$sample_names)

blueberry_ko_spearman_df_nsti2 <- blueberry_ko_spearman_df[which(blueberry_ko_spearman_df$cat == "NSTI=2"),]
rownames(blueberry_ko_spearman_df_nsti2) <- as.character(blueberry_ko_spearman_df_nsti2$sample_names)

blueberry_ec_metrics_df_nsti2 <- blueberry_ec_metrics_df[which(blueberry_ec_metrics_df$category == "NSTI=2"),]
rownames(blueberry_ec_metrics_df_nsti2) <- as.character(blueberry_ec_metrics_df_nsti2$sample)

blueberry_ko_metrics_df_nsti2 <- blueberry_ko_metrics_df[which(blueberry_ko_metrics_df$category == "NSTI=2"),]
rownames(blueberry_ko_metrics_df_nsti2) <- as.character(blueberry_ko_metrics_df_nsti2$sample)

blueberry_metaxa2_summary_percent_edit <- blueberry_metaxa2_summary_percent
colnames(blueberry_metaxa2_summary_percent_edit) <- gsub("^BB", "Bact", colnames(blueberry_metaxa2_summary_percent_edit))
colnames(blueberry_metaxa2_summary_percent_edit) <- gsub("_", ".", colnames(blueberry_metaxa2_summary_percent_edit))

# See how metrics vary with % eukaryotic.
par(mfrow=c(2,1))
plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ec_spearman_df_nsti2)]),
     blueberry_ec_spearman_df_nsti2$metric)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ko_spearman_df_nsti2)]),
     blueberry_ko_spearman_df_nsti2$metric)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ec_metrics_df_nsti2)]),
     blueberry_ec_metrics_df_nsti2$precision)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ec_metrics_df_nsti2)]),
     blueberry_ec_metrics_df_nsti2$recall)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ec_metrics_df_nsti2)]),
     blueberry_ec_metrics_df_nsti2$acc)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ko_metrics_df_nsti2)]),
     blueberry_ko_metrics_df_nsti2$precision)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ko_metrics_df_nsti2)]),
     blueberry_ko_metrics_df_nsti2$recall)

plot(as.numeric(blueberry_metaxa2_summary_percent_edit["Eukaryota", rownames(blueberry_ko_metrics_df_nsti2)]),
     blueberry_ko_metrics_df_nsti2$acc)


# Write out a table of eukaryotic UniRef50 ids as an abundance table so that these can be regrouped to KOs and ECs.
eukaryota_uniref50 <- rownames(uniref50_tax_subset)[which(uniref50_tax_subset$superkingdom == "Eukaryota")]
eukaryota_uniref50_out <- data.frame(eukaryota_uniref50, rep.int(1, length(eukaryota_uniref50)))
colnames(eukaryota_uniref50_out) <- c("# Gene Family", "eukaryotic_presence")

write.table(file = "/home/gavin/databases/test/eukaryota_uniref50_out.tsv", x = eukaryota_uniref50_out,
            quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


### Regrouped with these cmds:
# humann2_regroup_table -i eukaryota_uniref50_out.tsv -c /scratch/db/humann2/utility_mapping/map_pfam_uniref50.txt.gz -o eukaryota_PFAM_out.tsv &
#   humann2_regroup_table -i eukaryota_uniref50_out.tsv -c /scratch/db/humann2/utility_mapping/map_level4ec_uniref50.txt.gz -o eukaryota_level4ec_out.tsv &
#   humann2_regroup_table -ieukaryota_uniref50_out.tsv -c /scratch/db/humann2/utility_mapping/map_ko_uniref50.txt.gz -o eukaryota_KO_out.tsv &

ec_ref <- data.frame(t(read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/mean_func_tables/ec_18S.txt",
                     header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names=FALSE)))
ec_ref$rowSums <- rowSums(ec_ref)

ko_ref <- data.frame(t(read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/mean_func_tables/ko_18S.txt",
                                  header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names=FALSE)))
ko_ref$rowSums <- rowSums(ko_ref)

pfam_ref <- data.frame(t(read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/mean_func_tables/pfam_18S.txt",
                                  header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names=FALSE)))
pfam_ref$rowSums <- rowSums(pfam_ref)

ec_humann2 <- read.table("/home/gavin/databases/test/eukaryota_level4ec_out.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names=FALSE, comment.char="")
ec_humann2 <- ec_humann2[-1,, drop=FALSE]
rownames(ec_humann2) <- gsub("^", "EC:", rownames(ec_humann2))

ko_humann2 <- read.table("/home/gavin/databases/test/eukaryota_KO_out.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names=FALSE, comment.char="")
ko_humann2 <- ko_humann2[-1,, drop=FALSE]


pfam_humann2 <- read.table("/home/gavin/databases/test/eukaryota_PFAM_out.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE, row.names=1, check.names=FALSE, comment.char="")
pfam_humann2 <- pfam_humann2[-1,, drop=FALSE]

overlapping_ec <- rownames(ec_ref)[which(rownames(ec_ref) %in% rownames(ec_humann2))]
ec_humann2_overlap <- ec_humann2[overlapping_ec,, drop=FALSE]
ec_ref_overlap <- ec_ref[overlapping_ec,]

overlapping_ko <- rownames(ko_ref)[which(rownames(ko_ref) %in% rownames(ko_humann2))]
ko_humann2_overlap <- ko_humann2[overlapping_ko,, drop=FALSE]
ko_ref_overlap <- ko_ref[overlapping_ko,]

overlapping_pfam <- rownames(pfam_ref)[which(rownames(pfam_ref) %in% rownames(pfam_humann2))]
pfam_humann2_overlap <- pfam_humann2[overlapping_pfam,, drop=FALSE]
pfam_ref_overlap <- pfam_ref[overlapping_pfam,]

par(mfrow=c(1,3))
plot(ec_humann2_overlap$eukaryotic_presence, ec_ref_overlap$rowSums)
plot(ko_humann2_overlap$eukaryotic_presence, ko_ref_overlap$rowSums)
plot(pfam_humann2_overlap$eukaryotic_presence, pfam_ref_overlap$rowSums)

cor.test(ec_humann2_overlap$eukaryotic_presence, ec_ref_overlap$rowSums)
cor.test(ko_humann2_overlap$eukaryotic_presence, ko_ref_overlap$rowSums)
cor.test(pfam_humann2_overlap$eukaryotic_presence, pfam_ref_overlap$rowSums)
