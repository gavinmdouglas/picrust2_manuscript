tmp_meta <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/HMIWGS_healthy.csv", header=T, sep=",", comment.char="")
rownames(tmp_meta) <- tmp_meta$SRS.ID

tmp <- hmp_ko_picrust2_nsti2_gg_vs_mgs
rownames(tmp) <- tmp$sample_names
head(tmp[order(tmp$metric),])
tmp$mgs_depth <- colSums(hmp_infiles$all_kos_overlap$mgs_ko)[rownames(tmp)]
tmp$body_site <- tmp_meta[rownames(tmp), "Body.Site"]
plot(tmp$mgs_depth, tmp$metric, col=tmp$body_site, pch=16)