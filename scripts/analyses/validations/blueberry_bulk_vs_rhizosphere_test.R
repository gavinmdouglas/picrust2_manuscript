setwd("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/18S/picrust2_full_output/")

blue_18S_ec <- read.table("EC_metagenome_out_nsti_2.0/pred_metagenome_unstrat.tsv", header=T, sep="\t", row.names=1, check.names=FALSE)
blue_18S_pathabun <- read.table("pathways_out/path_abun_unstrat.tsv", header=T, sep="\t", row.names=1, check.names=FALSE)

blue_meta <- read.table("../../Blueberry_18S_metadata.tsv", header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, row.names=1)
blue_meta <- blue_meta[-1, ]
blue_meta <- blue_meta[colnames(blue_18S_ec),]

blue_meta <- blue_meta[which(blue_meta$Description %in% c("Bulk", "Rhizosphere")), ]


blue_18S_ec_relab <- data.frame(t(sweep(blue_18S_ec, 2, colSums(blue_18S_ec), `/`) * 100), check.names=TRUE)
blue_18S_pathabun_relab <- data.frame(t(sweep(blue_18S_pathabun, 2, colSums(blue_18S_pathabun), `/`) * 100), check.names=TRUE)

blue_18S_ec_relab$group <- blue_meta[rownames(blue_18S_ec_relab), "Description_1"]
blue_18S_pathabun_relab$group <- blue_meta[rownames(blue_18S_pathabun_relab), "Description_1"]

p_values <- c()
for(func in grep("PWY", colnames(blue_18S_pathabun_relab), value=TRUE)) {

  formula_tmp <- as.formula(paste(func, "~", "group"))

  wilcox_p <- wilcox.test(formula_tmp, data=blue_18S_pathabun_relab)$p.value

  p_values <- c(p_values, wilcox_p)
}


