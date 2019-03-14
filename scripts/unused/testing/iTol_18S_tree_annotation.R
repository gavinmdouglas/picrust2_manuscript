library(ape)
library(RColorBrewer)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/testing/")

euk_taxa <- read.table("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/db_taxa/all_RefSeq_taxa_levels_30Oct2018_eukaryota_no-metazoa.tsv",
                       header=T, sep="\t", comment.char="", stringsAsFactors = FALSE)
# Remove duplicated line.
euk_taxa <- euk_taxa[-118,]
rownames(euk_taxa) <- euk_taxa$taxid

refseq_summary <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/assembly_summary_refseq_30Oct2018.txt",
                       header=T, sep="\t", comment.char="", stringsAsFactors = FALSE, quote="", skip=1, row.names=1)

tree_in <- read.tree("/home/gavin/projects/picrust_pipeline/RefSeq_redownloaded/euk_alignments/raxml_18S/18S_ssu_align.eukarya.mask.derep.raxml.bestTree")
tip_taxid <- as.character(refseq_summary[gsub("_cluster", "", tree_in$tip.label), "taxid"])
tip_class <- euk_taxa[tip_taxid, "Class"]
tip_class_factor_level <- as.numeric(factor(tip_class))
tip_class_factor_level[which(is.na(tip_class_factor_level))] <- 34

# Output different colours for each genome within a given class Tips without a class should be set as NA.
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_selection <- c(sample(col_vector, 33), "#000000")

tip_class_factor_level_colours <- tolower(col_selection[tip_class_factor_level])

tip_class_factor_level_colours_pasted <- paste(tree_in$tip.label, tip_class_factor_level_colours, sep=" ")

fileConn <- file("18S_tree_class_iToL_colorstrip.txt")
writeLines(c("DATASET_COLORSTRIP",
              "SEPARATOR SPACE",
              "DATASET_LABEL 18S_class_colorstrip",
              "COLOR #ff0000",
              "COLOR_BRANCHES 1",
              "STRIP_WIDTH 25",
              "MARGIN 0",
              "BORDER_WIDTH 1",
              "BORDER_COLOR #000",
              "SHOW_INTERNAL 0",
              "DATA",
              tip_class_factor_level_colours_pasted),
           sep="\n", con=fileConn)

close(fileConn)

