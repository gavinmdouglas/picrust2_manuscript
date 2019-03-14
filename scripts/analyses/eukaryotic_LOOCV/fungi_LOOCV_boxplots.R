library(ggplot2)
library(cowplot)
library(reshape2)
library(Hmisc)

source("/home/gavin/gavin_backup/projects/picrust2_manuscript/scripts/picrust2_ms_functions.R")
setwd("/home/gavin/projects/picrust_pipeline/fungal_genomes/LOOCV")

make_tabular_LOOCV_metrics <- function(outfolder, num_subsets, reffile_path) {

  LOOCV_out <- data.frame(matrix(NA, nrow=num_subsets, ncol=7))
  colnames(LOOCV_out) <- c("level", "taxon", "mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")
  
  reffile <- data.frame(t(read.table(reffile_path, header=T, sep="\t", row.names=1)))
  
  LOOCV_subdirs <- list.dirs(path = outfolder, recursive = FALSE)

  rep_i = 1
  
  for (subdir in LOOCV_subdirs) {
    LOOCV_outfiles <- list.files(subdir, full.names = TRUE)
    
    subdir_clean <- sub(".*/", "", subdir)
  
    for(outfile in LOOCV_outfiles) {
      
      outfile_clean <- sub("_loocv.txt", "", outfile)
      outfile_clean <- sub(".*/", "", outfile_clean)
      
      loocv_pred <- data.frame(t(read.table(outfile, header=T, sep="\t", row.names=1)))
    
      # Remove NSTI row.
      nsti_row_i <- which(rownames(loocv_pred) == "metadata_NSTI")
      nsti_row <- as.numeric(loocv_pred[nsti_row_i, ])
      loocv_pred <- loocv_pred[-nsti_row_i, , drop=FALSE]
      
      reffile_subset <- reffile[, colnames(loocv_pred), drop=FALSE]
      
      rho_out <- cor_all_cols(tab1 = loocv_pred, tab2 = reffile_subset, cat_string="LOOCV", metric="spearman")
      
      acc_out <- calc_accuracy_metrics(reffile_subset, loocv_pred, category="LOOCV")
  
      LOOCV_out[rep_i, c("level", "taxon")] <- c(subdir_clean, outfile_clean)
      LOOCV_out[rep_i, c("mean_rho", "mean_acc", "mean_precision", "mean_recall", "mean_nsti")] <- c(mean(rho_out$metric),
                                                                                                      mean(acc_out$acc),
                                                                                                      mean(acc_out$precision),
                                                                                                      mean(acc_out$recall),
                                                                                                      mean(nsti_row))
      rep_i = rep_i + 1
    }
    
  }
  return(LOOCV_out)
}

# Note that the num_subsets is the number of lines in the "groupings" files like fungi_18S_taxa_groupings.tsv
LOOCV_out_18S <- make_tabular_LOOCV_metrics(outfolder="fungi_18S",
                           num_subsets = 761,
                           reffile_path="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_18S_counts.txt")

LOOCV_out_ITS <- make_tabular_LOOCV_metrics(outfolder="fungi_ITS",
                                            num_subsets = 672,
                                            reffile_path="/home/gavin/gavin_backup/projects/picrust2_manuscript/data/reference/mean_func_tables/ec_ITS_counts.txt")


LOOCV_out_18S$level <- capitalize(LOOCV_out_18S$level)

LOOCV_out_18S$level <- factor(LOOCV_out_18S$level, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Assembly"))

LOOCV_out_ITS$level <- capitalize(LOOCV_out_ITS$level)

LOOCV_out_ITS$level <- factor(LOOCV_out_ITS$level, levels=c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Assembly"))




LOOCV_out_18S_plot <- ggplot(LOOCV_out_18S, aes(x=level, y=mean_rho, fill=c("coral3"))) + geom_boxplot() +
xlab(c("Taxonomic Level")) + ggtitle("Fungi 18S EC Numbers") + guides(fill=FALSE) + ylim(c(0, 1.0)) +
  ylab(c("Mean Spearman's Correlation Coefficient")) + scale_fill_manual(values=c("coral3"))


LOOCV_out_ITS_plot <- ggplot(LOOCV_out_ITS, aes(x=level, y=mean_rho, fill=c("dodgerblue3"))) + geom_boxplot() +
  xlab(c("Taxonomic Level")) + ggtitle("Fungi ITS EC Numbers") + guides(fill=FALSE) + ylim(c(0, 1.0)) +
  ylab(c("Mean Spearman's Correlation Coefficient")) + scale_fill_manual(values=c("dodgerblue3"))

plot_grid(LOOCV_out_18S_plot,
          LOOCV_out_ITS_plot,
          labels=c("A", "B"))

wilcox.test(LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Assembly"), "mean_rho"],
            LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Phylum"), "mean_rho"])
#W = 977, p-value = 0.01862

wilcox.test(LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Species"), "mean_rho"],
            LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Phylum"), "mean_rho"])
# W = 919.5, p-value = 0.02256

wilcox.test(LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Genus"), "mean_rho"],
            LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Phylum"), "mean_rho"])
# W = 697.5, p-value = 0.018

wilcox.test(LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Family"), "mean_rho"],
            LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Phylum"), "mean_rho"])
#W = 478, p-value = 0.01932

wilcox.test(LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Order"), "mean_rho"],
            LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Phylum"), "mean_rho"])
#W = 252, p-value = 0.03496

wilcox.test(LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Class"), "mean_rho"],
            LOOCV_out_ITS[which(LOOCV_out_ITS$level=="Phylum"), "mean_rho"])
#W = 132, p-value = 0.07446



wilcox.test(LOOCV_out_18S[which(LOOCV_out_18S$level=="Assembly"), "mean_rho"],
            LOOCV_out_18S[which(LOOCV_out_18S$level=="Phylum"), "mean_rho"])
#W = 1412.5, p-value = 0.001243

wilcox.test(LOOCV_out_18S[which(LOOCV_out_18S$level=="Species"), "mean_rho"],
            LOOCV_out_18S[which(LOOCV_out_18S$level=="Phylum"), "mean_rho"])
# W = 1294.5, p-value = 0.001598

wilcox.test(LOOCV_out_18S[which(LOOCV_out_18S$level=="Genus"), "mean_rho"],
            LOOCV_out_18S[which(LOOCV_out_18S$level=="Phylum"), "mean_rho"])
# W = 996.5, p-value = 0.00201

wilcox.test(LOOCV_out_18S[which(LOOCV_out_18S$level=="Family"), "mean_rho"],
            LOOCV_out_18S[which(LOOCV_out_18S$level=="Phylum"), "mean_rho"])
#W = 685, p-value = 0.00146

wilcox.test(LOOCV_out_18S[which(LOOCV_out_18S$level=="Order"), "mean_rho"],
            LOOCV_out_18S[which(LOOCV_out_18S$level=="Phylum"), "mean_rho"])
#W = 361, p-value = 0.000912

wilcox.test(LOOCV_out_18S[which(LOOCV_out_18S$level=="Class"), "mean_rho"],
            LOOCV_out_18S[which(LOOCV_out_18S$level=="Phylum"), "mean_rho"])
# W = 170, p-value = 0.01568
