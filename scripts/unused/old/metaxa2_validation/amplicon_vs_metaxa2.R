### Comparisons of taxonomic composition from direct amplicons and rRNA sequences parsed from MGS.
### Note that only MGS rRNA that could be classified down to D_1 were retained (and that were in the expected domain).

library(stringr)
library(reshape2)
library(ggplot2)
library(cowplot)

split_7_level_tax <- function(in_df) {
  
  new_df <- data.frame(matrix(NA, nrow=nrow(in_df), ncol=8))
  colnames(new_df) <- c("D0", "D1", "D2", "D3", "D4", "D5", "D6", "confidence")
  rownames(new_df) <- in_df$Feature.ID
  new_df$confidence <- in_df$confidence
  
  tax_levels <- str_split(in_df$Taxon, ";")
  
  for(i in 1:nrow(in_df)) {
    out_taxa <- tax_levels[[i]]
    while(length(out_taxa) < 7) {
      out_taxa <- c(out_taxa, "Unclassified")
    }
    
    new_df[i, c("D0", "D1", "D2", "D3", "D4", "D5", "D6")] <- out_taxa
  }
  
  return(new_df)
  
}

amplicon_tax_16S <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/deblur_output/taxa/exported/taxonomy.tsv",
                               header = T, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
amplicon_biom_16S <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/16S/deblur_output_exported/blueberry_16S.biom.tsv",
                                header = T, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", skip=1, row.names=1)

mgs_tax_16S <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/metaxa2_out/picrust2_prep/exported/16S_taxonomy.tsv",
                               header = T, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
mgs_biom_16S <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/blueberry/mgs/metaxa2_out/picrust2_prep/pro_input.tsv",
                                header = T, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names=1)

# Rename columns to all start with "BB".
colnames(amplicon_biom_16S) <- gsub("Bact", "BB", colnames(amplicon_biom_16S))

# Split taxonomy into 7 different levels.
amplicon_tax_16S_split <- split_7_level_tax(amplicon_tax_16S)
mgs_tax_16S_split <- split_7_level_tax(mgs_tax_16S)

# Filter out ASVs if they are not classified as Bacteria or if they are unclassified at the order level.
amplicon_tax_16S_split_filt <- amplicon_tax_16S_split[grep("Bacteria", amplicon_tax_16S_split$D0) ,]
mgs_tax_16S_split_filt <- mgs_tax_16S_split[grep("Bacteria", mgs_tax_16S_split$D0) ,]

amplicon_tax_16S_split_filt <- amplicon_tax_16S_split_filt[grep("D_3_", amplicon_tax_16S_split_filt$D3) ,]
mgs_tax_16S_split_filt <- mgs_tax_16S_split_filt[grep("D_3_", mgs_tax_16S_split_filt$D3) ,]

# Convert raw counts to relative abundance.
amplicon_biom_16S_subset <- amplicon_biom_16S[rownames(amplicon_tax_16S_split_filt),]
amplicon_biom_16S_subset_relab <- sweep(amplicon_biom_16S_subset, 2,colSums(amplicon_biom_16S_subset), "/") * 100
amplicon_tax_16S_sum <- aggregate_relabun_by_7_levels(amplicon_biom_16S_subset_relab, amplicon_tax_16S_split_filt)

mgs_biom_16S_subset <- mgs_biom_16S[rownames(mgs_tax_16S_split_filt),]
mgs_biom_16S_subset_relab <- sweep(mgs_biom_16S_subset, 2,colSums(mgs_biom_16S_subset), "/") * 100
mgs_tax_16S_sum <- aggregate_relabun_by_7_levels(mgs_biom_16S_subset_relab, mgs_tax_16S_split_filt)

# Remove all taxa not found at abundance of at least 5% in 1 sample.
# Add an "other category".
amplicon_tax_16S_sum_D3 <- amplicon_tax_16S_sum$D3
amplicon_16S_D3_5per <- which(rowSums(amplicon_tax_16S_sum_D3[, -1] > 5) > 0)
amplicon_tax_16S_sum_D3_subset <- amplicon_tax_16S_sum_D3[amplicon_16S_D3_5per,]
amplicon_tax_16S_sum_D3_other <- data.frame(matrix(NA, nrow=1, ncol=ncol(amplicon_tax_16S_sum_D3)))
colnames(amplicon_tax_16S_sum_D3_other) <- colnames(amplicon_tax_16S_sum_D3)
amplicon_tax_16S_sum_D3_other$D3 <- "Other"
amplicon_tax_16S_sum_D3_other[,2:ncol(amplicon_tax_16S_sum_D3_other)] <- 100 - colSums(amplicon_tax_16S_sum_D3_subset[,-1])
amplicon_tax_16S_sum_D3_subset <- rbind(amplicon_tax_16S_sum_D3_subset, amplicon_tax_16S_sum_D3_other)
amplicon_tax_16S_sum_D3_subset_melt <- melt(amplicon_tax_16S_sum_D3_subset, id.vars = "D3")

mgs_tax_16S_sum_D3 <- mgs_tax_16S_sum$D3
mgs_16S_D3_5per <- which(rowSums(mgs_tax_16S_sum_D3[, -1] > 5) > 0)
mgs_tax_16S_sum_D3_subset <- mgs_tax_16S_sum_D3[mgs_16S_D3_5per,]
mgs_tax_16S_sum_D3_other <- data.frame(matrix(NA, nrow=1, ncol=ncol(mgs_tax_16S_sum_D3)))
colnames(mgs_tax_16S_sum_D3_other) <- colnames(mgs_tax_16S_sum_D3)
mgs_tax_16S_sum_D3_other$D3 <- "Other"
mgs_tax_16S_sum_D3_other[,2:ncol(mgs_tax_16S_sum_D3_other)] <- 100 - colSums(mgs_tax_16S_sum_D3_subset[,-1])
mgs_tax_16S_sum_D3_subset <- rbind(mgs_tax_16S_sum_D3_subset, mgs_tax_16S_sum_D3_other)
mgs_tax_16S_sum_D3_subset_melt <- melt(mgs_tax_16S_sum_D3_subset, id.vars = "D3")

qual_col <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
              "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
              "#5A0007", "#809693", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
              "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
              "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
              "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
              "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
              "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
              "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
              "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
              "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
              "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

# Get all unique D3 names so that factor levels can be made consistently.
amplicon_mgs_D3 <- sort(unique(c(amplicon_tax_16S_sum_D3_subset_melt$D3, mgs_tax_16S_sum_D3_subset_melt$D3)))

amplicon_tax_16S_sum_D3_subset_melt$D3 <- factor(amplicon_tax_16S_sum_D3_subset_melt$D3, levels=amplicon_mgs_D3)
mgs_tax_16S_sum_D3_subset_melt$D3 <- factor(mgs_tax_16S_sum_D3_subset_melt$D3, levels=amplicon_mgs_D3)

amplicon_tax_16S_sum_D3_subset_melt$variable <- gsub("\\.", "_", as.character(amplicon_tax_16S_sum_D3_subset_melt$variable))
mgs_tax_16S_sum_D3_subset_melt$variable <- as.character(mgs_tax_16S_sum_D3_subset_melt$variable)

amplicon_tax_16S_sum_D3_subset_plot <- ggplot(data = amplicon_tax_16S_sum_D3_subset_melt, mapping=aes(x=variable, y=value)) + geom_col(aes(fill=D3)) +
  scale_fill_manual(values=qual_col[as.numeric(sort(unique(amplicon_tax_16S_sum_D3_subset_melt$D3)))]) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative Abundance") + xlab("Sample")

mgs_tax_16S_sum_D3_subset_plot <- ggplot(data = mgs_tax_16S_sum_D3_subset_melt, mapping=aes(x=variable, y=value)) + geom_col(aes(fill=D3)) +
  scale_fill_manual(values=qual_col[as.numeric(sort(unique(mgs_tax_16S_sum_D3_subset_melt$D3)))]) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Relative Abundance") + xlab("Sample")

plot_grid(amplicon_tax_16S_sum_D3_subset_plot, mgs_tax_16S_sum_D3_subset_plot)



aggregate_relabun_by_7_levels <- function(in_biom, in_taxa_split) {
 
  # First subset to overlapping features only.
  overlapping_features <- rownames(in_taxa_split)[which(rownames(in_taxa_split) %in% rownames(in_biom))]
  in_biom <- in_biom[overlapping_features,]
  in_taxa_split <- in_taxa_split[overlapping_features,]

  D0_biom <- in_biom
  D0_biom$D0 <- in_taxa_split$D0
  D0_biom_sum <- aggregate(. ~ D0, FUN=sum, data=D0_biom)
  
  D1_biom <- in_biom
  D1_biom$D1 <- in_taxa_split$D1
  D1_biom_sum <- aggregate(. ~ D1, FUN=sum, data=D1_biom)
  
  D2_biom <- in_biom
  D2_biom$D2 <- in_taxa_split$D2
  D2_biom_sum <- aggregate(. ~ D2, FUN=sum, data=D2_biom)
  
  D3_biom <- in_biom
  D3_biom$D3 <- in_taxa_split$D3
  D3_biom_sum <- aggregate(. ~ D3, FUN=sum, data=D3_biom)
  
  D4_biom <- in_biom
  D4_biom$D4 <- in_taxa_split$D4
  D4_biom_sum <- aggregate(. ~ D4, FUN=sum, data=D4_biom)
  
  D5_biom <- in_biom
  D5_biom$D5 <- in_taxa_split$D5
  D5_biom_sum <- aggregate(. ~ D5, FUN=sum, data=D5_biom)
  
  D6_biom <- in_biom
  D6_biom$D6 <- in_taxa_split$D6
  D6_biom_sum <- aggregate(. ~ D6, FUN=sum, data=D6_biom)
  
  return(list(D0=D0_biom_sum,
              D1=D1_biom_sum,
              D2=D2_biom_sum,
              D3=D3_biom_sum,
              D4=D4_biom_sum,
              D5=D5_biom_sum,
              D6=D6_biom_sum))

  }
