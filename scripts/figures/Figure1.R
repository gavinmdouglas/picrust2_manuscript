### Figure 1A will be flowchart figure (made in powerpoint).
### Figure 1B will be percent of all ASVs/OTUs used for predictions when using denoising approach vs reference-based OTU picking.
### Figure 1C will be barplot of number of taxa.

rm(list=ls(all=TRUE))

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/data/")

library(ggplot2)
library(reshape2)
library(cowplot)
library(magick)
library(pdftools)

# Panel A.
flowchart_plot <- ggdraw() +
  draw_image("Fig1_Flowchart/gavinflowchart_mar11.png")

# Panel B.
OTU_counts <- read.table("working_tables/16S_validation_denovo_OTU_counts.tsv",
                         header=T, sep="\t", stringsAsFactors = FALSE)

# Get percents of OTUs and ASVs that passed PICRUSt2.
OTU_counts$OTU_percent_passed <- (OTU_counts$Ref.OTUs.passed.PICRUSt2/OTU_counts$Final.OTUs..pre.GG.filtering.)*100
OTU_counts$ASV_percent_passed <- (OTU_counts$Final.ASVs.passed.PICRUSt2/OTU_counts$Final.ASVs)*100

OTU_counts_subset <- OTU_counts[,c("Dataset", "OTU_percent_passed", "ASV_percent_passed")]

OTU_counts_subset_melt <- melt(OTU_counts_subset, variable.name="Datatype")

OTU_counts_subset_melt$Datatype <- as.character(OTU_counts_subset_melt$Datatype)

OTU_counts_subset_melt[which(OTU_counts_subset_melt$Datatype == "OTU_percent_passed"), "Datatype"] <- "Ref. OTUs"
OTU_counts_subset_melt[which(OTU_counts_subset_melt$Datatype == "ASV_percent_passed"), "Datatype"] <- "ASVs"
OTU_counts_subset_melt$Datatype <- factor(OTU_counts_subset_melt$Datatype,
                                          levels=c("Ref. OTUs", "ASVs"))

OTU_counts_subset_melt$Dataset <- gsub("Soil \\(Blueberry\\)", "Soil", OTU_counts_subset_melt$Dataset)

percent_remaining_plot <- ggplot(OTU_counts_subset_melt, aes(x=Datatype, y=value, fill=Dataset)) +
                                 geom_bar(stat="identity", position="dodge") +
                                 xlab("") +
                                 ylab("% ASVs/OTUs Remain") +
                                 scale_fill_manual(values=c("dark orange", "brown", "#f4a582", "cornflowerblue")) +
                                 guides(fill=guide_legend(ncol=2)) +
                                 theme(legend.position = c(0.02, 0.9),
                                        #legend.background = element_rect(color = "black", fill = "white", size = 0.4, linetype = "solid"),
                                        axis.text.x=element_text(angle=45, hjust=1, size=10),
                                        axis.text.y=element_text(size=10),
                                        axis.title.y=element_text(size=10),
                                        legend.text = element_text(colour="black", size = 8),
                                        legend.title = element_text(colour="black", size = 9))
# Panel C.
db_taxa_counts <- read.table("db_taxa/PICRUSt2_vs_PICRUSt1_taxa_counts.tsv",
                             header=T, sep="\t", stringsAsFactors = FALSE)
db_taxa_counts <- db_taxa_counts[, -which(colnames(db_taxa_counts) == "Superkingdom")]

colnames(db_taxa_counts)[1] <- "Database"

taxa_levels = colnames(db_taxa_counts)[2:ncol(db_taxa_counts)]

db_taxa_counts_melt <- melt(db_taxa_counts)

db_taxa_counts_melt$variable <- factor(db_taxa_counts_melt$variable, levels=taxa_levels)

taxa_barplot <- ggplot(db_taxa_counts_melt, aes(x=variable, y=value, fill=Database)) +
                        geom_bar(stat="identity", position="dodge") +
                        xlab("") +
                        ylab("Number of Taxa") +
                        scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
                        geom_text(aes(x=variable, y=value + 700, label = value),
                                  position = position_dodge(.9), hjust=.5, size=2.5) +
                        theme(legend.position = c(0.02, 0.9),
                              #legend.background = element_rect(color = "black", fill = "white", size = 0.4, linetype = "solid"),
                              axis.text.x=element_text(angle=45, hjust=1, size=10),
                              axis.text.y=element_text(size=10),
                              axis.title.y=element_text(size=10),
                              legend.text = element_text(colour="black", size = 8),
                              legend.title = element_text(colour="black", size = 9))

### Plot final figure.
bottom_row <- plot_grid(percent_remaining_plot, taxa_barplot,
                        labels = c('B', 'C'), ncol=2, nrow=1)
setTimeLimit(); setSessionTimeLimit()

#10x10
plot_grid(flowchart_plot,
          bottom_row,
          labels = c('A', ''),
          ncol=1,
          nrow=2,
          rel_heights=c(2, 1))
