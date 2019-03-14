### Comparing # pathways identified as a function of phylogenetic diversity of the samples.

hmp_KEGG_pathways <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime_pipeline/final_otu_tables/otu_table_gg-only_filt_norm_ko_L3.biom.txt",
                                skip=1, comment.char = "", quote="", stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

colnames(hmp_KEGG_pathways) <- gsub(".nonchimera.fasta", "", colnames(hmp_KEGG_pathways))

hmp_KEGG_pathways <- hmp_KEGG_pathways[, -which(colnames(hmp_KEGG_pathways) == "KEGG_Pathways")]

hmp_metacyc_pathways <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/picrust2_pipeline/picrust2_full_output_2.1.0-b/pathways_out/path_abun_unstrat.tsv",
                                   stringsAsFactors = FALSE, header=TRUE, sep="\t", row.names=1)

hmp_faith_pd <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/16S/qiime2_artifacts/diversity/faith_pd_diversity/exported/alpha-diversity.tsv",
                           header=T, sep='\t', row.names=1)


overlapping_col <- colnames(hmp_KEGG_pathways)[which(colnames(hmp_KEGG_pathways) %in% colnames(hmp_metacyc_pathways))]

hmp_faith_pd_subset <- hmp_faith_pd[overlapping_col,, drop=FALSE]

hmp_KEGG_pathways <- hmp_KEGG_pathways[, overlapping_col]
hmp_metacyc_pathways <- hmp_metacyc_pathways[, overlapping_col]

hmp_KEGG_pathways_relab <- data.frame(sweep(hmp_KEGG_pathways, 2, colSums(hmp_KEGG_pathways), '/')) * 100
hmp_metacyc_pathways_relab <- data.frame(sweep(hmp_metacyc_pathways, 2, colSums(hmp_metacyc_pathways), '/')) * 100

hmp_KEGG_pathways_relab <- hmp_KEGG_pathways_relab[-which(rowSums(hmp_KEGG_pathways_relab) == 0),]

par(mfrow=c(1,2))
plot(hmp_faith_pd_subset$faith_pd, colSums(hmp_KEGG_pathways_relab > 0), ylim=c(0, 350), xlim=c(0, 20), pch=16,
     main="PICRUSt1 (KEGG)", ylab="Number of Pathways", xlab="Faith's Phylogenetic Diversity")
plot(hmp_faith_pd_subset$faith_pd, colSums(hmp_metacyc_pathways_relab > 0), ylim=c(0, 350), xlim=c(0, 20), pch=16,
     main="PICRUSt2 (MetaCyc)", ylab="Number of Pathways", xlab="Faith's Phylogenetic Diversity")



### TESTING

library(reshape2)
library(ggplot2)
library(vegan)

hmp_map <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/HMIWGS_healthy.csv",
                      header=T, sep=",", stringsAsFactors = FALSE, comment.char = "")
hmp_map$Body.Site[which(hmp_map$Body.Site == "anterior_nares")] <- "Anterior nares"
hmp_map$Body.Site[which(hmp_map$Body.Site == "buccal_mucosa")] <- "Buccal mucosa"
hmp_map$Body.Site[which(hmp_map$Body.Site == "left_retroauricular_crease")] <- "Left retroauricular crease"
hmp_map$Body.Site[which(hmp_map$Body.Site == "right_retroauricular_crease")] <- "Right retroauricular crease"
hmp_map$Body.Site[which(hmp_map$Body.Site == "posterior_fornix")] <- "Posterior fornix"
hmp_map$Body.Site[which(hmp_map$Body.Site == "stool")] <- "Stool"
hmp_map$Body.Site[which(hmp_map$Body.Site == "supragingival_plaque")] <- "Supragingival plaque"
hmp_map$Body.Site[which(hmp_map$Body.Site == "tongue_dorsum")] <- "Tongue dorsum"

# One time to write out QIIME2 mapfile
# hmp_map <- hmp_map[,c("SRS.ID", "Body.Site")]
# colnames(hmp_map) <- c("sampleid", "bodysite")
# write.table(x = hmp_map, file="/home/gavin/projects/picrust_pipeline/data/validation/hmp/HMIWGS_healthy_qiime2.tsv",
#             sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

tol15rainbow_modified=c("white", "black", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")

bodysite_col =c("peachpuff3", "violetred4", "azure2",     "maroon3",    "goldenrod2", "chocolate3", "pink3",      "khaki3")


hmp_KEGG_pathways_relab_dist <- vegdist(t(hmp_KEGG_pathways_relab), method = "bray")
hmp_KEGG_pathways_relab_dist_NMDS <- metaMDS(hmp_KEGG_pathways_relab_dist)
hmp_KEGG_pathways_relab_dist_NMDS_df <- data.frame(NMDS1=hmp_KEGG_pathways_relab_dist_NMDS$points[,1],
                                                   NMDS2=hmp_KEGG_pathways_relab_dist_NMDS$points[,2],
                                                   bodysite=hmp_map[rownames(hmp_KEGG_pathways_relab_dist_NMDS$points), "Body.Site"])

ggplot(data=hmp_KEGG_pathways_relab_dist_NMDS_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour=bodysite, size=1.5)) +
  theme_minimal() +
  scale_colour_manual(values=bodysite_col) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=6)))


hmp_metacyc_pathways_relab_dist <- vegdist(t(hmp_metacyc_pathways_relab), method = "bray")
hmp_metacyc_pathways_relab_dist_NMDS <- metaMDS(hmp_metacyc_pathways_relab_dist)
hmp_metacyc_pathways_relab_dist_NMDS_df <- data.frame(NMDS1=hmp_metacyc_pathways_relab_dist_NMDS$points[,1],
                                                   NMDS2=hmp_metacyc_pathways_relab_dist_NMDS$points[,2],
                                                   bodysite=hmp_map[rownames(hmp_metacyc_pathways_relab_dist_NMDS$points), "Body.Site"])

ggplot(data=hmp_metacyc_pathways_relab_dist_NMDS_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour=bodysite, size=1.5)) +
  theme_minimal() +
  scale_colour_manual(values=bodysite_col) +
  guides(size=FALSE, fill = guide_legend(override.aes = list(size=6)))








hmp_KEGG_pathways_relab$pathway <- rownames(hmp_KEGG_pathways_relab)
hmp_KEGG_pathways_relab_melt <- melt(hmp_KEGG_pathways_relab, factorsAsStrings=TRUE)
hmp_KEGG_pathways_relab_melt$variable <- as.character(hmp_KEGG_pathways_relab_melt$variable)
hmp_KEGG_pathways_relab_melt$body_site <- hmp_map[hmp_KEGG_pathways_relab_melt$variable, "Body.Site"]
hmp_KEGG_pathways_relab_melt <- hmp_KEGG_pathways_relab_melt[, -which(colnames(hmp_KEGG_pathways_relab_melt) == "variable")]
hmp_KEGG_pathways_relab_melt_by_bodysite <- aggregate(value ~ body_site + pathway, data=hmp_KEGG_pathways_relab_melt, FUN=mean)


hmp_metacyc_pathways_relab$pathway <- rownames(hmp_metacyc_pathways_relab)
hmp_metacyc_pathways_relab_melt <- melt(hmp_metacyc_pathways_relab, factorsAsStrings=TRUE)
hmp_metacyc_pathways_relab_melt$variable <- as.character(hmp_metacyc_pathways_relab_melt$variable)
hmp_metacyc_pathways_relab_melt$body_site <- hmp_map[hmp_metacyc_pathways_relab_melt$variable, "Body.Site"]
hmp_metacyc_pathways_relab_melt <- hmp_metacyc_pathways_relab_melt[, -which(colnames(hmp_metacyc_pathways_relab_melt) == "variable")]
hmp_metacyc_pathways_relab_melt_by_bodysite <- aggregate(value ~ body_site + pathway, data=hmp_metacyc_pathways_relab_melt, FUN=mean)

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

ggplot(hmp_metacyc_pathways_relab_melt_by_bodysite, aes(x=body_site, y=value)) +
  geom_bar(aes(fill = pathway, colour="black"), stat = "identity") +
  guides(fill=FALSE, colour=FALSE) +
  scale_fill_manual(values=sample(color, length(unique(hmp_metacyc_pathways_relab_melt_by_bodysite$pathway))))

ggplot(hmp_KEGG_pathways_relab_melt_by_bodysite, aes(x=body_site, y=value)) +
  geom_bar(aes(fill = pathway, colour="black"), stat = "identity") +
  guides(fill=FALSE, colour=FALSE) +
  scale_fill_manual(values=sample(color, length(unique(hmp_KEGG_pathways_relab_melt_by_bodysite$pathway))))


library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


pie(rep(1,n), col=sample(color, n))

table(hmp_map[colnames(hmp_metacyc_pathways), "Body.Site"])




SRS011584  0.1529007624 -0.0595797718            Posterior fornix
SRS016553  0.0524581788  0.0928000084              Anterior nares
SRS024596  0.0659348171  0.0645333829  Left retroauricular crease
SRS065347  0.1223518241 -0.0674690392            Posterior fornix
SRS024310  0.1218017969 -0.0692571212            Posterior fornix

