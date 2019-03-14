library(ggplot2)
library(cowplot)

setwd("/home/gavin/gavin_backup/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# HMP
hmp_rho <- readRDS(file = "hmp_ko_spearman_df.rds")
hmp_rho_nsti2 <- hmp_rho[which(hmp_rho$cat == "NSTI=2"),]
rownames(hmp_rho_nsti2) <- hmp_rho_nsti2$sample_names

hmp_mgs <- read.table("../../data/mgs_validation/hmp/humann2_ko_unstrat.tsv", header=T, sep="\t", row.names=1)
hmp_mgs <- hmp_mgs[-which(rownames(hmp_mgs) %in% c("UNMAPPED", "UNGROUPED")), rownames(hmp_rho_nsti2)]
hmp_rho_nsti2$mgs_RPK <- colSums(hmp_mgs)

# In the case of the HMP dataset colour the dots by bodysite.
hmp_map <- read.table("/home/gavin/projects/picrust_pipeline/data/validation/hmp/HMIWGS_healthy.csv",
                      header=T, sep=",", stringsAsFactors = FALSE, comment.char = "")

rownames(hmp_map) <- hmp_map$SRS.ID

hmp_rho_nsti2$body_site <- hmp_map[as.character(hmp_rho_nsti2$sample_names), "Body.Site"]
hmp_rho_nsti2$body_site <- gsub("_", " ", hmp_rho_nsti2$body_site)
hmp_rho_nsti2$body_site <- sapply(hmp_rho_nsti2$body_site, capwords)

qual_col <- c("brown",
"chartreuse4",
"orange",
"black",
"dark grey",
"blue")

hmp_rho_depth_plot <- ggplot(hmp_rho_nsti2, aes(x=mgs_RPK, y=metric, color=body_site)) +
  geom_point(size=3) + scale_color_manual(name="Body Site", values=qual_col) +
  xlab("Number of Reads Mapped to KOs per Kilobase") +
  ylab("Spearman Correlation Coefficient") +
  ggtitle("HMP") +
  xlim(200000, 8000000) +
  ylim(0.7, 1) + theme(legend.position="none")


# Mammal
mammal_rho <- readRDS(file = "mammal_ko_spearman_df.rds")
mammal_rho_nsti2 <- mammal_rho[which(mammal_rho$cat == "NSTI=2"),]
rownames(mammal_rho_nsti2) <- mammal_rho_nsti2$sample_names

mammal_mgs <- read.table("../../data/mgs_validation/mammal/humann2_ko_unstrat.tsv", header=T, sep="\t", row.names=1)
mammal_mgs <- mammal_mgs[-which(rownames(mammal_mgs) %in% c("UNMAPPED", "UNGROUPED")), rownames(mammal_rho_nsti2)]
mammal_rho_nsti2$mgs_RPK <- colSums(mammal_mgs)

mammal_rho_depth_plot <- ggplot(mammal_rho_nsti2, aes(x=mgs_RPK, y=metric)) +
  geom_point(size=3) + #scale_color_manual(name="Body Site", values=qual_col) +
  xlab("Number of Reads Mapped to KOs per Kilobase") +
  ylab("Spearman Correlation Coefficient") +
  ggtitle("Mammal") +
  xlim(200000, 700000) +
  ylim(0.7, 1)


# Ocean
ocean_rho <- readRDS(file = "ocean_ko_spearman_df.rds")
ocean_rho_nsti2 <- ocean_rho[which(ocean_rho$cat == "NSTI=2"),]
rownames(ocean_rho_nsti2) <- ocean_rho_nsti2$sample_names

ocean_mgs <- read.table("../../data/mgs_validation/ocean/humann2_ko_unstrat.tsv", header=T, sep="\t", row.names=1)
ocean_mgs <- ocean_mgs[-which(rownames(ocean_mgs) %in% c("UNMAPPED", "UNGROUPED")), rownames(ocean_rho_nsti2)]
ocean_rho_nsti2$mgs_RPK <- colSums(ocean_mgs)

ocean_rho_depth_plot <- ggplot(ocean_rho_nsti2, aes(x=mgs_RPK, y=metric)) +
  geom_point(size=3) + #scale_color_manual(name="Body Site", values=qual_col) +
  xlab("Number of Reads Mapped to KOs per Kilobase") +
  ylab("Spearman Correlation Coefficient") +
  ggtitle("Ocean") +
  xlim(3000000, 12000000) +
  ylim(0.7, 1)


# Soil (Blueberry)
blueberry_rho <- readRDS(file = "blueberry_ko_spearman_df.rds")
blueberry_rho_nsti2 <- blueberry_rho[which(blueberry_rho$cat == "NSTI=2"),]
rownames(blueberry_rho_nsti2) <- blueberry_rho_nsti2$sample_names

blueberry_mgs <- read.table("../../data/mgs_validation/blueberry/humann2_ko_unstrat.tsv", header=T, sep="\t", row.names=1)
blueberry_mgs <- blueberry_mgs[-which(rownames(blueberry_mgs) %in% c("UNMAPPED", "UNGROUPED")), rownames(blueberry_rho_nsti2)]
blueberry_rho_nsti2$mgs_RPK <- colSums(blueberry_mgs)

blueberry_rho_depth_plot <- ggplot(blueberry_rho_nsti2, aes(x=mgs_RPK, y=metric)) +
  geom_point(size=3) + #scale_color_manual(name="Body Site", values=qual_col) +
  xlab("Number of Reads Mapped to KOs per Kilobase") +
  ylab("Spearman Correlation Coefficient") +
  ggtitle("Soil (Blueberry)") +
  xlim(0, 400000) +
  ylim(0.7, 1)


plot_grid(hmp_rho_depth_plot, mammal_rho_depth_plot,
          ocean_rho_depth_plot, blueberry_rho_depth_plot,
          labels=c("A", "B", "C", "D"))
