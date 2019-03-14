### Plot 16S vs MGS validations.

setwd("/home/gavin/projects/picrust2_manuscript/saved_RDS/16S_vs_MGS_metrics/")

library(sinaplot)
library(ggplot2)

# Plot KO validations 16S vs MGS.
hmp_ko_cosine <- readRDS("hmp_ko_cosine_df.rds")
hmp_ko_spearman <- readRDS("hmp_ko_spearman_df.rds")
hmp_ko_cosine_gg97 <- readRDS("hmp_ko_picrust2_gg97_vs_mgs_cosine.rds")
hmp_ko_spearman_gg97 <- readRDS("hmp_ko_picrust2_gg97_vs_mgs_spearman.rds")

mammal_ko_cosine <- readRDS("mammal_ko_cosine_df.rds")
mammal_ko_spearman <- readRDS("mammal_ko_spearman_df.rds")
mammal_ko_cosine_gg97 <- readRDS("mammal_ko_picrust2_gg97_vs_mgs_cosine.rds")
mammal_ko_spearman_gg97 <- readRDS("mammal_ko_picrust2_gg97_vs_mgs_spearman.rds")

ocean_ko_cosine <- readRDS("ocean_ko_cosine_df.rds")
ocean_ko_spearman <- readRDS("ocean_ko_spearman_df.rds")
ocean_ko_cosine_gg97 <- readRDS("ocean_ko_picrust2_gg97_vs_mgs_cosine.rds")
ocean_ko_spearman_gg97 <- readRDS("ocean_ko_picrust2_gg97_vs_mgs_spearman.rds")

soil_ko_cosine <- readRDS("soil_ko_cosine_df.rds")
soil_ko_spearman <- readRDS("soil_ko_spearman_df.rds")
soil_ko_cosine_gg97 <- readRDS("soil_ko_picrust2_gg97_vs_mgs_cosine.rds")
soil_ko_spearman_gg97 <- readRDS("soil_ko_picrust2_gg97_vs_mgs_spearman.rds")

combined_pathabun_cosine <- readRDS("combined_pathabun_cosine_df.rds")
combined_pathabun_spearman <- readRDS("combined_pathabun_spearman_df.rds")
combined_pathcov_cosine <- readRDS("combined_pathcov_cosine_df.rds")
combined_pathcov_spearman <- readRDS("combined_pathcov_spearman_df.rds")

# Subset to KO table: NULL, PICRUSt2 and GG97 only.
hmp_ko_cosine_gg97_tab <- rbind(hmp_ko_cosine[which(hmp_ko_cosine$cat == "null"),],
                                hmp_ko_cosine[which(hmp_ko_cosine$cat == "PICRUSt2"),],
                                hmp_ko_cosine_gg97)
hmp_ko_cosine_gg97_tab$cat <- factor(hmp_ko_cosine_gg97_tab$cat)

hmp_ko_spearman_gg97_tab <- rbind(hmp_ko_spearman[which(hmp_ko_spearman$cat == "null"),],
                                hmp_ko_spearman[which(hmp_ko_spearman$cat == "PICRUSt2"),],
                                hmp_ko_spearman_gg97)
hmp_ko_spearman_gg97_tab$cat <- factor(hmp_ko_spearman_gg97_tab$cat)

mammal_ko_cosine_gg97_tab <- rbind(mammal_ko_cosine[which(mammal_ko_cosine$cat == "null"),],
                                mammal_ko_cosine[which(mammal_ko_cosine$cat == "PICRUSt2"),],
                                mammal_ko_cosine_gg97)
mammal_ko_cosine_gg97_tab$cat <- factor(mammal_ko_cosine_gg97_tab$cat)

mammal_ko_spearman_gg97_tab <- rbind(mammal_ko_spearman[which(mammal_ko_spearman$cat == "null"),],
                                  mammal_ko_spearman[which(mammal_ko_spearman$cat == "PICRUSt2"),],
                                  mammal_ko_spearman_gg97)
mammal_ko_spearman_gg97_tab$cat <- factor(mammal_ko_spearman_gg97_tab$cat)

ocean_ko_cosine_gg97_tab <- rbind(ocean_ko_cosine[which(ocean_ko_cosine$cat == "null"),],
                                ocean_ko_cosine[which(ocean_ko_cosine$cat == "PICRUSt2"),],
                                ocean_ko_cosine_gg97)
ocean_ko_cosine_gg97_tab$cat <- factor(ocean_ko_cosine_gg97_tab$cat)

ocean_ko_spearman_gg97_tab <- rbind(ocean_ko_spearman[which(ocean_ko_spearman$cat == "null"),],
                                  ocean_ko_spearman[which(ocean_ko_spearman$cat == "PICRUSt2"),],
                                  ocean_ko_spearman_gg97)
ocean_ko_spearman_gg97_tab$cat <- factor(ocean_ko_spearman_gg97_tab$cat)

soil_ko_cosine_gg97_tab <- rbind(soil_ko_cosine[which(soil_ko_cosine$cat == "null"),],
                                soil_ko_cosine[which(soil_ko_cosine$cat == "PICRUSt2"),],
                                soil_ko_cosine_gg97)
soil_ko_cosine_gg97_tab$cat <- factor(soil_ko_cosine_gg97_tab$cat)

soil_ko_spearman_gg97_tab <- rbind(soil_ko_spearman[which(soil_ko_spearman$cat == "null"),],
                                  soil_ko_spearman[which(soil_ko_spearman$cat == "PICRUSt2"),],
                                  soil_ko_spearman_gg97)
soil_ko_spearman_gg97_tab$cat <- factor(soil_ko_spearman_gg97_tab$cat)



par(mfrow=c(2,2))
boxplot(hmp_ko_cosine$metric ~ hmp_ko_cosine$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="HMP")
sinaplot(metric ~ cat, data=hmp_ko_cosine,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(hmp_ko_cosine$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(hmp_ko_cosine$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(hmp_ko_cosine$cat), xpd = TRUE)

boxplot(mammal_ko_cosine$metric ~ mammal_ko_cosine$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Mammal")
sinaplot(metric ~ cat, data=mammal_ko_cosine,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(mammal_ko_cosine$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(mammal_ko_cosine$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(mammal_ko_cosine$cat), xpd = TRUE)

boxplot(ocean_ko_cosine$metric ~ ocean_ko_cosine$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Ocean")
sinaplot(metric ~ cat, data=ocean_ko_cosine,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(ocean_ko_cosine$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(ocean_ko_cosine$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(ocean_ko_cosine$cat), xpd = TRUE)

boxplot(soil_ko_cosine$metric ~ soil_ko_cosine$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Soil")
sinaplot(metric ~ cat, data=soil_ko_cosine,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(soil_ko_cosine$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(soil_ko_cosine$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(soil_ko_cosine$cat), xpd = TRUE)

### Same but for spearman:
boxplot(hmp_ko_spearman$metric ~ hmp_ko_spearman$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="HMP")
sinaplot(metric ~ cat, data=hmp_ko_spearman,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(hmp_ko_spearman$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(hmp_ko_spearman$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(hmp_ko_spearman$cat), xpd = TRUE)

boxplot(mammal_ko_spearman$metric ~ mammal_ko_spearman$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Mammal")
sinaplot(metric ~ cat, data=mammal_ko_spearman,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(mammal_ko_spearman$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(mammal_ko_spearman$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(mammal_ko_spearman$cat), xpd = TRUE)

boxplot(ocean_ko_spearman$metric ~ ocean_ko_spearman$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Ocean")
sinaplot(metric ~ cat, data=ocean_ko_spearman,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(ocean_ko_spearman$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(ocean_ko_spearman$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(ocean_ko_spearman$cat), xpd = TRUE)

boxplot(soil_ko_spearman$metric ~ soil_ko_spearman$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Soil")
sinaplot(metric ~ cat, data=soil_ko_spearman,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(soil_ko_spearman$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(soil_ko_spearman$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(soil_ko_spearman$cat), xpd = TRUE)

par(mfrow=c(2,2))
# Plots for GG97:
boxplot(hmp_ko_cosine_gg97_tab$metric ~ hmp_ko_cosine_gg97_tab$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="HMP")
sinaplot(metric ~ cat, data=hmp_ko_cosine_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(hmp_ko_cosine_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(hmp_ko_cosine_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(hmp_ko_cosine_gg97_tab$cat), xpd = TRUE)

boxplot(mammal_ko_cosine_gg97_tab$metric ~ mammal_ko_cosine_gg97_tab$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Mammal")
sinaplot(metric ~ cat, data=mammal_ko_cosine_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(mammal_ko_cosine_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(mammal_ko_cosine_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(mammal_ko_cosine_gg97_tab$cat), xpd = TRUE)

boxplot(ocean_ko_cosine_gg97_tab$metric ~ ocean_ko_cosine_gg97_tab$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Ocean")
sinaplot(metric ~ cat, data=ocean_ko_cosine_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(ocean_ko_cosine_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(ocean_ko_cosine_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(ocean_ko_cosine_gg97_tab$cat), xpd = TRUE)

boxplot(soil_ko_cosine_gg97_tab$metric ~ soil_ko_cosine_gg97_tab$cat, outline=FALSE,
        ylab="cosine (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Soil")
sinaplot(metric ~ cat, data=soil_ko_cosine_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(soil_ko_cosine_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(soil_ko_cosine_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(soil_ko_cosine_gg97_tab$cat), xpd = TRUE)

boxplot(hmp_ko_spearman_gg97_tab$metric ~ hmp_ko_spearman_gg97_tab$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="HMP")
sinaplot(metric ~ cat, data=hmp_ko_spearman_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(hmp_ko_spearman_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(hmp_ko_spearman_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(hmp_ko_spearman_gg97_tab$cat), xpd = TRUE)

boxplot(mammal_ko_spearman_gg97_tab$metric ~ mammal_ko_spearman_gg97_tab$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Mammal")
sinaplot(metric ~ cat, data=mammal_ko_spearman_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(mammal_ko_spearman_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(mammal_ko_spearman_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(mammal_ko_spearman_gg97_tab$cat), xpd = TRUE)

boxplot(ocean_ko_spearman_gg97_tab$metric ~ ocean_ko_spearman_gg97_tab$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Ocean")
sinaplot(metric ~ cat, data=ocean_ko_spearman_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(ocean_ko_spearman_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(ocean_ko_spearman_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(ocean_ko_spearman_gg97_tab$cat), xpd = TRUE)

boxplot(soil_ko_spearman_gg97_tab$metric ~ soil_ko_spearman_gg97_tab$cat, outline=FALSE,
        ylab="spearman (16S vs MGS)", las=2, ylim=c(0, 1), col=c("tomato3"), xaxt = "n", main="Soil")
sinaplot(metric ~ cat, data=soil_ko_spearman_gg97_tab,add=TRUE, pch=19, cex=0.5, xaxt="n", yaxt="n", adjust=1, maxwidth=0.8)
# Set x-axis ticks, but no labels
axis(1, labels = FALSE, at=1:length(levels(soil_ko_spearman_gg97_tab$cat)))
# Plot x labs at default x position
text(x =  seq_along(levels(soil_ko_spearman_gg97_tab$cat)), y = par("usr")[3] - .05, srt = 45, adj = 1,
     labels = levels(soil_ko_spearman_gg97_tab$cat), xpd = TRUE)


# Plot pathway abundance validations
ggplot(combined_pathabun_spearman, aes(x=dataset, y=metric, fill=cat)) +
  geom_boxplot() +
  ylab("Spearman correlation") +
  xlab("Dataset") +
  ylim(low=0, high=1) +
  scale_fill_manual(values=c("grey", "tomato3")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(combined_pathabun_cosine, aes(x=dataset, y=metric, fill=cat)) +
  geom_boxplot() +
  ylab("cosine correlation") +
  xlab("Dataset") +
  ylim(low=0, high=1) +
  scale_fill_manual(values=c("grey", "tomato3")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Pathway coverage plot:
plot_combined_pathcov_spearman <- ggplot(combined_pathcov_spearman, aes(x=dataset, y=metric, fill=cat)) +
  geom_boxplot() +
  ylab("Spearman correlation") +
  xlab("Dataset") +
  ylim(low=0, high=1) +
  scale_fill_manual(values=c("grey", "tomato3")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# calculate p-values
wilcox.test(hmp_ko_spearman[which(hmp_ko_spearman$cat=="null"), "metric"],
            hmp_ko_spearman[which(hmp_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(hmp_ko_spearman[which(hmp_ko_spearman$cat=="PICRUSt1"), "metric"],
            hmp_ko_spearman[which(hmp_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(hmp_ko_spearman[which(hmp_ko_spearman$cat=="PanFP"), "metric"],
            hmp_ko_spearman[which(hmp_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(hmp_ko_spearman[which(hmp_ko_spearman$cat=="Piphillin"), "metric"],
            hmp_ko_spearman[which(hmp_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(hmp_ko_spearman[which(hmp_ko_spearman$cat=="Tax4Fun"), "metric"],
            hmp_ko_spearman[which(hmp_ko_spearman$cat=="PICRUSt2"), "metric"])


wilcox.test(mammal_ko_spearman[which(mammal_ko_spearman$cat=="null"), "metric"],
            mammal_ko_spearman[which(mammal_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(mammal_ko_spearman[which(mammal_ko_spearman$cat=="PICRUSt1"), "metric"],
            mammal_ko_spearman[which(mammal_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(mammal_ko_spearman[which(mammal_ko_spearman$cat=="PanFP"), "metric"],
            mammal_ko_spearman[which(mammal_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(mammal_ko_spearman[which(mammal_ko_spearman$cat=="Piphillin"), "metric"],
            mammal_ko_spearman[which(mammal_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(mammal_ko_spearman[which(mammal_ko_spearman$cat=="Tax4Fun"), "metric"],
            mammal_ko_spearman[which(mammal_ko_spearman$cat=="PICRUSt2"), "metric"])

wilcox.test(ocean_ko_spearman[which(ocean_ko_spearman$cat=="null"), "metric"],
            ocean_ko_spearman[which(ocean_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(ocean_ko_spearman[which(ocean_ko_spearman$cat=="PICRUSt1"), "metric"],
            ocean_ko_spearman[which(ocean_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(ocean_ko_spearman[which(ocean_ko_spearman$cat=="PanFP"), "metric"],
            ocean_ko_spearman[which(ocean_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(ocean_ko_spearman[which(ocean_ko_spearman$cat=="Piphillin"), "metric"],
            ocean_ko_spearman[which(ocean_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(ocean_ko_spearman[which(ocean_ko_spearman$cat=="Tax4Fun"), "metric"],
            ocean_ko_spearman[which(ocean_ko_spearman$cat=="PICRUSt2"), "metric"])

wilcox.test(soil_ko_spearman[which(soil_ko_spearman$cat=="null"), "metric"],
            soil_ko_spearman[which(soil_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(soil_ko_spearman[which(soil_ko_spearman$cat=="PICRUSt1"), "metric"],
            soil_ko_spearman[which(soil_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(soil_ko_spearman[which(soil_ko_spearman$cat=="PanFP"), "metric"],
            soil_ko_spearman[which(soil_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(soil_ko_spearman[which(soil_ko_spearman$cat=="Piphillin"), "metric"],
            soil_ko_spearman[which(soil_ko_spearman$cat=="PICRUSt2"), "metric"])
wilcox.test(soil_ko_spearman[which(soil_ko_spearman$cat=="Tax4Fun"), "metric"],
            soil_ko_spearman[which(soil_ko_spearman$cat=="PICRUSt2"), "metric"])