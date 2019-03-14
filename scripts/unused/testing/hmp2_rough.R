





### TESTING ###


# CLR transform and filter out features found in less than 20% of samples.
clr_and_filt <- function(in_df) {
  in_df <- in_df + 1
  in_df_clr <- data.frame(apply(in_df, 2, function(x){log(x) - mean(log(x))}))
  low_freq_rows <- which(rowSums(in_df > 1) < 0.2*ncol(in_df))
  if(length(low_freq_rows)) {
    return(in_df_clr[-low_freq_rows , ])
  }
}


# Define function to read in dataframe, scramble last column, and run RF.
run_scrambled_RF <- function(in_df, num_tree=501) {
  in_df[,ncol(in_df)] <- sample(in_df[,ncol(in_df)])
  in_df_rf <- return(randomForest(x=in_df[,1:(ncol(in_df)-1)],
                                  y=in_df[ , ncol(in_df)],
                                  ntree=501))
}

parse_varImp <- function(rf_output) {
  
  rf_imp <- as.data.frame(rf_output$importance)
  
  rf_imp$features <- rownames(rf_imp)
  
  rf_imp_sorted <- arrange(rf_imp, desc(MeanDecreaseAccuracy))
  
  rownames(rf_imp_sorted) <- rf_imp_sorted$features
  
  rf_imp_sorted <- rf_imp_sorted[, -which(colnames(rf_imp_sorted) %in% "features")]
  
  return(rf_imp_sorted)
}


run_quick_rf <- function(in_df, mapfile, num_tree=1000) {
  
  in_df_prep <- data.frame(t(in_df[,rownames(mapfile) ]))
  in_df_prep$diagnosis <- as.factor(mapfile$diagnosis)
  
  RF_in_df_prep <- randomForest( x=in_df_prep[,1:(ncol(in_df_prep)-1)],
                                 y=in_df_prep[ , ncol(in_df_prep)],
                                 ntree=num_tree, importance=TRUE, proximity=TRUE)
  
}

# Get mean % contribution by a particular clade. Returns list for each function with % contributed by that clade.
df_clade_mean_contrib <- function(strat_func_by_tax_rel, clade_string, num_cores) {
  
  funcs <- unique(gsub("\\|.+$", "", rownames(strat_func_by_tax_rel)))
  
  func_mean_clade_contrib <- mclapply(funcs, function(x) { 
    single_func_mean_contrib(x, strat_func_by_tax_rel, clade_string)
    
  }, mc.cores=num_cores)
  
  names(func_mean_clade_contrib) <- funcs
  
  return(func_mean_clade_contrib)
  
}

single_func_mean_contrib <- function(func, func_df, clade_string) {
  
  func_rows <- func_df[grep(func, rownames(func_df)),]
  
  if(length(grep(clade_string, rownames(func_rows))) > 0) {
    
    clade_func_abun <- as.numeric(func_rows[grep(clade_string, rownames(func_rows)),])
    func_abun <- colSums(func_rows) + 0.00001
    
    return(mean(clade_func_abun / func_abun) * 100)
    
  } else {
    
    return(0)
    
  }
  
}

# Read in mapfile.
hmp2_map <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_16S_reformat.tsv", sep="\t", header=T, 
                       comment.char="", quote="", stringsAsFactors = FALSE)

# Add "HMP.2." to the start of all external ids.
hmp2_map$External.ID <- paste(rep_len("HMP.2.", length(hmp2_map$External.ID)), hmp2_map$External.ID, sep="")
hmp2_map <- hmp2_map[-which(duplicated(hmp2_map$External.ID)),]
rownames(hmp2_map) <- hmp2_map$External.ID

hmp2_map_full <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_col_subset.txt",
                            sep="\t", header=T, comment.char="", quote="", stringsAsFactors = FALSE)

hmp2_map_metabolome <- hmp2_map_full[which(hmp2_map_full$data_type == "metabolomics"),]
rownames(hmp2_map_metabolome) <- hmp2_map_metabolome$External.ID

hmp2_map_reformat <- read.table("/home/gavin/gavin_backup/projects/hmp2_ibd_product/full_metadata/hmp2_metadata_rnaseq_reformat.tsv",
                                sep="\t", header=T, comment.char="", quote="", stringsAsFactors = FALSE)
hmp2_map_reformat$matching_16S_id <- paste("HMP.2.", hmp2_map_reformat$matching_16S_id, sep="")
hmp2_map_reformat <- hmp2_map_reformat[-which(hmp2_map_reformat$matching_16S_id == "HMP.2.NA"),]
hmp2_map_reformat <- hmp2_map_reformat[-which(duplicated(hmp2_map_reformat$matching_16S_id)),]
rownames(hmp2_map_reformat) <- hmp2_map_reformat$External.ID


metabolome_samples <- colnames(metabolome)[which(colnames(metabolome) %in% rownames(hmp2_map_metabolome))]
hmp2_map_metabolome_subset <- hmp2_map_metabolome[metabolome_samples,]


hmp2_map_full[colnames(metabolome) , "matching_16S_id"]



# Read in KO, pathway abundances, mapfile, and ASV abundances.
ko_in <- read.table("16S/April2018_redownload/picrust2_full_output_2.1.0-b/KO_metagenome_out/pred_metagenome_unstrat.tsv",
                    header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", quote="", row.names=1)
ko_in_strat <- read.table("16S/April2018_redownload/picrust2_full_output_2.1.0-b/KO_metagenome_out/pred_metagenome_strat.tsv",
                          header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", quote="")
rownames(ko_in_strat) <- paste(ko_in_strat$function., ko_in_strat$sequence, sep="|")

pathabun_in <- read.table("16S/April2018_redownload/picrust2_full_output_2.1.0-b/pathways_out/path_abun_unstrat.tsv",
                          header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", quote="", row.names=1)
pathabun_in_strat <- read.table("16S/April2018_redownload/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                header=T, sep="\t", stringsAsFactors = FALSE, comment.char="", quote="")
rownames(pathabun_in_strat) <- paste(pathabun_in_strat$pathway, pathabun_in_strat$sequence, sep="|")

metabolome <- read.csv("/home/gavin/projects/hmp2_ibd_working/HMP2_metabolomics_proteomics/HMP2_metabolomics.csv", header=T, sep=",", comment.char="", stringsAsFactors = FALSE)
metabolome_raw <- metabolome
# metabolome <- metabolome_raw
rownames(metabolome) <- metabolome$Compound
metabolome <- data.frame(t(metabolome[, metabolome_samples]), check.names=FALSE)
metabolome[is.na(metabolome)] <- 0

metabolome$participant <- hmp2_map_metabolome[metabolome_samples, "Participant.ID"]

metabolome_mean <- aggregate(. ~ participant, data=metabolome, FUN=mean)
metabolome_mean_raw <- metabolome_mean
rownames(metabolome_mean) <- metabolome_mean$participant
metabolome_mean <- metabolome_mean[, -which(colnames(metabolome_mean) == "participant")]

# Remove unnecessary columns.
col2remove <- c("function", "pathway", "description", "sequence")
ko_in_strat <- ko_in_strat[,-which(colnames(ko_in_strat) %in% col2remove)]
pathabun_in_strat <- pathabun_in_strat[,-which(colnames(pathabun_in_strat) %in% col2remove)]

# Split into samples overlapping with metabolome for ileum and rectum separately.
map_Rectum <- hmp2_map[which(hmp2_map$biopsy_location == "Rectum"), ]
map_Ileum <- hmp2_map[which(hmp2_map$biopsy_location == "Ileum"), ]

map_Rectum_nodup <- map_Rectum[-which(duplicated(map_Rectum$Participant.ID)),]
map_Ileum_nodup <- map_Ileum[-which(duplicated(map_Ileum$Participant.ID)),]

named_compounds <- metabolome_raw[which(metabolome_raw$Metabolite != ""), "Compound"]

pathabun_in_patient_Rectum <- data.frame(t(pathabun_in))
pathabun_in_patient_Rectum <- pathabun_in_patient_Rectum[-which(! rownames(pathabun_in_patient_Rectum) %in% rownames(map_Rectum_nodup)),]
rownames(pathabun_in_patient_Rectum) <- map_Rectum_nodup[rownames(pathabun_in_patient_Rectum), "Participant.ID"]

pathabun_in_patient_Ileum <- data.frame(t(pathabun_in))
pathabun_in_patient_Ileum <- pathabun_in_patient_Ileum[-which(! rownames(pathabun_in_patient_Ileum) %in% rownames(map_Ileum_nodup)),]
rownames(pathabun_in_patient_Ileum) <- map_Ileum_nodup[rownames(pathabun_in_patient_Ileum), "Participant.ID"]

metabolome_mean_Rectum_overlap <- metabolome_mean[rownames(pathabun_in_patient_Rectum),]
metabolome_mean_Ileum_overlap <- metabolome_mean[rownames(pathabun_in_patient_Ileum),]

metabolome_mean_Rectum_overlap_relab <- data.frame(sweep(metabolome_mean_Rectum_overlap, 1, rowSums(metabolome_mean_Rectum_overlap), '/') * 100)
metabolome_mean_Rectum_overlap_relab_filt <- metabolome_mean_Rectum_overlap_relab[, named_compounds]
metabolome_mean_Rectum_overlap_relab_filt <- metabolome_mean_Rectum_overlap_relab_filt[, -which(colSums(metabolome_mean_Rectum_overlap_relab_filt > 0) < 22)]

pathabun_in_patient_Rectum_relab <- data.frame(sweep(pathabun_in_patient_Rectum, 1, rowSums(pathabun_in_patient_Rectum), '/') * 100)
pathabun_in_patient_Rectum_relab_filt <- pathabun_in_patient_Rectum_relab[, -which(colSums(pathabun_in_patient_Rectum_relab > 0) < 22)]

metabolome_mean_Ileum_overlap_relab <- data.frame(sweep(metabolome_mean_Ileum_overlap, 1, rowSums(metabolome_mean_Ileum_overlap), '/') * 100)
metabolome_mean_Ileum_overlap_relab_filt <- metabolome_mean_Ileum_overlap_relab[, named_compounds]
metabolome_mean_Ileum_overlap_relab_filt <- metabolome_mean_Ileum_overlap_relab_filt[, -which(colSums(metabolome_mean_Ileum_overlap_relab_filt > 0) < 22)]

pathabun_in_patient_Ileum_relab <- data.frame(sweep(pathabun_in_patient_Ileum, 1, rowSums(pathabun_in_patient_Ileum), '/') * 100)
pathabun_in_patient_Ileum_relab_filt <- pathabun_in_patient_Ileum_relab[, -which(colSums(pathabun_in_patient_Ileum_relab > 0) < 22)]


# Get correlations between named compounds and MetaCyc pathways relative abundances in the rectum.
metabolome_vs_pathabun_Rectum_df <- data.frame(matrix(NA, nrow=ncol(metabolome_mean_Rectum_overlap_relab_filt)*ncol(pathabun_in_patient_Rectum_relab_filt), ncol=4))
colnames(metabolome_vs_pathabun_Rectum_df) <- c("metabolite", "pathway", "rho", "p")


i = 0
for(metabolite in colnames(metabolome_mean_Rectum_overlap_relab_filt)) {
  for(pathway in colnames(pathabun_in_patient_Rectum_relab_filt)) {
    i = i + 1
    
    cor_out <- cor.test(as.numeric(metabolome_mean_Rectum_overlap_relab_filt[, metabolite]), as.numeric(pathabun_in_patient_Rectum_relab_filt[, pathway]), method="spearman")
    metabolome_vs_pathabun_Rectum_df[i, c("metabolite", "pathway")] <- c(metabolite, pathway) 
    metabolome_vs_pathabun_Rectum_df[i, c("rho", "p")] <- c(cor_out$estimate, cor_out$p.value)
  }
}

metabolome_vs_pathabun_Ileum_df$fdr <- p.adjust(metabolome_vs_pathabun_Ileum_df$p, "fdr")

metabolome_vs_pathabun_Ileum_df <- data.frame(matrix(NA, nrow=ncol(metabolome_mean_Ileum_overlap_relab_filt)*ncol(pathabun_in_patient_Ileum_relab_filt), ncol=4))
colnames(metabolome_vs_pathabun_Ileum_df) <- c("metabolite", "pathway", "rho", "p")


i = 0
for(metabolite in colnames(metabolome_mean_Ileum_overlap_relab_filt)) {
  for(pathway in colnames(pathabun_in_patient_Ileum_relab_filt)) {
    i = i + 1
    
    cor_out <- cor.test(as.numeric(metabolome_mean_Ileum_overlap_relab_filt[, metabolite]), as.numeric(pathabun_in_patient_Ileum_relab_filt[, pathway]), method="spearman")
    metabolome_vs_pathabun_Ileum_df[i, c("metabolite", "pathway")] <- c(metabolite, pathway) 
    metabolome_vs_pathabun_Ileum_df[i, c("rho", "p")] <- c(cor_out$estimate, cor_out$p.value)
  }
}

metabolome_vs_pathabun_Ileum_df$fdr <- p.adjust(metabolome_vs_pathabun_Ileum_df$p, "fdr")
metabolome_vs_pathabun_Ileum_df_sig <- metabolome_vs_pathabun_Ileum_df[which(metabolome_vs_pathabun_Ileum_df$fdr < 0.05),]


compound2metabolite <- metabolome_raw[, c("Compound", "Metabolite")]
rownames(compound2metabolite) <- compound2metabolite$Compound

metabolome_vs_pathabun_Ileum_df_sig$metabolite_name <- compound2metabolite[metabolome_vs_pathabun_Ileum_df_sig$metabolite, "Metabolite"]



plot(as.numeric(metabolome_mean_Ileum_overlap_relab_filt[, "C18n_QI73"]),
     as.numeric(pathabun_in_patient_Ileum_relab_filt[, "ANAEROFRUCAT.PWY"]))


tmp_input <- data.frame(pathway=hmp2_pathabun_ileum_CD_filt[, "ANAEROFRUCAT-PWY"],
                        gene=host_rna_rpm_ileum_CD_subset_log[, "APOA1"],
                        site=rnaseq_meta_ileum_CD$site_name,
                        consent_age=rnaseq_meta_ileum_CD$consent_age)

pcor.test(tmp_input$pathway, tmp_input$gene, tmp_input$gene)

pcor.test(x=tmp_input$pathway,
          y=tmp_input$gene,
          z=tmp_input[, c("site", "consent_age")],
          method = "spearman")


tmp3 <- anova(tmp1, tmp2)

overlapping_rectum <- rownames(host_rna_rpm_rectum)[which(rownames(host_rna_rpm_rectum) %in% colnames(hmp2_pathabun))]
hmp2_pathabun_rectum <- data.frame(t(hmp2_pathabun[, overlapping_rectum]))
host_rna_rpm_rectum_subset <- host_rna_rpm_rectum[overlapping_rectum,]



# partial correlation between "hl" and "disp" given "deg" and "BC"
pcor.test(y.data$hl,y.data$disp,y.data[,c("deg","BC")])
pcor.test(y.data[,1],y.data[,2],y.data[,c(3:4)])
pcor.test(y.data[,1],y.data[,2],y.data[,-c(1:2)])



p_values <- c()
for(pathway in colnames(hmp2_pathabun_ileum_CD_filt)) {
  print(pathway)
  p_values <- c(p_values, cor.test(as.numeric(hmp2_pathabun_ileum_CD_filt[, pathway]),
                                   as.numeric(host_rna_rpm_ileum_subset[CD_ileum_samples, "APOA1"]),
                                   method="spearman")$p.value)
}

p_values_fdr <- p.adjust(p_values, "fdr")

plot(hmp2_pathabun_ileum_CD_filt$P562.PWY, host_rna_rpm_ileum_subset[CD_ileum_samples, "APOA1"],
     xlab="P562-PWY - myo-inositol degradation I",
     ylab="APOA1 reads per million", pch=16)

hmp2_pathabun_strat <- read.table("/home/gavin/projects/hmp2_ibd_working/16S/April2018_redownload/picrust2_full_output_2.1.0-b/pathways_out/path_abun_strat.tsv",
                                  header=T, sep="\t")

hmp2_pathabun_strat[, 3:ncol(hmp2_pathabun_strat)] <- data.frame(sweep(hmp2_pathabun_strat[, 3:ncol(hmp2_pathabun_strat)], 2,
                                                                       colSums(hmp2_pathabun_strat[, 3:ncol(hmp2_pathabun_strat)], na.rm = TRUE), '/') * 100)

hmp2_pathabun_strat_P562 <- hmp2_pathabun_strat[which(hmp2_pathabun_strat$pathway == "P562-PWY"),]

rownames(hmp2_pathabun_strat_P562) <- hmp2_pathabun_strat_P562$sequence

hmp2_pathabun_strat_P562 <- hmp2_pathabun_strat_P562[, -c(1, 2)]

colnames(hmp2_pathabun_strat_P562) <- gsub("HMP\\.2\\.", "", colnames(hmp2_pathabun_strat_P562))

hmp2_pathabun_strat_P562 <- hmp2_pathabun_strat_P562[, matching_16S_ids]
colnames(hmp2_pathabun_strat_P562) <- rnaseq_meta_subset[matching_16S_ids, "External.ID"]

hmp2_pathabun_strat_P562_ileum_CD <- hmp2_pathabun_strat_P562[, CD_ileum_samples]
hmp2_pathabun_strat_P562_ileum_CD <- hmp2_pathabun_strat_P562_ileum_CD[-which(rowSums(is.na(hmp2_pathabun_strat_P562_ileum_CD)) > 0),]

library(reshape2)
library(ggplot2)

hmp2_pathabun_strat_P562_ileum_CD_tmp <- hmp2_pathabun_strat_P562_ileum_CD
hmp2_taxa <- read.table("/home/gavin/projects/hmp2_ibd_working/16S/April2018_redownload/deblur_output_final/hmp2_ibd_16S_taxonomy_unfilt.tsv", header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
hmp2_taxa_levels <- add_tax_cols(hmp2_taxa)

hmp2_pathabun_strat_P562_ileum_CD_tmp$genus <- hmp2_taxa_levels[rownames(hmp2_pathabun_strat_P562_ileum_CD_tmp), "Genus"]

hmp2_pathabun_strat_P562_ileum_CD_genus_sum <- aggregate(. ~ genus, data=hmp2_pathabun_strat_P562_ileum_CD_tmp, FUN=sum)

hmp2_pathabun_strat_P562_ileum_CD_genus_sum <- hmp2_pathabun_strat_P562_ileum_CD_genus_sum[-which(rowSums(hmp2_pathabun_strat_P562_ileum_CD_genus_sum[, -1]) == 0),]

tmp_melt <- melt(hmp2_pathabun_strat_P562_ileum_CD_genus_sum)
tmp_melt$variable <- factor(tmp_melt$variable, levels=names(sort(colSums(hmp2_pathabun_strat_P562_ileum_CD))))

# Select colours based on same phylum.
actino_col <- c("indianred", "indianred1", "lightpink2", "firebrick1")
firmicutes_col <- c("skyblue1", "skyblue3", "slateblue", "slateblue4", "turquoise1", "turquoise3")
proteo_col <- c("seagreen", "seagreen1", "springgreen", "greenyellow", "olivedrab3")

col2use <- c(actino_col, firmicutes_col, proteo_col)

ggplot(tmp_melt, aes(x=variable, y=value)) +
  geom_bar(aes(fill = genus), stat = "identity") +
  ylab("P562-PWY - myo-inositol degradation I") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values=col2use)


# Test for correlations with 3 significant genes with Phyla relab.

hmp2_biom <- readRDS("prepped_tables/hmp2_biom_Ileum.rds")
asv_tax_in <- read.table("hmp2_asv_taxa.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
asv_tax_in_levels <- add_tax_cols(asv_tax_in)
rownames(asv_tax_in_levels) <- asv_tax_in_levels$Feature.ID

hmp2_biom$Phylum <- asv_tax_in_levels[rownames(hmp2_biom), "Phylum"]

hmp2_biom_Phylum <- aggregate(. ~ Phylum, FUN=sum, data=hmp2_biom)

rownames(hmp2_biom_Phylum) <- hmp2_biom_Phylum$Phylum
hmp2_biom_Phylum <- hmp2_biom_Phylum[, -which(colnames(hmp2_biom_Phylum) == "Phylum")]


host_rna_rpm_ileum_CD_subset <- host_rna_rpm_ileum_CD[, c("DUOX2", "MMP3", "NAT8")]

hmp2_biom_Phylum_t <- data.frame(t(hmp2_biom_Phylum), check.names = FALSE)
hmp2_biom_Phylum_t <- hmp2_biom_Phylum_t[rownames(host_rna_rpm_ileum_CD_subset),]
hmp2_biom_Phylum_t <- hmp2_biom_Phylum_t[, -which(colSums(hmp2_biom_Phylum_t) == 0)]

Phylum_vs_gene_spearman <- data.frame(matrix(NA, ncol=4, nrow=ncol(hmp2_biom_Phylum_t) * ncol(host_rna_rpm_ileum_CD_subset)))

colnames(Phylum_vs_gene_spearman) <- c("Phylum", "Gene", "rho", "p")

i = 1

for(Phylum in colnames(hmp2_biom_Phylum_t)) {
  
  for(gene in colnames(host_rna_rpm_ileum_CD_subset)) {
    
    Phylum_vs_gene_spearman[i, c("Phylum", "Gene")] <- c(Phylum, gene)
    
    cor_out <- cor.test(hmp2_biom_Phylum_t[, Phylum], host_rna_rpm_ileum_CD_subset[, gene], method="spearman")
    
    Phylum_vs_gene_spearman[i, c("rho", "p")] <- c(cor_out$estimate, cor_out$p.value)
    
    i = i + 1
    
  }
  
}


