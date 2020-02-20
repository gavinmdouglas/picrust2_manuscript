library(sinaplot)
library(parallel)
library(stringi)

asinTransform <- function(p) { 2 * asin(sqrt(p)) }

compare_func_ratios_by_asv_groups <- function(strat_table, ASVs_of_interest, sample_group1, sample_group2) {
  
  # Set first column name to be "func" (even if it's pathway) and restrict to samples to those of interest.
  orig_colnames <- colnames(strat_table)
  colnames(strat_table) <- c("func", orig_colnames[2:length(orig_colnames)])
  strat_table <- strat_table[, c("func", "sequence", sample_group1, sample_group2)]
  
  # Remove rows that are all 0s.
  zero_rows <- which(rowSums(strat_table[, 3:ncol(strat_table)]) == 0)
  if(length(zero_rows) > 0) {
    strat_table <- strat_table[-zero_rows,]
  }
  
  # Split table by ASV groups.
  strat_table_asv_interest <- strat_table[which(strat_table$sequence %in% ASVs_of_interest),]
  strat_table_asv_other <- strat_table[which(! strat_table$sequence %in% ASVs_of_interest),]
  
  # Only keep functions that are found at least once in the ASVs of interest.
  func2keep <- strat_table_asv_interest$func
  
  if(length(which(duplicated(func2keep))) > 0) {
    func2keep <- func2keep[-which(duplicated(func2keep))]
  }
  
  strat_table_asv_interest <- strat_table_asv_interest[which(strat_table_asv_interest$func %in% func2keep), ]
  strat_table_asv_other <- strat_table_asv_other[which(strat_table_asv_other$func %in% func2keep), ]
  
  # Remove sequence column.
  strat_table_asv_interest <- strat_table_asv_interest[, -which(colnames(strat_table_asv_interest) == "sequence")]
  strat_table_asv_other <- strat_table_asv_other[, -which(colnames(strat_table_asv_other) == "sequence")]
  
  # Get summed function counts by function.
  strat_table_asv_interest_summed <- aggregate(. ~ func, data=strat_table_asv_interest, FUN=sum)
  strat_table_asv_other_summed <- aggregate(. ~ func, data=strat_table_asv_other, FUN=sum)
  
  # Set functions as rownames and remove column.
  rownames(strat_table_asv_interest_summed) <- strat_table_asv_interest_summed$func
  rownames(strat_table_asv_other_summed) <- strat_table_asv_other_summed$func
  strat_table_asv_interest_summed <- strat_table_asv_interest_summed[, -1]
  strat_table_asv_other_summed <- strat_table_asv_other_summed[, -1]
  
  # Add in pathways that might be missing from each table.
  all_funcs <- c(rownames(strat_table_asv_interest_summed), rownames(strat_table_asv_other_summed))
  
  if(length(which(duplicated(all_funcs))) > 0) {
    all_funcs <- all_funcs[-which(duplicated(all_funcs))]
  }
  
  strat_table_asv_interest_summed_NoMiss <- add_missing_funcs(strat_table_asv_interest_summed, all_funcs)
  strat_table_asv_other_summed_NoMiss <- add_missing_funcs(strat_table_asv_other_summed, all_funcs)  
  
  # Sort dataframes to be in same order and add pseudocount.
  strat_table_asv_interest_summed_NoMiss <- strat_table_asv_interest_summed_NoMiss + 1
  strat_table_asv_other_summed_NoMiss <- strat_table_asv_other_summed_NoMiss[rownames(strat_table_asv_interest_summed_NoMiss),] + 1
  
  # Get table of ratios between the ASV groups.
  out_ratio <- strat_table_asv_interest_summed_NoMiss/strat_table_asv_other_summed_NoMiss
  
  # Test for different ratios between each funcs between the sample groups.
  ratio_wilcoxon_raw_p <- c()
  ratio_wilcoxon_raw_w <- c()
  path_mean_diff <- c()
  
  for(func in rownames(out_ratio)) {
    
    func_wilcox <- wilcox.test(as.numeric(out_ratio[func, sample_group1]),
                               as.numeric(out_ratio[func, sample_group2]))
    
    ratio_wilcoxon_raw_p <- c(ratio_wilcoxon_raw_p, func_wilcox$p.value)
    ratio_wilcoxon_raw_w <- c(ratio_wilcoxon_raw_w, func_wilcox$statistic)
    
    path_mean_diff <- c(path_mean_diff, mean(as.numeric(out_ratio[func, sample_group1])) - mean(as.numeric(out_ratio[func, sample_group2])))
    
  }
  
  ratio_wilcoxon_raw <- data.frame(feature=rownames(out_ratio),
                                   wilcox_p=ratio_wilcoxon_raw_p,
                                   wilcox_w=ratio_wilcoxon_raw_w,
                                   mean_diff=path_mean_diff)

  return(list(interest=strat_table_asv_interest_summed_NoMiss,
              other=strat_table_asv_other_summed_NoMiss,
              ratio=out_ratio,
              wilcox=ratio_wilcoxon_raw))
  
}


# Function to subset mapfile by specified column and value and will take the first row for any
# "Participant.ID" ids that are in multiple rows (replicates).
subset_hmp_map <- function(mapfile, col_interest, col_value) {
  mapfile <- mapfile[which(mapfile[, col_interest] == col_value), , drop=FALSE]
  
  dup_samples <- mapfile[which(duplicated(mapfile$Participant.ID)), "Participant.ID"]
  
  if(length(dup_samples) > 0) {
    print("These are the duplicate participant rows in mapfile:")
    print(mapfile[which(mapfile$Participant.ID %in% dup_samples), ])
    
    mapfile <- mapfile[-which(duplicated(mapfile$Participant.ID)) , ]
    
  }
  
  return(mapfile)
}


plot_rpm_boxplot_and_sinaplot <- function(rpm_table, gene_col, factor_col, tissue=NULL, deseq2_results=NULL, ylab_line=4) {

  par(mar=c(5.1, 5.1, 4.1, 2.1))
  
  in_table <- substitute(rpm_table)
  
  # Set main title to be gene column name.
  main_title <- gene_col

  # Add tissue and deseq2 results P-value as well if specified.
  
  if(! is.null(tissue)) { 
    main_title <- paste(main_title, tissue)
  }

  if(! is.null(deseq2_results)) {
    adj_p <- deseq2_results[gene_col, "padj"]
    
    # Add asterix if adjusted P-value is less than 0.05.
    if(adj_p < 0.05) {
      main_title <- paste(main_title, "*, P=", formatC(adj_p, format = "e", digits = 2), sep="")
    } else {
      main_title <- paste(main_title, ", P=", formatC(adj_p, format = "e", digits = 2), sep="")
    }
  }
  
  gene_col <- as.name(gene_col)
  
  eval(substitute(boxplot(gene_col~factor_col,
                          data=in_table,
                          outline=FALSE,
                          ylim=c(0, (max(rpm_table$gene_col) + max(rpm_table$gene_col)*0.2)),
                          ylab="",
                          main=main_title,
                          las=1)))
  
  title(ylab="Reads per million", line=ylab_line, cex.lab=1)
  
  eval(substitute(sinaplot(gene_col~factor_col,
                           data=in_table,
                           add=TRUE,
                           yaxt='n')))
  
  par(mar=c(5.1,4.1,4.1,2.1))
}


read_biom_tsv <- function(filename) {
 
  return(read.table(filename,
                    header = T,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    comment.char = "",
                    skip = 1,
                    row.names=1))
}


# Function to split up single Taxon column output by QIIME2 into separate columns.
# Note that "unclassified" values, e.g. "s__", will be added even in cases where
# the classifier didn't reach lower levels.
add_tax_cols <- function(taxa_table) {
  
  # Make empty columns for taxa levels.
  taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  for(taxa_level in taxa_levels) { 
    taxa_table[, taxa_level] <- NA 
  }

  unclassified_labels <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  
  # Loop over all ASVs and parse the taxa information into these individual columns.
  for(i in 1:nrow(taxa_table)) {
    
    # Get that ASV's full taxonomy.
    # This string will be altered at each loop of a taxonomic level below.
    seq_taxonomy <- taxa_table[i, "Taxon"]
    
    # Some taxonomy assignments wont go down to species so need to skip 
    # those by checking how many levels are classified (will do this based on length below).
    exp_length <- 7
    
    for(tax_level in rev(taxa_levels)) {
      
      # Go to next iteration if the lowest taxonomic level is higher than current tax_level.
      exp_length_diff <- exp_length - length(stri_split(seq_taxonomy, regex = "; ")[[1]])
      if(exp_length_diff != 0) {
        taxa_table[i, tax_level] <- paste(c(seq_taxonomy, 
                                           unclassified_labels[(exp_length - exp_length_diff + 1):exp_length]), 
                                         collapse="; ")
        exp_length <- exp_length - 1
        next
      }
      
      # Set taxa level for this sequence.
      taxa_table[i, tax_level] <- seq_taxonomy
      
      # Remove the last level from the taxonomy.
      seq_taxonomy <- sub("; [^;]*$", "", seq_taxonomy)
      
      exp_length <- exp_length - 1
    }
  }
  
  return(taxa_table)
}


# Function that reads in table outputted by add_tax_cols and BIOM Table and 
# calculates abundance tables at each taxonomic level. Will return named list.
abun_by_taxa_level <- function(taxa_table_w_levels, in_biom) {
  
  taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  summed_tables <- list()
  
  for(taxa_level in taxa_levels) {
  
    in_biom_tax <- in_biom
    in_biom_tax[, taxa_level] <- factor(taxa_table_w_levels[rownames(in_biom), taxa_level])
    
    summed_tab <- aggregate(x=in_biom[, 1:(ncol(in_biom_tax)-1)], by=list(in_biom_tax[, ncol(in_biom_tax)]), sum)
    
    # Set rownames to be first column and reverse this column before adding to list.
    rownames(summed_tab) <- summed_tab[, 1]
    summed_tables[[taxa_level]] <- summed_tab[, -1]
  }
  
  return(summed_tables)
}


cor_df_cols <- function(df_in, col_set1, col_set2) {
  
  out_df <- data.frame(matrix(NA, nrow=length(col_set1)*length(col_set2), ncol=4))
  
  colnames(out_df) <- c("feature1", "feature2", "estimate", "p.value")
  
  out_i <- 1
  
  for(col1 in col_set1) {
    
    for(col2 in col_set2) {
      
      cor_test_out <- cor.test(df_in[,col1],
                               df_in[,col2],
                               method="spearman",
                               exact=FALSE)
      
      out_df[out_i,] <- c(col1, col2, cor_test_out$estimate, cor_test_out$p.value)
      
      out_i <- out_i + 1 
    }
    
  }
  
  return(out_df)
  
}


# Function to return the odds ratios for a set number of permutations.
permute_func_dist <- function(func_tab, subset_size, num_permutations=1000, num_cores=1) {
  
  odds_ratios_list <- mclapply(c(1:num_permutations), function(x) {
    
    row_subset <- sample(1:nrow(func_tab), size = subset_size, replace = FALSE)
    
    func_subset <- func_tab[row_subset,]
    
    func_other <- func_tab[-row_subset,]
    
    return(((colSums(func_subset) + 1)/sum(func_subset))/((colSums(func_other) + 1)/sum(func_other)))
  }, mc.cores=num_cores)
  
  odds_ratios <- data.frame(matrix(unlist(odds_ratios_list), nrow=num_permutations, byrow=TRUE))
  
  colnames(odds_ratios) <- colnames(func_tab)
  
  return(odds_ratios)
}


# Function to return the odds ratio of the abundance of a function in a subset of genomes compared
# to the rest of genomes.
func_odds_ratios <- function(func_tab, row_subset) {
  
  row_subset_i <- which(rownames(func_tab) %in% row_subset)
  
  func_subset <- func_tab[row_subset_i,]
  
  func_other <- func_tab[-row_subset_i,]
  
  odds_ratios <- ((colSums(func_subset) + 1)/sum(func_subset))/((colSums(func_other) + 1)/sum(func_other))
  
  return(odds_ratios)
}

# Get proportion of null values more extreme than observed
enriched_funcs_p <- function(null_OR, observed_OR) {
  
  enrich_p <- unlist(lapply(names(observed_OR), function(x) {
                     length(which(null_OR[, x] >= observed_OR[x]))/nrow(null_OR)
  }))
  
  names(enrich_p) <- names(observed_OR)
  
  return(enrich_p)
}

create_procrustes_bioplot <- function(procrustes_out, category1, category2) {
  library(ggplot2)
  library(grid)
  
  data_coor <- data.frame(rda1=procrustes_out$Yrot[,1],
                          rda2=procrustes_out$Yrot[,2],
                          xrda1=procrustes_out$X[,1],
                          xrda2=procrustes_out$X[,2])
  
  return(ggplot(data_coor) +
           geom_point(aes(x=rda1, y=rda2, colour=category1),  size=4) +
           geom_point(aes(x=xrda1, y=xrda2, colour=category2),  size=4) +
           geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2),arrow=arrow(length=unit(0.2,"cm"))) +
           xlab("Dimension 1") + ylab("Dimension 2")+ theme_bw() + 
           theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                 legend.title = element_blank(), text = element_text(size=20)))
}
