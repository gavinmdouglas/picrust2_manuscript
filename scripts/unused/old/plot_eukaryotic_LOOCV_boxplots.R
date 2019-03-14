setwd("/home/gavin/projects/picrust2_manuscript/data/taxa_LOOCV/")

source("/home/gavin/projects/picrust2_manuscript/analyses/picrust2_ms_functions.R")

# Calculate null distribution for each E.C. database.
ec_18S <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_18S.txt",
                     header=T, row.names=1, stringsAsFactors = FALSE, sep="\t", check.names=FALSE)

ec_ITS <- read.table("/home/gavin/projects/picrust_pipeline/RefSeq_18S_ITS/mean_func_tables/ec_ITS.txt",
                     header=T, row.names=1, stringsAsFactors = FALSE, sep="\t", check.names=FALSE)

ec_18S_rand <- ec_18S
rownames(ec_18S_rand) <- sample(rownames(ec_18S))
ec_ITS_rand <- ec_ITS
rownames(ec_ITS_rand) <- sample(rownames(ec_ITS))

tmp1 <- data.frame(t(ec_18S_rand))
tmp2 <- data.frame(t(ec_18S))



null_18S_cosine <- cor_all_cols(tab1 = t(ec_18S_rand), tab2 = t(ec_18S), cat_string="Null", metric="cosine")
null_18S_spearman <- cor_all_cols(tab1 = t(ec_18S_rand), tab2 = t(ec_18S), cat_string="Null", metric="spearman")
null_ITS_cosine <- cor_all_cols(tab1 = t(ec_ITS_rand), tab2 = t(ec_ITS), cat_string="Null", metric="cosine")
null_ITS_spearman <- cor_all_cols(tab1 = t(ec_ITS_rand), tab2 = t(ec_ITS), cat_string="Null", metric="spearman")

ec_class_18S <- read.table("18S_Class_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_order_18S <- read.table("18S_Order_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_family_18S <- read.table("18S_Family_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_genus_18S <- read.table("18S_Genus_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_species_18S <- read.table("18S_Species_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_assembly_18S <- read.table("18S_assembly_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)

ec_class_ITS <- read.table("ITS_Class_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_order_ITS <- read.table("ITS_Order_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_family_ITS <- read.table("ITS_Family_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_genus_ITS <- read.table("ITS_Genus_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_species_ITS <- read.table("ITS_Species_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)
ec_assembly_ITS <- read.table("ITS_assembly_metrics.tsv", header=T, sep="\t", stringsAsFactors = FALSE)


par(mfrow=c(1,2))
boxplot(null_18S_spearman$metric, ec_class_18S$spearman, ec_order_18S$spearman, ec_family_18S$spearman, ec_genus_18S$spearman, ec_species_18S$spearman, ec_assembly_18S$spearman,
        main="18S Spearman", ylab="Spearman correlation", names=c("Null", "Class", "Order", "Family", "Genus", "Species", "Assembly"), ylim=c(0,1), col="grey")
boxplot(null_ITS_spearman$metric, ec_class_ITS$spearman, ec_order_ITS$spearman, ec_family_ITS$spearman, ec_genus_ITS$spearman, ec_species_ITS$spearman, ec_assembly_ITS$spearman,
        main="ITS Spearman", ylab="Spearman correlation", names=c("Null", "Class", "Order", "Family", "Genus", "Species", "Assembly"), ylim=c(0,1), col="grey")
boxplot(null_ITS_cosine$metric, ec_class_ITS$cosine, ec_order_ITS$cosine, ec_family_ITS$cosine, ec_genus_ITS$cosine, ec_species_ITS$cosine, ec_assembly_ITS$cosine,
        main="ITS Cosine", ylab="Cosine similarity", names=c("Null", "Class", "Order", "Family", "Genus", "Species", "Assembly"), ylim=c(0,1), col="grey")
boxplot(null_18S_cosine$metric, ec_class_18S$cosine, ec_order_18S$cosine, ec_family_18S$cosine, ec_genus_18S$cosine, ec_species_18S$cosine, ec_assembly_18S$cosine,
        main="18S Cosine", ylab="Cosine similarity", names=c("Null", "Class", "Order", "Family", "Genus", "Species", "Assembly"), ylim=c(0,1), col="grey")
