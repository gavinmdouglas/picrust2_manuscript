### Simple script to regroup a TSV-formatted BIOM file.
library("Biostrings")
library("optparse")

option_list <- list(
  
  make_option(c("--BIOM"), type="character", default=NULL,
              help="Tab-delimited BIOM file." , metavar="path"),

  make_option(c("--BLAST"), type="character", default=NULL,
              help="BLAST output in BLAST6out format." , metavar="path"),
  
  make_option(c("--db_match"), type="character", default=NULL,
              help="FASTA of database sequences hit by query sequences." , metavar="path"),
  
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Output re-grouped tab-delimited BIOM file." , metavar="path")
)

opt_parser <- OptionParser(
  option_list=option_list,

  usage = "%prog [options] --BIOM PATH --BLAST PATH --db_match PATH --out PATH",

  description = "Re-regroup BIOM file based on matching ids in BLAST output. Will also run some sanity checks like making sure all the re-grouped ids are found in a FASTA file of matching reference ids (that is expected to have been created along with the BLAST output table)."
)

opt <- parse_args(opt_parser)

BIOM_infile <- opt$BIOM
BLAST_outfile <- opt$BLAST
ref_FASTA_outfile <- opt$db_match
regrouped_outfile <- opt$out

if(is.null(BIOM_infile) || is.null(BLAST_outfile) || is.null(ref_FASTA_outfile) || is.null(regrouped_outfile)) {
  stop("All four arguments must be set.")
}
  
# First figure out if the table has an extra line ("# Constructed from biom file") at the top or not to figure out how to read in the entire table.
biom_head <- read.table(BIOM_infile, header=FALSE, sep="\t", nrows = 1, comment.char="", stringsAsFactors = FALSE)

if(biom_head$V1[1] == "# Constructed from biom file") {
 biom_in <- read.table(BIOM_infile, header=TRUE, sep="\t", comment.char="", skip = 1, row.names=1)
} else if(biom_head$V1[1] == "#OTU ID") {
  biom_in <- read.table(BIOM_infile, header=TRUE, sep="\t", comment.char="", row.names=1)
} else {
 stop("Stopping - expected table to begin with either \"# Constructed from biom file\" or \"#OTU ID\".") 
}

# Read in BLAST output and make sure only 1 reference hit taken per query sequence.
BLAST_out <- read.table(BLAST_outfile, header=FALSE, sep="\t", stringsAsFactors = FALSE)

BLAST_out$V1 <- as.character(BLAST_out$V1)
BLAST_out$V2 <- as.character(BLAST_out$V2)

if(length(which(duplicated(BLAST_out$V1))) > 0) {
  stop("More than 1 match per query sequence in file - clean-up file first before running this command.")
}

rownames(BLAST_out) <- BLAST_out$V1

# Remove sequences with no reference hits.
missing_seqs <- which(! rownames(biom_in) %in% rownames(BLAST_out))

if(length(missing_seqs) > 0) {
  biom_in <- biom_in[-missing_seqs, ]
}

biom_in$other <- BLAST_out[rownames(biom_in), "V2"]
biom_in_other <- aggregate(. ~ other, data=biom_in, FUN=sum)
rownames(biom_in_other) <- biom_in_other$other
biom_in_other <- biom_in_other[, -which(colnames(biom_in_other) == "other")]


# Check that all reference ids are also present in the FASTA file that was created outside of this program.
ref_FASTA_seq_ids = names(readDNAStringSet(ref_FASTA_outfile))
missing_FASTA_ids <- ref_FASTA_seq_ids[which(! ref_FASTA_seq_ids %in% rownames(biom_in_other))]

if(length(missing_FASTA_ids) > 0) {
  print(missing_FASTA_ids)
  stop("Stopping - not all ids in FASTA found in regrouped BIOM table (see above).") 
}

# Write out re-grouped BIOM table.
sample_ids <- colnames(biom_in_other)
biom_in_other[, "#OTU ID"] <- rownames(biom_in_other)
biom_in_other <- biom_in_other[, c("#OTU ID", sample_ids)]

write.table(file = regrouped_outfile, x = biom_in_other, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")


### Run sanity checks on a few mappings to make sure regrouping done correctly.
### This was done when regrouping the HMP1 files.
### 1 ASV mapped to 885560
### 3 ASVs mapped to 888466
### 4 ASVs mapped to 122517

# biom_in_885560 <- biom_in[which(biom_in$other == "885560"), -which(colnames(biom_in) == "other")]
# biom_in_888466 <- biom_in[which(biom_in$other == "888466"), -which(colnames(biom_in) == "other")]
# biom_in_122517 <- biom_in[which(biom_in$other == "122517"), -which(colnames(biom_in) == "other")]
# 
# identical(colSums(biom_in_885560), colSums(biom_in_other["885560", ]))
# identical(colSums(biom_in_888466), colSums(biom_in_other["888466", ]))
# identical(colSums(biom_in_122517), colSums(biom_in_other["122517", ]))
