# Get basic counts of number of centroid sequences per cluster.
# Also generate distribution of all sequence lengths.
# Get summary statistics to report in text as well.

# Read in plotrix, which can be used to put breaks in axis.
library(plotrix)
library(Biostrings)

# Read in file containing clusters.
clusters <- read.delim("/home/gavin/projects/picrust_pipeline/IMG_pipeline_prep/ssu_align_pipeline/new_cluster_files/ssu_align_out.bacteria.mask_clusters_no-low-qual.txt",
                       header=F, sep="\t", stringsAsFactors=FALSE)


# Need to remove the 6 eukaryotic clusters.

euk_seqs <- c("2503982003", "2170459028", "2518645510",
              "2524614788", "2524614793", "2751185563")

clusters <- clusters[-which(clusters$V1 %in% euk_seqs),, drop=FALSE]

dim(clusters)


### 20000 clusters

cluster_sizes <- sapply(clusters$V1, function(x) { length(strsplit(x, " ")[[1]]) })
length(which(cluster_sizes > 1))

### Only 3002 clusters are of more than 1 sequence

summary(cluster_sizes)
### Range of 1 to 1379 sequences per cluster and mean of 2.096.

# Histogram of sequences per cluster (of clusters with more than 1 sequence).
cluster_hist <- hist(cluster_sizes[which(cluster_sizes > 1)], breaks=100, plot = FALSE)

gap.barplot(cluster_hist$counts, gap=c(200, 2590), xlab="Number of sequences in clusters",
            ylab="Number of clusters", main="", ytics=c(seq(from=0, to=150, by=50), seq(from=2600, to=3000, by=100)),
            xtics=cluster_hist$breaks)

### Read in FASTA.
cluster_seqs <- readDNAStringSet("/home/gavin/projects/picrust_pipeline/IMG_pipeline_prep/ssu_align_pipeline/ssu_align_out.bacteria.mask.derep.no-low-qual.afa")

# Remove gap characters.
cluster_seqs_nogap <- lapply(cluster_seqs, function(x){ gsub("-", "", x)})

cluster_seqs_nogap_lengths <- sapply(cluster_seqs_nogap, nchar)

### mean: 1489.563
### sd: 65.81238

hist(cluster_seqs_nogap_lengths, xlim=c(1200, 1650), ylim=c(0, 6250), main="", xlab="Sequence lengths", col="dark blue")
