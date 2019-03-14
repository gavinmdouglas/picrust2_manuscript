### Exploring whether the outlier KOs were genomes or individual KOs.

ec <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/past_tables/ec.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

ko <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ko.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

cog <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/cog.txt.gz"),
                 row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

pfam <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/pfam.txt.gz"),
                  row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

tigrfam <- read.table(gzfile("/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/tigrfam.txt.gz"),
                   row.names=1, header=T, sep="\t", stringsAsFactors = FALSE, check.names=FALSE)

# Set 10 to be max in all tables.
ec[ec > 10] <- 10
ec_orig_col <- colnames(ec)
ec$assembly <- rownames(ec)
ec <- ec[, c("assembly", ec_orig_col)]
write.table(x = ec, file = "/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ec.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

ko[ko > 10] <- 10
ko_orig_col <- colnames(ko)
ko$assembly <- rownames(ko)
ko <- ko[, c("assembly", ko_orig_col)]
write.table(x = ko, file = "/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/ko.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

cog[cog > 10] <- 10
cog_orig_col <- colnames(cog)
cog$assembly <- rownames(cog)
cog <- cog[, c("assembly", cog_orig_col)]
write.table(x = cog, file = "/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/cog.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

pfam[pfam > 10] <- 10
pfam_orig_col <- colnames(pfam)
pfam$assembly <- rownames(pfam)
pfam <- pfam[, c("assembly", pfam_orig_col)]
write.table(x = pfam, file = "/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/pfam.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

tigrfam[tigrfam > 10] <- 10
tigrfam_orig_col <- colnames(tigrfam)
tigrfam$assembly <- rownames(tigrfam)
tigrfam <- tigrfam[, c("assembly", tigrfam_orig_col)]
write.table(x = tigrfam, file = "/home/gavin/github_repos/picrust_repos/picrust2/default_files/prokaryotic/tigrfam.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

