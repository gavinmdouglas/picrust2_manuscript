A number of people have had trouble pairing the 16S and shotgun metagenomic (MGS) data per sample in the Human Microbiome Project (HMP) repository. [Vince Maffei](https://github.com/vmaffei) has put together a README for how this pairing can be done. The majority of the text below was written by him with minor modifications by me. You should first look here for an explanation for how the HMP samples are paired: http://hmpdacc.org/micro_analysis/microbiome_analyses.php.

### Summary
HMP shotgun MGS samples are given "SRS" IDs, which correspond to 16S sample "SRR" IDs in the V3V5_map file provided by HMP DACC. Note that rows in the V3V5 map table are individual 16S samples. The "SRS_SampleID" column in the V3V5 map table will be used to reference shotgun MGS samples only, which pair to the remaining 16S sample information in a given row. SRS IDs are assigned URLs to downloadable raw shotgun MGS read files in the HMP DACC file browser. SRR IDs will be used to download 16S read files from the NCBI Sequence Read Archive (SRA).

### Step 1: Modify the V3V5 map file (ppAll_V35_map.txt) for R parsing.

Download the file with wget:

```
wget http://hmpdacc.org/doc/ppAll_V35_map.txt
```

The column "HMPBodySiteHMPBodySubsite" in `ppAll_V35_map.txt` should be "HMPBodySite[tab]HMPBodySubsite". Adding a tab fixes a column name shift:

```
sed 's/HMPBodySubsiteHMPBodySite/HMPBodySubsite\tHMPBodySite/g' ppAll_V35_map.txt > 16SV35_meta_map.txt
```
### Step 2: Modify the MGS map file (HMIWGS_healthy.csv) for R parsing.

First go to HMP DACC (ex: http://hmpdacc.org/HMIWGS/healthy/) and click "Save as CSV" to download the tabled data. Name this `HMIWGS_healthy.csv` for the purposes of this tutorial.

Convert to tab delimited, remove quotes, and remove carriage returns.

```
sed 's/,/\t/g' HMIWGS_healthy.csv | sed 's/"//g' | sed 's/^M/\n/g'  > pre_processed_dl_list.txt
```

Note: you may need to re-type the carriage return character ^M manually by typing ctrl+V and then ctrl+M within the command line prompt.

### Step 3: Get overlapping ids, download 16S SRA files, and get MGS fastq URLs.

In R, filter `16SV35_meta_map.txt` for the SRS IDs listed in `pre_processed_dl_list.txt`, output master URLs of MGS fastqs, and download the 16S SRA files directly.

```R
# Load package required for working with SRA metadata.
library(SRAdb)

# Load package required for interacting with database.
library(DBI)

# Set the destination for the SRA db and the output directory for raw data.
srafile <- "SRAmetadb.sqlite"
destDir <- "raw_16S"

# Download metadata db if it doesn't already exist.
# (1894.2 MB when I downloaded it with this timestamp: 2017-06-16 23:12:31).
if(!file.exists(srafile)) srafile <<- getSRAdbFile()

# Create destination directory if it doesn't exist.
if(!file.exists(destDir)) dir.create(destDir)

# Read in map files.
meta_map <- read.table("pre_processed_dl_list.txt", header=TRUE, row.names=NULL,
                       sep='\t',stringsAsFactors=FALSE)

V3V5_map <- read.table("16SV35_meta_map.txt", header=TRUE, row.names=NULL,
                       sep='\t',stringsAsFactors=FALSE)

# Get column subsets and rename to be more descriptive.
meta_map_df <- data.frame("Meta_SRS_Sample_ID"=meta_map$SRS.ID,
                          "URL"=meta_map$Reads.File.Location,
                          stringsAsFactors=FALSE)

V3V5_16S_map_df <- data.frame("V3V5_Run"=V3V5_map$Run,
                              "V3V5_NAP_Sample_ID"=V3V5_map$NAP,
                              "Subject_ID"=V3V5_map$RSID,
                              "SRS_Sample_ID"=V3V5_map$SRS_SampleID,
                              "Sex"=V3V5_map$Sex, "Site"=V3V5_map$HMPBodySite,
                              "Subsite"=V3V5_map$HMPBodySubsite,
                              stringsAsFactors=FALSE)

# Get overlapping SRS ids.
paired <- intersect(meta_map_df$Meta_SRS_Sample_ID,
                    V3V5_16S_map_df$SRS_Sample_ID)

# Subset rows of each dataframe where SRS ids are within this overlapping set.
meta_map_df_sub <- meta_map_df[meta_map_df$Meta_SRS_Sample_ID %in% paired,]
V3V5_16S_map_df_sub <- V3V5_16S_map_df[V3V5_16S_map_df$SRS_Sample_ID %in% paired,]

# Get download links for each file.
meta_dl <- paste("ftp://public-ftp.hmpdacc.org",
                 gsub("/data","",meta_map_df_sub$URL), sep="")

V3V5_dl <- V3V5_16S_map_df_sub$V3V5_Run

# Write out download links for MGS files.
write.table(meta_dl, file="meta_final_url.txt", col.names=FALSE,
            row.names=FALSE, quote=FALSE, sep='\t')

# Now download the 16S SRR files from SRA using SRAdb in R.

# Connect to downloaded metadata database.
con = dbConnect(RSQLite::SQLite(), srafile)

# Loop over all ids and download 16S SRA files.
sapply(V3V5_dl, function(x) try(getSRAfile(x, sra_con=con, destDir=destDir),
                                silent=FALSE))

```
### Step 4: Download MGS files.

Loop over each MGS datafile and download with wget.

```sh
cat meta_final_url.txt | while read f; do wget ${f}; done;
```

### Step 5: Convert 16S .sra files to FASTQs.

```
fastq-dump --split-files -F --gzip --skip-technical -O ../fastq_files *sra
```
