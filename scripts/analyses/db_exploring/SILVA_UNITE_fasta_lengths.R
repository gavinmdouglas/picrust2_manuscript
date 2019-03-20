### Get range of 18S and ITS sequences in SILVA and UNITE databases respectively.

### SILVA database used was: v132 SSURef (2,090,668 seqs)
### UNITE database used was: v7.2 - https://doi.org/10.15156/BIO/587475 (30,695 seqs)

library(Biostrings)

silva_in <- readDNAStringSet("/home/gavin/tmp/18S_ITS_tmp_db/SILVA_132_SSURef_tax_silva.fasta")
unite_in <- readDNAStringSet("/home/gavin/tmp/18S_ITS_tmp_db/sh_general_release_dynamic_01.12.2017.fasta")

silva_in_lengths <- sapply(silva_in, nchar)
unite_in_lengths <- sapply(unite_in, nchar)

# Get min and max for each database.
summary(silva_in_lengths)
# 605 and 3076

summary(unite_in_lengths)
# 146 and 2570
