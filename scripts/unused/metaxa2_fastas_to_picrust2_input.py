#!/usr/bin/env python

from collections import defaultdict
from picrust2.util import read_fasta, write_fasta
import argparse
import hashlib
import os

__author__ = "Gavin Douglas"

parser = argparse.ArgumentParser(
                        description="Reads in FASTAs of rRNA parsed by metaxa2. "
                        "Each separate FASTA is assumed to be for a different sample. "
                        "Will dereplicate seqs and return single FASTA and ASV abundance TABLE"
                        "md5sum hashes of sequences will be output headers and sample name are assumed "
                        "to be first field after splitting basename of files by . and removing \"_metaxa2_out\"",
                        formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('FASTAs', metavar='FASTA', type=str, nargs='+',
                    help='metaxa2 FASTA outputs')

parser.add_argument("-o", "--output", help="Prefix of output files.",
                    required=True, type=str)


def main():

    args = parser.parse_args()

    # Initialize dict of all different seqs and set of all sample names.
    all_seqs = defaultdict(dict)
    samples = set()

    # Loop through all FASTAs and get dereplicated counts of sequences.
    for fastafile in args.FASTAs:
        
        samplename = os.path.basename(fastafile).split(".")[0]
        samplename = samplename.replace("_metaxa2_out", "")

        samples.add(samplename)

        seqs_in = read_fasta(fastafile).values()

        for sequence in seqs_in:
            
            if sequence not in all_seqs:
                all_seqs[sequence] = defaultdict(int)

            all_seqs[sequence][samplename] += 1

    # Get md5sum of each sequence and write out FASTA with md5sums as headers.
    # Also write out sequence abundance table with ids of the hashes as well.

    out_fasta = open(args.output + ".fasta", "w")
    out_table = open(args.output + ".tsv", "w")

    samples_ordered = list(samples)

    # Write table header.
    tab_header = ["sequence"] + samples_ordered
    print("\t".join(tab_header), file=out_table)

    for sequence in all_seqs.keys():
        seq_hash = hashlib.md5()
        seq_hash.update(sequence.encode())
        seq_hash = seq_hash.hexdigest()

        print(">" + seq_hash, file=out_fasta)
        print(sequence, file=out_fasta)

        # Get list of abundances for all samples for this sequence.
        table_line = [seq_hash]
        for sample in samples_ordered:
            table_line.append(str(all_seqs[sequence][sample]))

        print("\t".join(table_line), file=out_table)


    out_fasta.close()
    out_table.close()


if __name__ == "__main__":
    main()
