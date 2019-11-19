#!/usr/bin/python3

import argparse
import re
import sys
from os import path
import biom
from picrust2.util import read_fasta, make_output_dir, biom_to_pandas_df


def main():

    parser = argparse.ArgumentParser(

        description="Creates output FASTA for each sample with each ASV repeated for every count in that sample.",

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to full FASTA file.", required=True)

    parser.add_argument("-b", "--biom", metavar="BIOM", type=str,
                        help="Path to BIOM table.", required=True)

    parser.add_argument("-o", "--outdir", metavar="PATH", type=str,
                        help="Name of folder to make for output files.", required=True)

    args = parser.parse_args()

    in_fasta = read_fasta(args.fasta)

    in_table = biom_to_pandas_df(biom.load_table(args.biom))

    # If no sequences in file then stop job.
    if not in_fasta:
        sys.exit("Stopping - no sequences in file.")

    make_output_dir(args.outdir)

    for sample in in_table.columns:
        sample_outfile = args.outdir + "/" + sample + ".fasta"

        sample_outfh = open(sample_outfile, 'wt')

        for asv in in_table.index.values:
            asv_count = in_table.loc[asv, sample]
            if asv_count > 0:
                for i in range(int(asv_count)):
                    print(">" + asv + "_" + sample + "_" + str(i), file=sample_outfh)
                    print(in_fasta[asv], file=sample_outfh)

        sample_outfh.close()


if __name__ == '__main__':
    main()
