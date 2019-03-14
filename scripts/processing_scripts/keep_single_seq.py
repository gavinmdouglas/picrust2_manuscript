#!/usr/bin/python3

import argparse
import re
import sys
from os import path
from picrust2.util import read_fasta, write_fasta


def main():

    parser = argparse.ArgumentParser(

        description="Reads in FASTA file and keeps single sequence (or "
                    "possibly no sequence). Works by first screening out all "
                    "sequences of length less or greater than the lower and "
                    "upper bounds given. Will then screen out sequences with "
                    "greater than a set percent of Ns. Then preferentially keeps "
                    "sequences with lower proportion of Ns. Finally if there is "
                    "still a tie, will just choose a sequence randomly. Output file only "
                    "created if there is a final sequence to be written.",

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to FASTA file", required=True)

    parser.add_argument("-l", "--lower_length", metavar="SIZE", type=int,
                        required=True, help="Min length of each seq in FASTA")

    parser.add_argument("-u", "--upper_length", metavar="SIZE", type=int,
                        required=True, help="Max length of each seq in FASTA")

    parser.add_argument("-p", "--prop_n", metavar="SIZE", type=float,
                        required=True, help="Proportion of N characters "
                                            "permitted per-sequence")

    parser.add_argument("-o", "--output_dir", metavar="PATH", type=str,
                        required=True, 
                        help="Output directory to write final FASTA of single "
                             "sequence IF there is a sequence left to write "
                             "after filtering.")

    parser.add_argument("--rename_seq", action="store_true",
                        help="Flag to indicate that the sequence header should "
                             "be renamed to be the the first 2 fields of the "
                             "filenames after delimiting by \'_\'")

    parser.add_argument("--rename_seq_full", action="store_true",
                        help="Flag to indicate that the sequence header should "
                             "be renamed to be the full filename.")

    args = parser.parse_args()

    in_fasta = read_fasta(args.fasta)

    # If no sequences in file then stop job.
    if not in_fasta:
        sys.exit("Stopping - no sequences in file.")

    # Remove all sequences with length outside cut-off range or with greater
    # than the specified proportion of N characters.
    seq2remove = set()

    seq_N_pro = {}

    for seq_id, sequence in in_fasta.items():

        seq_len = len(sequence)

        N_pro = sequence.upper().count("N")/seq_len

        if seq_len < args.lower_length or seq_len > args.upper_length or N_pro > args.prop_n:
            seq2remove.add(seq_id)
        else:
            seq_N_pro[seq_id] = N_pro

    # Remove the specified sequences.
    for seq_id in seq2remove:
        del in_fasta[seq_id]

    # If no sequences in file then stop job.
    if not in_fasta:
        sys.exit("Stopping - no sequences left after filtering.")

    # Of remaining sequences figure out which has the lower proportion of Ns.
    # If there is a tie then the first sequence is taken (since dictionary keys
    # are unordered this results in a random selection).
    best_seq = None
    lowest_pro_N = 1.1

    for seq_id, pro_N in seq_N_pro.items():
        if pro_N < lowest_pro_N:
            best_seq = seq_id

    out_basename = path.splitext(path.basename(args.fasta))[0]

    outfile = path.join(args.output_dir, out_basename + "_best.fna")

    # Add the best sequence to a dictionary so it can be output easily.
    out_seq = {}

    # If rename_seq option set then replace current header with first 2 fields
    # of filename after delimiting by "_".
    if args.rename_seq:
        file_split = out_basename.split("_")
        seqname = file_split[0] + "_" + file_split[1]
    elif args.rename_seq_full:
        seqname = out_basename
    else:
        seqname = best_seq

    out_seq[seqname] = in_fasta[best_seq].upper()



    write_fasta(out_seq, outfile)

if __name__ == '__main__':
    main()
