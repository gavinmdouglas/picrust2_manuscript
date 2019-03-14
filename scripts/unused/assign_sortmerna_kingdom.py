#!/usr/bin/python3

import argparse
import sys
import os


def main():

    parser = argparse.ArgumentParser(
        description="Determine Kingdom of reads extracted by SortMeRNA, "
                    "based on the BLAST outputfile. Will return mixed if "
                    "multiple kingdoms are hit.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-b", "--blast", metavar="PATH",
                        type=str, help="Path to SortMeRNA BLAST outfile.",
                        required=True)

    parser.add_argument("-t", "--taxa", metavar="PATH",
                        type=str, help="Table of reference sequence taxonomy.",
                        required=True)

    args = parser.parse_args()

    # Dictionary to map reference sequences to kingdoms.
    ref_kingdom = {}

    # Read in Kingdom per reference sequence.
    with open(args.taxa, 'r') as taxa_in:
        for line in taxa_in:
            line_split = line.split()
            
            seq_name = line_split[0]
            kingdom = line_split[1].split(";")[0]

            ref_kingdom[seq_name] = kingdom

    # Dictionary for keeping track of kingdom of query sequences.
    query_kingdom = {}

    with open(args.blast, 'r') as blast_in:
        for line in blast_in:
            line_split = line.split()

            read_name = line_split[0]
            kingdom = ref_kingdom[line_split[1]]

            # If last line also corresponded to same read then check both hits
            # are in same kingdom.
            if read_name in query_kingdom:
                if query_kingdom[read_name] != kingdom:
                    query_kingdom[read_name] = "mixed"
            else:
                query_kingdom[read_name] = kingdom


    outfile = os.path.basename(args.blast) + ".kingdom.tsv"
    out_read = open(outfile, 'w')

    total_count = {"Archaea": 0,
                   "Bacteria": 0,
                   "Eukaryota": 0,
                   "mixed": 0} 

    out_read.write("query\tkingdom\n")

    for query, kingdom in query_kingdom.items():
        out_read.write(query + "\t" + kingdom + "\n")
        total_count[kingdom] += 1

    out_read.close()

    print(args.blast, total_count["Archaea"], total_count["Bacteria"],
          total_count["Eukaryota"], total_count["mixed"], sep="\t")

if __name__ == '__main__':
    main()
