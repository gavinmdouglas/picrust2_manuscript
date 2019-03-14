#!/usr/bin/python3

import argparse
import os
import re

def main():

    parser = argparse.ArgumentParser(
        description="Script to parse IMG genome ids and species names from "
                    "FASTA (from IMG)",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PATH",
                        type=str, help="Input fasta",
                        required=True)

    args = parser.parse_args()

    past_ids = set()

    print("IMG_id" + "\t" + "Species")

    with open(args.input, "r") as in_fasta:

            for line in in_fasta:

                line = line.rstrip()

                if line[0] == ">":

                    tax_start = line.index("[")
                    tax_end = line.index("]")

                    species = line[tax_start + 1:tax_end]

                    line_split = line.split("_")

                    genome_id = line_split[1]

                    if genome_id not in past_ids:

                        print(genome_id + "\t" + species)

                        past_ids.add(genome_id)

if __name__ == '__main__':
    main()
