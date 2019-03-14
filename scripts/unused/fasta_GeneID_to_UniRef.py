#!/usr/bin/python3

import argparse
import os
import pandas as pd
import pickle
import re
from collections import defaultdict
import pprint

# Define class that contains counts of all functions for a given genome.
class genome_function_counts():

    def __init__(self, genome):
        self.genome = genome
        self.__functions = defaultdict(int)

    def add_function_count(self, function):
        self.__functions[function] += 1

    def get_function_count(self, function):
        return self.__functions[function]

    def get_function_names(self):
        return self.__functions.keys()


def main():

    parser = argparse.ArgumentParser(
        description="Reads in pickled dictionary mapping GeneID to UniRef "
                    "ids. Then reads in directory of FASTAs corresponding to "
                    "CDS regions from genomes that contain RefSeq GeneIDs in "
                    "headerlines. Will output table of counts of each type of "
                    "UniRef id per genome.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-d", "--dict", metavar="DICT", type=str,
                        help="Path to pickled dictionary file.", required=True)

    parser.add_argument("-i", "--input", metavar="DIRECTORY", type=str,
                        help="Path to directory of CDS FASTAs to parse.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="FILE", type=str,
                        help="Path to output table.", required=True)

    args = parser.parse_args()

    # Read in pickled dictionary.
    with open(args.dict, "rb") as input_dict:
        geneid2uniref = pickle.load(input_dict)

    # Read through all input FASTAs and get GeneIDs from headerlines.
    in_fastas = os.listdir(args.input)

    # Dict to keep track of genome_function_counts classes for each genome.
    genome_counts = defaultdict(genome_function_counts)
    genome_ids = []

    for fasta in in_fastas:

        # Get genome id from FASTA file name.
        genome_id_match = re.search(r'(GCF_[^_]+)_', os.path.basename(fasta))
        genome_id = genome_id_match.group(1)
        genome_counts[genome_id] = genome_function_counts(genome_id)
        genome_ids.append(genome_id)

        # Read through FASTA, identify GeneIDs, and get matching UniRef ids.
        with open(os.path.join(args.input, fasta), "r") as input_fasta:
            for line in input_fasta:
                if line[0] == ">":

                    # If header-line then parse out GeneID.
                    re_match = re.search(r'GeneID:(\d+)]', line)

                    if re_match is None:
                        continue

                    try:
                        GeneID_match = re_match.group(1)
                    except IndexError as error:
                        continue

                    if GeneID_match in geneid2uniref:
                        matching_uniref = list(geneid2uniref[GeneID_match])
                        for uniref in matching_uniref:
                            print(GeneID_match + " " + uniref)
                            genome_counts[genome_id].add_function_count(uniref)

    all_uniref = []
    for genome in genome_counts.keys():
        all_uniref += list(genome_counts[genome].get_function_names())
    all_uniref = list(set(all_uniref))

    output_file = open(args.output, "w")

    # Print header-line.
    header_line = ["genome"] + genome_ids
    print("\t".join(header_line), file=output_file)

    for function in all_uniref:
        out = [function]
        for genome in genome_ids:
            out.append(str(genome_counts[genome].get_function_count(function)))
        print("\t".join(out), file=output_file)

    output_file.close()

if __name__ == '__main__':
    main()
