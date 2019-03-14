#!/usr/bin/python3

import argparse
import sys
import os
import gzip
from collections import defaultdict

# Define class that is counter of functional ids for a given genome.
class genome_function_counts():

    def __init__(self, genome):
        self.genome = genome
        self.__counts = defaultdict(int)

    def add_function_count(self, function_id):
        self.__counts[function_id] += 1

    def return_function_count(self, function_id):
        return self.__counts[function_id]


def main():

    parser = argparse.ArgumentParser(
    description="Script to parse out all EC number counts for a fungi genome",

    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PATH",
                        type=str, required=True,
                        help="Path to folder with all count tables in subdirectories.")

    parser.add_argument("-g", "--genome_ids", metavar="PATH",
                        type=str, required=True,
                        help="Path to file with public genome ids in 2nd column")

    args = parser.parse_args()

    public_ids = set()

    # Read in all public genome ids.
    with open(args.genome_ids, "rt") as genome_ids_file:
        for line in genome_ids_file:
            public_ids.add(line.split()[-1])

    # Initialize set containing all functions of interest and a dictionary
    # containing genome_function_counts classes that keep track of the counts
    # of each category for each genome.
    unique_functions = set()
    genome2counts = {}

    for public_id in public_ids:

        genome2counts[public_id] = genome_function_counts(public_id)

        expected_dir = os.path.join(args.input, public_id)

        if not os.path.isdir(expected_dir):
            print("Expected subdirectory not found: " + expected_dir,
                  file=sys.stderr)
            continue

        subdir_files = os.listdir(expected_dir)

        if len(subdir_files) != 1:
            sys.exit("Error - multiple files in " +
                     os.path.join(args.input, public_id))

        infile = os.path.join(args.input, public_id, subdir_files[0])

        lc = 0
        ec_col = None
        protein_col = None

        # Set of past EC numbers and protein combos so there aren't repeats.
        past_protein_ec = set()

        with gzip.open(infile, "rt") as file_in:
            for line in file_in:
                if lc == 0:
                    header_split = line.split("\t")
                    ec_i = [i for i, s in enumerate(header_split) if 'ecNum' in s]

                    if len(ec_i) != 1:
                        sys.exit("No single EC column identified for this header: " + line)

                    ec_col = ec_i[0]
                    lc += 1
                    continue

                line_split = line.split("\t")

                protein_id = line_split[0]

                ec_num = line_split[ec_col]

                if "-" in ec_num:
                    continue

                # Skip if already come across this protein / EC # combination.
                combo = protein_id + "|" + ec_num
                if combo in past_protein_ec:
                    continue
                else:
                    past_protein_ec.add(combo)

                # Add "EC:" to front of id.
                ec_num = "EC:" + ec_num

                genome2counts[public_id].add_function_count(ec_num)
                unique_functions.add(ec_num)

    # Print header line.
    print("assembly", "\t".join(list(unique_functions)), sep="\t")

    # Loop through all genomes and output counts of each function id.
    for public_id in public_ids:

        output_line = [public_id]

        for function in unique_functions:

            output_line.append(str(genome2counts[public_id].return_function_count(
                                                                    function)))

        print("\t".join(output_line))


if __name__ == '__main__':
    main()
