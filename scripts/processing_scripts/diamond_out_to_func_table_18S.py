#!/usr/bin/python3

import argparse
import sys
from os import listdir
from os.path import isfile, join, basename
import gzip
from collections import defaultdict


# Define class that is counter of hits for a given diamond output.
class diamond_hit_counts():

    def __init__(self, name):
        self.name = name
        self.__counts = defaultdict(int)

    def add_hit_counts(self, hit_id):
        self.__counts[hit_id] += 1

    def return_hit_count(self, hit_id):
        return self.__counts[hit_id]


def main():

    parser = argparse.ArgumentParser(description="Script to parse directory \
of diamond output files (WITH ONLY THE TOP HIT PER QUERY) and builds \
abundance table of hits. Used to parse 18S diamond output.", 

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-d", "--diamond_out", metavar="IN_FOLDER",
                        type=str, help="Path to folder with diamond outfiles.",
                        required=True)

    parser.add_argument("--id", metavar="PERCENT_IDENTITY", type=float,
                        help="Min percent identity", required=True)

    args = parser.parse_args()

    diamond_basenames = [f for f in listdir(args.diamond_out) 
                        if isfile(join(args.diamond_out, f))]

    diamond_outfiles = [join(args.diamond_out, f) for f in diamond_basenames]

    file2counts = {}

    unique_hits = set()

    for i, diamond_outfile in enumerate(diamond_outfiles):

        file2counts[diamond_basenames[i]] = diamond_hit_counts(diamond_basenames[i])

        previous_query = set()

        with open(diamond_outfile, "rt") as diamond_outfile_read:
            for line in diamond_outfile_read:

                # Remove line terminator from end of line and split on tabs.
                line = line.rstrip("\r\n")
                line_split = line.split("\t")

                query_name = line_split[0]
                
                percent_id = line_split[3]

                # Skip match if it does not pass the identity setting.
                if float(percent_id) < args.id:
                    continue

                # Skip query sequence if it has already been parsed.
                if query_name in previous_query:
                    continue
                else:
                    previous_query.add(query_name)                 

                hit_name = line_split[1].split("|")[0]

                unique_hits.add(hit_name)

                # Add hit name as count for diamond outfile.
                file2counts[diamond_basenames[i]].add_hit_counts(hit_name)
          

    # Print header line.
    header_out = ["assembly"]
    for outfile in diamond_basenames:
        header_out += [outfile.replace("_protein.faa_diamond_out.txt", "")]
    print("\t".join(header_out))

    for hit_name in unique_hits:
        output_line = [hit_name]
        for outfile in diamond_basenames:
            output_line.append(str(file2counts[outfile].return_hit_count(
                                                       hit_name)))
        print("\t".join(output_line))


if __name__ == '__main__':
    main()
