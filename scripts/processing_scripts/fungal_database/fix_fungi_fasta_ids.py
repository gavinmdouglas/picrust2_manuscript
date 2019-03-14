#!/usr/bin/python3

import sys
import os
import argparse

def main():

    parser = argparse.ArgumentParser(
        description="Reads in FASTA file, file with genome ids, and file with "
                    "unusual ids mapped to genome ids. Will correct all FASTA "
                    "headers to be correct genome ids.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to FASTA file", required=True)

    parser.add_argument("--genome_ids", metavar="FILE", type=str,
                        help="Path to file with genome ids in 2nd column.",
                        required=True)

    parser.add_argument("--unusual", metavar="FILE", type=str,
                        help="Path to file with unusual ids mapped to genome "
                             "ids", required=True)

    args = parser.parse_args()

    public_ids = set()
    unusual_ids = {}

    # Read in public ids.
    with open(args.genome_ids, "rt") as in_public_ids:
        for line in in_public_ids:
            line = line.rstrip()
            line_split = line.split()
            public_ids.add(line_split[-1])

    # Read in mapping of unusual ids to genome ids.
    with open(args.unusual, "rt") as in_unusual:
        for line in in_unusual:
            line = line.rstrip()
            line_split = line.split()
            unusual_ids[line_split[0]] = line_split[1]

    # Read through FASTA and fix ids.
    with open(args.fasta, "rt") as in_fasta:
        for line in in_fasta:

            line = line.rstrip()

            # Print out non-header lines.
            if line[0] != ">":
                print(line)
                continue
            
            header = line[1:]

            header_match_id = match_str_subsets_return_match(header, public_ids)

            if header_match_id:
                print(">" + header_match_id)
            else:
                header_alt_match = match_str_subsets_return_match(header, set(list(unusual_ids.keys())))

                if header_alt_match:
                    print(">" + unusual_ids[header_alt_match])
                else:
                    sys.exit("No matched if for this line: " + line)




def match_str_subsets_return_match(in_str, in_set):

    in_str = in_str.replace(".", "_")
    str_split = in_str.split("_")

    current = str_split[0]
    if current in in_set:
        return(current)

    i = 1
    while i <= len(str_split) - 1:
        current = current + "_" + str_split[i]
        if current in in_set:
            return(current)
        i += 1


    return(False)


if __name__ == '__main__':
    main()
