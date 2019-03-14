#!/usr/bin/python3

import argparse
import sys
import os
import gzip
from collections import defaultdict


def file2set(filename):

    '''Reads in a text file and returns a set containing each unique line.'''

    lines2return = set()

    with open(filename, "r") as file_in:
        for line in file_in:

            # If first char is comment then skip.
            if line[0] == "#":
                continue

            # Remove line terminator from end of line.
            line = line.rstrip("\r\n")

            lines2return.add(line)

    return lines2return


def read_mapfile(filename, ignore_2nd=False, higher2lower=True):

    '''Reads in wide format mapping file (tab-delimited), which can be gzipped.
    Will return dictionary mapping first column to list of all other columns
    per row if higher2lower=True. If False then will return a dictionary with
    the other columns as keys and the first column as values. The second
    column will be removed if ignore_2nd=True.'''

    filename_str, extension = os.path.splitext(filename)

    mapping = defaultdict(list)

    def parse_mapfile_line(line_in):

        line_in = line_in.rstrip("\n\r")
        line_split = line_in.split("\t")

        # Remove second element if specified.
        if ignore_2nd:
            del line_split[1]

        # Add ids to mapping dictionary.
        if higher2lower:
            mapping[line_split.pop(0)] = line_split
        else:
            higher_func = line_split.pop(0)
            for mapped_id in line_split:
                mapping[mapped_id] += [higher_func]

    if extension == ".gz":

        with gzip.open(filename, "rt") as file_in:

            for line in file_in:

                parse_mapfile_line(line)

    else:

        with open(filename, "rt") as file_in:

            for line in file_in:

                parse_mapfile_line(line)

    return mapping


def parse_gff_function_ids(line_info, category):
    '''Takes in 9th column of gff3 when split on tabs and the functional
    category that should be returned. Returns a list all functional
    ids that the annotated line matches as a list.'''

    # If category is EC_number than expect it to be elsewhere in info column.
    if category == "EC_number":

        # Check whether EC_number is in expected position.
        if(line_info[0:9] != "EC_number"):
            print("EC_number not first characters", file=sys.stderr)
            return([])

        # Exclude any EC_number containing an "n" (preliminary) or "-"
        # (higher-order) character.
        ec_nums = []
        for ec in line_info.split(";")[0][10:].split(","):
            if "n" not in ec and "-" not in ec:
                ec_nums.append(ec)
        return(ec_nums)

    else:

        # Split string on ";" and find section that starts with "Dbxref=".
        info_split = line_info.split(";")

        ids = []

        for section in info_split:

            # Find section starting with Dbxref=.
            if section[0:7] == "Dbxref=":

                # Once found, separate all function by splitting on ",".
                for function in section[7:].split(","):
                    if category in function:

                        # Finally, if function is category of interest then
                        # remove category label and add to list that will be
                        # returned. Note that KEGG ids have ":" in the id name
                        # so return function name joined on ":".

                        function_split = function.split(":")
                        function_split.pop(0)

                        ids.append(":".join(function_split))

        # Return unique ids.
        return(list(set(ids)))


# Define class that is counter of functional ids for a given gff file.
class gff_function_counts():

    def __init__(self, filename):
        self.filename = filename
        self.__counts = defaultdict(int)

    def add_function_count(self, function_id):
        self.__counts[function_id] += 1

    def return_function_count(self, function_id):
        return self.__counts[function_id]


def main():

    parser = argparse.ArgumentParser(description="Script to parse directory \
of gff3 files and parse out counts of functions into table based on \
annotated CDS information.", epilog='''
Usage example:

python3 gff3s_to_func_table.py -a /path/to/gffs -f uniprot -m \
uniprot2uniref90.tsv.gz > uniref90_counts.tsv

 ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-a", "--annotation_path", metavar="annotation_path",
                        type=str, help="Path to folder with all annotations.",
                        required=True)

    parser.add_argument("-f", "--function_category", metavar="category",
                        type=str, help="Functional category that should \
be parsed. One of: EC_number, uniprot, Pfam, TIGRFAM, GO, eggNOG",
                        choices=["EC_number", "uniprot", "Pfam", "TIGRFAM",
                                 "GO", "eggNOG", "KEGG"], required=True)

    parser.add_argument("-g", "--gff2keep", metavar="gff2keep", type=str,
                        default=None, help="File containing gffs to keep. \
This script assumes that all filenames have the form \"id.gff\". Specify \
only the id if gffs are in subfolders otherwise specify the gff filenames \
(without full path)")

    parser.add_argument("--in_subfolder", action="store_true",
                        help="Flag to indicate that all gff3 files are in \
subdirectories in the form\"id/id.gff\".")

    parser.add_argument("--transpose", action="store_true",
                        help="Flag to indicate that the output should be \
genomes as columns and functions as rows.")

    parser.add_argument("-m", "--mapping", metavar="mapping", type=str,
                        default=None, help="If you would like to convert \
the category in the gff files to a different category, then you can specify \
the mapping file here. The mapping file needs to be a wide format mapping \
file to get category of interest.")

    parser.add_argument("--ignore_2nd_col", action="store_true",
                        help="Flag to ignore second column in mapping file \
(file specified with --mapping option).")

    args = parser.parse_args()

    # Get dictionary mapping category in gffs to category of interest if
    # mapping file specified.
    if args.mapping:

        if args.ignore_2nd_col:
            ignore_second_col = True
        else:
            ignore_second_col = False

        category_map = read_mapfile(args.mapping,
                                    ignore_2nd=ignore_second_col,
                                    higher2lower=False)

    # If gff2keep file specified then read in files or ids.
    # Otherwise get list of all files in directory specified.
    if args.gff2keep:
        gff_files = list(file2set(args.gff2keep))
    else:
        gff_files = os.listdir(args.annotation_path)

    # If "--in-subfolder" option specified then alter list to point to files
    # in format "folder/folder.gff". Otherwise just add full path to filename.
    if args.in_subfolder:
        gff_full_path = [args.annotation_path + "/" + folder + "/" + folder +
                         ".gff" for folder in gff_files]
    else:
        gff_full_path = [args.annotation_path + "/" + file
                         for file in gff_files]

    # Initialize list containing all functions of interest and a dictionary
    # containing gff_function_counts classes that keep track of the counts
    # of each category for each gff.
    unique_functions = set()
    gff2counts = {}

    # Loop through gff files and keep count of category of interest for each
    # one.
    for i, gff in enumerate(gff_full_path):

        # Initialize gff's counter class.
        gff2counts[gff_files[i]] = gff_function_counts(gff_files[i])

        with open(gff, "rt") as gff3:
            for line in gff3:

                # If first char is comment then skip.
                if line[0] == "#":
                    continue

                # Remove line terminator from end of line and split on tabs.
                line = line.rstrip("\r\n")
                line_split = line.split("\t")

                # Skip line if it isn't a CDS annotation.
                if line_split[2] != "CDS":
                    continue

                # Skip line if category of interest isn't in info column.
                if args.function_category not in line_split[8]:
                    continue

                # Parse out function ids from 9th column.
                # Note that there will usually only be one id returned.
                function_ids = parse_gff_function_ids(line_split[8],
                                                      args.function_category)

                # If mapping file specified then convert ids to category of
                # interest (taking unique ids only).
                if args.mapping:

                    converted_ids = []
                    for function_id in function_ids:
                        if function_id in category_map:
                            converted_ids += category_map[function_id]

                    function_ids = list(set(converted_ids))

                # Loop through functional ids and add count to gff.
                for function_id in function_ids:

                    # Add function count.
                    gff2counts[gff_files[i]].add_function_count(
                                            function_id)

                    # Add function to set of all functions parsed so far.
                    unique_functions.add(function_id)

    if not args.transpose:

        # Print header line
        print("assembly", "\t".join(list(unique_functions)), sep="\t")

        # Loop through all gff files and output counts of each function id.
        for gff in gff_files:

            # Add assembly id to output line and make sure ".gff" is removed.
            output_line = [gff.replace(".gff", "")]

            for function in unique_functions:

                output_line.append(str(gff2counts[gff].return_function_count(
                                                       function)))

            print("\t".join(output_line))

    elif args.transpose:

        # Print header line.
        header_out = ["assembly"]
        for gff in gff_files:
            header_out += [gff.replace(".gff", "")]

        print("\t".join(header_out))

        for function in unique_functions:

            # Add assembly id to output line and make sure ".gff" is removed.
            output_line = [function]

            for gff in gff_files:

                output_line.append(str(gff2counts[gff].return_function_count(
                                                       function)))
            print("\t".join(output_line))


if __name__ == '__main__':
    main()
