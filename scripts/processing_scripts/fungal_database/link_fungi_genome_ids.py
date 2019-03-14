#!/usr/bin/python3

import sys
import os

# Note this script was written quickly for single usage and
# contains hard-coded paths. The purpose is to identify
# links between 1000 fungi genome ids and the filenames.
# This is necessary since the files are not named consistently.

def main():

    public_ids = set()

    possible_naming_to_ids = {}

    missing_ids = set()

    # Read in public ids.
    with open("/home/gavin/projects/picrust_pipeline/fungal_genomes/fungi_genome_assembly_id_links_no_mito_public.txt", "rt") as in_public_ids:
        for line in in_public_ids:
            line = line.rstrip()
            line_split = line.split()
            public_ids.add(line_split[-1])

            spec_name = line_split[0][0:3] + line_split[1][0:2]
            possible_naming_to_ids[spec_name] = line_split[-1]

    path_prefix = "/home/gavin/projects/picrust_pipeline/fungal_genomes/all_fungi"

    folders_to_check = ["1000_fungi_AA_all", "1000_fungi_EC_all",
                        "1000_fungi_genomes_all"]

    # Loop over all public ids and identify sub dir in all 3 folders.
    for public_id in public_ids:

        for folder in folders_to_check:

            expected_dir = os.path.join(path_prefix, folder, public_id)

            # Some folders not found across all 3 datatypes - which will be
            # ignored.
            if not os.path.isdir(expected_dir):
                print("This folder not found: " + expected_dir, file=sys.stderr)
                missing_ids.add(public_id)
                continue

            containing_files = os.listdir(expected_dir)

            for filename in containing_files:

                if filename == "tmp":
                    continue

                if not match_filename_subsets_w_set(filename, public_ids):

                    filename = filename.replace(".", "_")
                    file_split = filename.split("_")
                    Spec = file_split[0][0:3] + file_split[1][0:2]

                    if Spec in possible_naming_to_ids:
                        print("No obvious public id for " + filename + ", but "
                              "possibly " + possible_naming_to_ids[Spec] +
                              " in " + folder)
                    else:
                        print("No obvious public id for " + filename + " in " +
                              folder)


def match_filename_subsets_w_set(filename, in_set):

    filename = filename.replace(".", "_")
    file_split = filename.split("_")

    current = file_split[0]
    i = 1
    while i <= len(file_split) - 1:
        if current in in_set:
            return(True)
        current = current + "_" + file_split[i]
        i += 1


    return(False)


if __name__ == '__main__':
    main()
