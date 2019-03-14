#!/usr/bin/python3

import argparse
from collections import defaultdict
import pickle
import sys


def main():

    parser = argparse.ArgumentParser(
        description="Reads in UniRef mapfile of UniRef ids to NCBI GeneIDs. "
                    "Will output picked dictionary of GeneID mappings to "
                    "UniRef ids, which can be either UniRef50 or UniRef90.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-m", "--map", metavar="MAP", type=str,
                        help="Path to UniRef mapfile.", required=True)

    parser.add_argument("-c", "--category", choices=['UniRef50', 'UniRef90'],
                        type=str, help="Path to UniRef mapfile.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        help="Path to output pickled file.", required=True)

    args = parser.parse_args()

    # Inititalize dictionaries for keeping track of links between uniref ids.
    uniparc2uniref = {}

    # Initialize final dictionary that will be pickled.
    geneid2uniref = defaultdict(set)

    # Inititalize set for keeping track of UniParc ids that do not seem to be
    # linked with a UniRef id. These ids will be used for a sanity check at the
    # end.
    missing_uniref = set()

    # Read through mapfile. The key assumption is that UniParc maps with
    # GeneIDs will always be below UniParc maps with UniRef ids.
    with open(args.map, "r") as mapfile:
        for line in mapfile:
            line = line.split()

            uniparc = line[0]
            category = line[1]
            category_value = line[2]

            if category == args.category:
                uniparc2uniref[uniparc] = category_value
            elif category == "GeneID":
                if uniparc not in uniparc2uniref:
                    missing_uniref.add(uniparc)
                    continue

                geneid2uniref[category_value].add(uniparc2uniref[uniparc])
            else:
                continue

    # Double-check that all UniParc ids that did not seem to have a UniRef link
    # still are not present in the uniparc2uniref dictionary (i.e. check that
    # this mapping was not simply after the GeneID mapping).
    for uniparc in missing_uniref:
        if uniparc in uniparc2uniref:
            sys.exit("Error - UniParc id " + uniparc + " is now in the dict, "
                     "which means that the mapping to UniRef occurred after "
                     "the GeneID mapping.")

    # Write out dictionary as pickle file.
    with open(args.output, 'wb') as pickle_out:
        pickle.dump(geneid2uniref, pickle_out, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
