#!/usr/bin/python3

import argparse
import os

def main():

    parser = argparse.ArgumentParser(
        description="Parse all metaxa2 taxonomy files in a folder and output "
                    "table of taxa counts.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input_folder", metavar="PATH",
                        type=str, help="Path to folder with output folders.",
                        required=True)

    args = parser.parse_args()

    infiles = os.listdir(args.input_folder)

    taxa_counts = {}

    unique_taxa = set()

    samples = []

    for f in infiles:
        if "taxonomy" not in f:
            continue

        sample = f.replace("_metaxa2_out.taxonomy.txt", "")

        taxa_counts[sample] = {}

        samples.append(sample)

        full_file = args.input_folder + "/" + f

        with open(full_file, 'r') as in_tax:
            for line in in_tax:
                line_split = line.split('\t')
                taxa_split = line_split[1].split(';')
                if len(taxa_split) == 0:
                    tax_string = "Unknown"
                elif taxa_split[0] == '':
                    tax_string = "Unknown"
                elif len(taxa_split) == 1 or taxa_split[1] == '':
                    tax_string = taxa_split[0] + '|' + "Unknown"
                else:
                    tax_string = taxa_split[0] + '|' + taxa_split[1]

                unique_taxa.add(tax_string)

                if tax_string in taxa_counts[sample]:
                    taxa_counts[sample][tax_string] = taxa_counts[sample][tax_string] + 1
                else:
                    taxa_counts[sample][tax_string] = 1

    unique_taxa_sorted = sorted(list(unique_taxa))
    headerline = "\t".join(["sample"] + unique_taxa_sorted)
    print(headerline)

    for sample in samples:

        out_line = [sample]

        for taxon in unique_taxa_sorted:

            if taxon in taxa_counts[sample]:
                out_line.append(str(taxa_counts[sample][taxon]))
            else:
                out_line.append(str(0))

        print("\t".join(out_line))

if __name__ == '__main__':
    main()
