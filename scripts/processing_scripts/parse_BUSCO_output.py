#!/usr/bin/python3

import argparse
import os

def main():

    parser = argparse.ArgumentParser(
        description="Script to parse BUSCO output folders and return one "
                    "table of all short summary tables.",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input_folder", metavar="PATH",
                        type=str, help="Path to folder with output folders.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="PATH",
                        type=str, help="Path to output table.",
                        required=True)

    args = parser.parse_args()

    in_folders = os.listdir(args.input_folder)

    # Open output table and write out header.
    out_table = open(args.output, "w")
    print("File\tC\tS\tD\tF\tM\tTotal", file=out_table)

    # Loop over all input folders and parse short summary files.
    for f in in_folders:
        basename = f.replace("run_", "short_summary_")
        summary_file = basename + ".txt"

        with open(os.path.join(args.input_folder, f, summary_file), "r") as in_summary:

            for line in in_summary:

                if "Complete BUSCOs (C)" in line:
                    C_val = line.split()[0]
                elif "Complete and single-copy BUSCOs (S)" in line:
                    S_val = line.split()[0]
                elif "Complete and duplicated BUSCOs (D)" in line:
                    D_val = line.split()[0]
                elif "Fragmented BUSCOs (F)" in line:
                    F_val = line.split()[0]
                elif "Missing BUSCOs (M)" in line:
                    M_val = line.split()[0]
                elif "Total BUSCO groups searched" in line:
                    T_val = line.split()[0]

        print(basename, C_val, S_val, D_val, F_val, M_val, T_val, sep="\t",
              file=out_table)

    out_table.close()


if __name__ == '__main__':
    main()
