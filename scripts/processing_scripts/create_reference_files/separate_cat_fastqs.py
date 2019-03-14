#!/usr/bin/python3

import argparse
import gzip
import sys
import os


def main():

    parser = argparse.ArgumentParser(
        description="Splits GZIPPED FASTQ of R1 and R2 sequences into separate FASTQs",

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="PATH",
                        type=str, help="Path to gzipped FASTQ.",
                        required=True)

    parser.add_argument("-o", "--output", metavar="PATH",
                        type=str, help="Folder to write output files.",
                        required=True)

    args = parser.parse_args()

    # Make output folder if it does not exist.
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Determine output file names and open output filehandles.
    out_prefix = os.path.basename(args.input).replace(".fastq.gz", "")
    R1_outfile = os.path.join(args.output, out_prefix + "_R1.fastq.gz")
    R2_outfile = os.path.join(args.output, out_prefix + "_R2.fastq.gz")

    R1_out = gzip.open(R1_outfile, 'wb')
    R2_out = gzip.open(R2_outfile, 'wb')

    # Variable to keep track of whether last read was R1 or R2.
    read_type = None

    # Counter to determine when each 4th line is reached.
    line_i = 0

    # Read through FASTQ, determine whether read is R1 or R2, and write to
    # appropriate output FASTQ.
    with gzip.open(args.input, 'rb') as fastq_in:
        for line in fastq_in:

            line = line.rstrip()

            # Check if this is first line of read.
            if line_i == 0:
                if line[-2:] == "/1":
                    read_type = "R1"
                elif line[-2:] == "/2":
                    read_type = "R2"
                else:
                    sys.exit("Error - read type not in FASTQ line: " + line)
            line_i += 1

            # Print line out to appropriate outfile.
            if read_type == "R1":
                R1_out.write(line + "\n")
            elif read_type == "R2":
                R2_out.write(line + "\n")

            # Reset counter if last line of read.
            if line_i == 4:
                line_i = 0

    R1_out.close()
    R2_out.close()

if __name__ == '__main__':
    main()
