#!/usr/bin/python3

import argparse
import skbio.sequence
import os
import sys

def read_fasta(filename, cut_header=True):

    '''Read in FASTA file and return dictionary with each independent sequence
    id as a key and the corresponding sequence string as the value.
    '''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    with open(filename, "r") as fasta:

        for line in fasta:

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":

                if cut_header:
                    name = line.split()[0][1:]
                else:
                    name = line[1:]

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    return seq


def main():

    parser = argparse.ArgumentParser(description="Reads in rRNA coordinates\
 from gff3 file and slices out sequence from assembly FASTA. User can specify\
 other rRNAs to output besides 16S if needed. Note that other coordinates can\
 be in the gff3 as well, but only elements annotated as \"rRNA\" will be\
 parsed", epilog='''Example of usage:

python3 rRNA_from_gff3.py -p 16s_product GFF3 REF_GENOME

 ''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("gff", metavar="annotation_file", type=str,
                        help="gff3 file containing rRNA coordinates")

    parser.add_argument("ref", metavar="reference_genome", type=str,
                        help="fasta file with reference sequence")

    parser.add_argument("-p", "--product_str", metavar="rRNA_product_string",
                        type=str, help="String matching product of rRNA line\
 (default=16s_rRNA)", default="16s_rRNA")

    args = parser.parse_args()

    ref = read_fasta(args.ref)

    # Get reference genome basename without extension.
    ref_name, ref_name_extension = os.path.splitext(
                                            os.path.basename(args.ref))

    rRNA_counter = 0

    # Read through gff3 line-by-line.
    with open(args.gff, "r") as gff3:
        for line in gff3:

            # If first character is comment then skip.
            if line[0] == "#":
                continue

            # Remove line terminator/newline character.
            line = line.rstrip("\r\n")

            # Split line on tab characters.
            line_split = line.split("\t")

            # Check if element is defined as rRNA and matches product string.
            if line_split[2] == "rRNA" and args.product_str in line_split[8]:

                rRNA_counter += 1

                scaffold = line_split[0]
                start = int(line_split[3])
                stop = int(line_split[4])
                strand = line_split[6]

                # Get nucleotide sequence of rRNA based on coordinates.
                rRNA_slice = skbio.sequence.DNA(sequence=ref[scaffold]
                                                [start-1:stop].upper())

                # Define name of sequence based on genome and genomic
                # coordinates.

                name = "|".join([ref_name, scaffold, str(start), str(stop),
                                 strand])

                # Take reverse complement if gene is on negative strand.
                if strand == "-":
                    rRNA_slice = rRNA_slice.reverse_complement()

                # Print rRNA sequence in FASTA format.
                print(">" + name)
                print(rRNA_slice)

    if rRNA_counter == 0:
        print("Warning: file " + args.gff +
              " did not contain any rRNA matches.", file=sys.stderr)


if __name__ == '__main__':
    main()
