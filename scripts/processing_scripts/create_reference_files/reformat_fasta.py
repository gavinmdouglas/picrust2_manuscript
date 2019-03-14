#!/usr/bin/python3

import argparse
import re
from gff3s_to_func_table import file2set


def main():

    parser = argparse.ArgumentParser(description="Reads in FASTA file and \
simplify header lines to be first element after splitting by specified \
delimiter (\"|\" by default). You can also optionally specify a different \
field. Sequences that are centroids of clusters as \
identified by having a size > 1 in the headerline will have a user-specified \
string appended to the name (\"_cluster\" by default). Instead of parsing the \
headerlines for cluster sizes users also have the option of inputting a file \
containing the ids that should have the cluster label appended to them. If only \
a subset of sequences are wanted then an optional file containing the ids of \
these sequences (after all other processing steps), one per line, can be input",
epilog='''Usage example:

python3 reformat_fasta.py -d \";\" -c \"_cluster\" input.fasta > output.fasta

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to FASTA file", required=True)

    parser.add_argument("-d", "--delimiter", metavar="STRING", type=str,
                        default="|", required=False,
                        help="Delimiter to split headerline")

    parser.add_argument("-c", "--centroid", metavar="STRING", type=str,
                        default="_cluster", required=False,
                        help="String to append to end of centroid names")

    parser.add_argument("--remove_gaps", action="store_true",
                        help="Flag to indicate that gaps characters should \
be removed.")

    parser.add_argument("--ids2keep", metavar="FILE", type=str,
                        default=None, required=False,
                        help="Optional file containing ids of sequences to \
retain. These ids should be in the post-process format (i.e. after splitting \
by the delimiter and appending the cluster label if approporiate). There \
should be one id per line.")

    parser.add_argument("--ids2label", metavar="FILE", type=str,
                        default=None, required=False,
                        help="Optional file containing ids of sequences to \
re-label with the cluster label (--centroid option) appended. There should \
be one id per line. The ids should be in the post-process format (see \
--ids2keep option description). Note that headerlines with a size > 1 will \
still be appended with the cluster label when this option is set.")

    parser.add_argument("--field_index", metavar="INT", type=int,
                        default=0, required=False,
                        help="Index of field containing field after "
                             "splitting by specified delimiter. Note this is "
                             "zero-based")

    args = parser.parse_args()

    # Read in ids to keep and ids to label as a cluster if each respective
    # is set. Also define boolean variable to keep track of whether a \
    # sequence should be printed out or not (this will always be True) \
    # if the --ids2keep option isn't set.
    if args.ids2keep:
        ids2keep = file2set(args.ids2keep)
        keep_seq = False
    else:
        keep_seq = True

    if args.ids2label:
        ids2label = file2set(args.ids2label)


    # Read in FASTA line-by-line.
    with open(args.fasta, "r") as fasta:

        for line in fasta:
            
            # Remove newline characters from end.
            line = line.rstrip("\n\r")

            # If not header line then print out and go to next loop
            # (unless gaps need to be removed).
            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] != ">":

                # Skip line if not in specified set of ids to retain.
                if not keep_seq:
                    continue

                # Remove gap characters if flag set.
                if args.remove_gaps:
                    line = line.replace("-", "")
                    # Print out line after removing gaps if non-empty.
                    if line:
                        print(line)

                else:
                    print(line)

                continue

            # If the number of seqs is in line, then figure out if this
            # is a centroid or not. Otherwise just print out simplified
            # header.

            new_header = line.split(args.delimiter)[args.field_index]

            if new_header[0] != ">":
                new_header = ">" + new_header

            if "size=" in line:

                size = int(re.match(r'.+;size=(\d+);', line).group(1))

                if size > 1:
                    new_header = new_header + args.centroid

            # Alternatively re-label id if it is listed as a cluster id.
            elif args.ids2label and new_header[1:] in ids2label:
                new_header = new_header + args.centroid

            # If ids2keep option set then check if id is in set.
            # Otherwise print out new header automatically.
            if args.ids2keep:

                if new_header[1:] in ids2keep:
                    keep_seq = True
                else:
                    keep_seq = False
                    continue

            print(new_header)


if __name__ == '__main__':
    main()
