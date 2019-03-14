#!/usr/bin/python3

import argparse


def read_fasta(filename, cut_header=False):

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

                name = name.rstrip("\r\n")

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    return seq


def main():

    parser = argparse.ArgumentParser(description="Reads in FASTA file and "
                                                 "identified identical "
                                                 "sequences. Will output file "
                                                 "listing identical sequences "
                                                 "on each line and a dereplicated "
                                                 "fasta as well. Will optionally "
                                                 "throw out sequences smaller "
                                                 "than a given length as well.",


formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        help="Path to FASTA file", required=True)

    parser.add_argument("--out_lists", metavar="OUTFILE", type=str,
                        help="Path to outfile of identical seqs",
                        required=True)

    parser.add_argument("--out_derep", metavar="OUTFILE", type=str,
                        help="Path to outfile of dereplicated seqs",
                        required=True)

    parser.add_argument("--min_length", metavar="INT", type=int, default=None,
                        help="Min length of seqs (after removing gaps and Ns)",
                        required=False)

    parser.add_argument("--out_small", metavar="OUTFILE", type=str,
                        help="Path to outfile of ids of sequences below min_length.",
                        required=False, default="short_seq_ids.txt")

    args = parser.parse_args()

    in_fasta = read_fasta(args.fasta)

    # Loop through all sequences and identify identical sequences.
    identical_seqs = {}
    small_seqs = set()

    if args.min_length:
        short_out = open(args.out_small, "w")

    for seq_id, seq in in_fasta.items():

        # Convert seq to uppercase.
        seq = seq.upper()

        # Remove gap and N characters.
        seq_nogap = seq.replace("-", "")
        seq_nogap = seq_nogap.replace("N", "")

        if args.min_length:
            if len(seq_nogap) < args.min_length:
                small_seqs.add(seq_id)
                print(seq_id, file=short_out)
                continue

        if seq_nogap in identical_seqs:
            identical_seqs[seq_nogap] += [seq_id]
        else:
            identical_seqs[seq_nogap] = [seq_id]

    if args.min_length:
        short_out.close()

    lists_outfile = open(args.out_lists, "w")
    derep_outfile = open(args.out_derep, "w")

    for seq, ids in identical_seqs.items():

        # Print out ids matching identical sequence to file.
        print(" ".join(ids), file=lists_outfile)

        if len(ids) > 1:
            print(">" + str(ids[0]) + "-cluster", file=derep_outfile)
        else:
            print(">" + str(ids[0]), file=derep_outfile)

        print(in_fasta[ids[0]], file=derep_outfile)


    lists_outfile.close()
    derep_outfile.close()


if __name__ == '__main__':
    main()
