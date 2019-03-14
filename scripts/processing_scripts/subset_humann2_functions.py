#!/usr/bin/python3

import argparse

def main():

    parser = argparse.ArgumentParser(

        description="Subset unstratified humann2 table to only those functions in list of functions.",

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--ids2keep", metavar="PATH", type=str,
                        help="Path to file with ids 2 keep", required=True)

    parser.add_argument("-i", "--input", metavar="PATH", type=str,
                        help="Path to unstratified humann2 table",
                        required=True)


    args = parser.parse_args()

    ids2keep = set()

    with open(args.ids2keep, 'r') as ids_infile:
        for line in ids_infile:
            line = line.rstrip()
            ids2keep.add(line)

    lc = 0
    with open(args.input, 'r') as input_humann2:
        for line in input_humann2:
            line = line.rstrip()
            if lc == 0:
                print(line)
            lc += 1

            line_split = line.split(' ')
            line_split[0] = line_split[0].replace(':', '')
            if line_split[0] in ids2keep:
                print(line)

if __name__ == '__main__':
    main()
