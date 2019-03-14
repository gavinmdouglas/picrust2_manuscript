import argparse
import numpy as np
import pandas as pd
from gff3s_to_func_table import file2set

def main():

    parser = argparse.ArgumentParser(
                description="Will calculate the mean trait values for all "
                            "rows in a cluster. Requires an input trait "
                            "table, a file containing clustered sequences, "
                            "and optionally a file containing which rows to "
                            "exclude. Mean values will be rounded UP to the "
                            "nearest integer. Any functions that are all 0s "
                            "will be removed.",
                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--func_tab", metavar="FUNC_TABLE",
                        type=str, help="Input table of function abundances - "
                                       "functions are colunns.", required=True)

    parser.add_argument("-o", "--output", metavar="OUT_TABLE",
                        type=str, help="Output table.", required=True)

    parser.add_argument("-c", "--clusters", metavar="CLUSTERS_FILE",
                        type=str, help="File containing links between "
                                       "sequences in a cluster (output of "
                                       "derep_fasta.py).", required=True)

    parser.add_argument("-e", "--cluster2exclude", metavar="EXCLUDE_FILE",
                        type=str, help="File containing clusters to exclude.",
                        required=False, default=None)

    parser.add_argument("--function2keep", metavar="FUNCTION_TO_KEEP",
                        type=str, required=False, default=None,
                        help="File containing functions that should be kept. "
                             "Expected format is 1 function id per line (and "
                             "should be the first field after delmiting by "
                             "whitespace.")

    args = parser.parse_args()

    rows2rm = []
    # Read in rows to rm if applicable.
    if args.cluster2exclude:
        rows2rm = list(file2set(args.cluster2exclude))

    if args.function2keep:
        raw_cols2keep = pd.read_table(args.function2keep, sep="\t", header=None)
        cols2keep = list(raw_cols2keep.iloc[:, 0])
        cols2keep.sort()

    in_table = pd.read_table(args.func_tab, sep="\t", index_col=0)

    print("Read in input functional table.")

    # Convert rownames to strings (in case they are all integers).
    in_table.index = [str(x) for x in in_table.index.values]

    # If exclude options set then drop the corresponding rows and columns
    if args.cluster2exclude:
        rows2rm = [row for row in rows2rm if row in in_table.index]
        in_table.drop(rows2rm, axis=0, inplace=True)

    print('Rows excluded.')

    if args.function2keep:
        in_table = in_table[cols2keep]

    print('Specified columns retained.')

    # Dict with cluster name as key and seqs as values (only for clusters
    # with more than 1 sequence.
    clust = {}

    # List of singleton sequences.
    singletons = []

    with open(args.clusters, "r") as clusters_in:
        for line in clusters_in:
            line = line.rstrip()
            line_split = line.split(" ")

            if line_split[0] in rows2rm:
                continue
            
            # If multiple seqs in cluster then add -"cluster" to cluster name.
            if len(line_split) > 1:
                cluster_name = line_split[0] + "-cluster"
                clust[cluster_name] = line_split
            else:
                singletons += [line_split[0]]

    print('Read through clusters file.')

    # Subset rows that correspond to singletons (i.e. no need to take mean).
    in_table_single = in_table.loc[singletons, :]

    # Subset different df to all rows EXCEPT singletons.
    in_table_clusters = in_table.drop(singletons, axis=0)

    # First take off 0.000000001 so that 0.5 values are always
    # rounded down.
    in_table_clusters = in_table_clusters - 0.000000001

    # Add column that will contain clusters to this table (and set all values
    # to NAN)
    in_table_clusters["clusters"] = np.nan

    # Loop through all cluster ids and add cluster id to all rows for that
    # cluster.
    for cluster_id, cluster_seqs in clust.items():
        in_table_clusters.loc[cluster_seqs, "clusters"] = cluster_id

    print('Made the two different dataframes for singletons and clusters.')

    # Get mean abundances for all rows in the same cluster. Then round each
    # value to the nearest integer.
    in_table_clusters = pd.pivot_table(in_table_clusters, index="clusters",
                                       aggfunc=np.mean).round()

    in_table_clusters = in_table_clusters.astype(int)

    print('Calculated mean for rows in same cluster.')

    out_table = in_table_single.append(in_table_clusters)

    print('Appended single genomes and clusters together.')

    cols2drop = (out_table == 0).all(axis=0)
    out_table = out_table.loc[:, (out_table != 0).any(axis=0)]

    print('Removed functions that are zero across all assemblies/clusters. The following columns were removed (if true):')

    print(cols2drop)

    # Append the dataframes of singletons and clusters together and write out.
    out_table.to_csv(args.output, sep="\t", index_label="assembly")


if __name__ == '__main__':
    main()
