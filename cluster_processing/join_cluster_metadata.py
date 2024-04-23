import numpy as np
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cell_metadata", dest="subcluster_cell_metadata", nargs='+', default=[], required=True,
                        help="All of the metadata files generated for each subcluster")
    args = parser.parse_args()

    subcluster_dfs = [pd.read_table(file, index_col=0) for file in args.subcluster_cell_metadata]
    per_cell_metadata = pd.concat(subcluster_dfs, axis=0)
    #Info by cluster
    # atac data needs to have n_fragment col in obs, coming in.
    per_cluster_metadata = (
        per_cell_metadata.reset_index()
        .groupby("CellClusterID")
        .agg(
            nCells=("index", "count"),
            MeanRNAUMIsPerCell=("RNA_total_counts", "mean"),
            MeanATACFragmentsPerCell=("ATAC_n_fragment", "mean"),
        )
    )

    # save metadata tables
    print("Saving total metadata tables.")
    per_cell_metadata.to_csv("barcode_level_metadata.tsv", sep="\t", header=True)
    per_cluster_metadata.to_csv("cluster_level_metadata.tsv", sep="\t", header=True)


if __name__ == "__main__":
    main()
