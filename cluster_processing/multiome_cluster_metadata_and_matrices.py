import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as snap
import argparse


def main():
    # Read in args
    parser = argparse.ArgumentParser()
    parser.add_argument("expression_h5", type=str)
    parser.add_argument("fragments_tsv", type=str)
    parser.add_argument("cluster_labels", type=str)
    parser.add_argument("--min_num_fragments", type=int,
                        help='optional argument for min_fragments for reading in fragments')
    args = parser.parse_args()

    print("Reading in scRNA counts.")
    rna_counts = sc.read_10x_h5(args.expression_h5)

    print("Reading in scATAC fragments. This may take awhile.")
    atac_counts = snap.pp.import_data(
        args.fragments_tsv, chrom_sizes=snap.genome.hg38, sorted_by_barcode=False,
        min_num_fragments=200 if args.min_num_fragments is None else args.min_num_fragments
    )

    # cluster labels must have cell-barcodes as index, and cluster labels in first column (and likely only column)
    cell_clusters = pd.read_table(args.cluster_labels, index_col=0)
    # for writing out to big_wigs, need in string format
    # this matches the file type evie gave, with index = cell names,
    # and first (and should be only) column corresponding to cluster names, regardless of col name
    if (cell_clusters.iloc[:, 0].dtypes == "int" or cell_clusters.iloc[:, 0].dtypes == "float"):
        cell_clusters["cluster_string"] = "cluster_" + cell_clusters.iloc[:, 0].astype(str)
    else:
        cell_clusters.rename(
            columns={cell_clusters.columns[0]: "cluster_string"}, inplace=True
        )

    # Add cell cluster info to each file
    clusters_series = cell_clusters["cluster_string"]
    rna_counts.obs["CellClusterID"] = clusters_series
    atac_counts.obs["CellClusterID"] = clusters_series

    # drop cells with na for cluster
    atac_counts = atac_counts[atac_counts.obs.dropna().index, :]
    rna_counts = rna_counts[rna_counts.obs.dropna().index, :]

    print("Generate cell qc metrics.")
    # mito and ribo info, and total counts, maybe we want users to do this beforehand
    # but we at least need "total_counts" and "n_fragment" columns in rna, atac, respectively
    rna_counts.var["mito"] = rna_counts.var_names.str.startswith("MT-")
    rna_counts.var["ribo"] = rna_counts.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        rna_counts,
        qc_vars=["mito", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    rna_counts.obs.drop(
        columns=["total_counts_mito", "total_counts_ribo"], inplace=True
    )

    # copy objects metadata
    rna_obs = rna_counts.obs.copy()
    atac_obs = atac_counts.obs.copy()
    rna_obs.columns = "RNA_" + rna_obs.columns
    atac_obs.columns = "ATAC_" + atac_obs.columns
    # Info by cell
    per_cell_metadata = rna_obs.merge(
        atac_obs, left_index=True, right_index=True, how="outer"
    )
    # we use the original assignments, because some may not be present in both objects
    per_cell_metadata['CellClusterID'] = clusters_series

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

    # Get CPM information, normalized within each cluster group
    for cluster_name in rna_counts.obs.CellClusterID.unique():
        tot = rna_counts[rna_counts.obs.CellClusterID == cluster_name, :].X.A.sum(
            axis=0
        )
        cpm = (tot / tot.sum()) * 1e6
        rna_counts.var[f"{cluster_name}_CPM"] = cpm

    # Run these loops separately to be safe. some clusters may not be in both.
    print("Saving per-cluster atac files.")
    for cluster_name in atac_counts.obs.CellClusterID.unique():
        atac_cluster_only = atac_counts[
            atac_counts.obs["CellClusterID"] == cluster_name
        ]
        atac_cluster_only.write(
            (
                f"atac_{cluster_name}.h5ad"
                if "cluster" in cluster_name
                else f"atac_cluster_{cluster_name}.h5ad"
            ),
            compression="gzip",
        )

    print("Saving per-cluster rna files.")
    for cluster_name in rna_counts.obs.CellClusterID.dropna().unique():
        rna_cluster_only = rna_counts[rna_counts.obs["CellClusterID"] == cluster_name]
        rna_cluster_only.write_h5ad(
            (
                f"rna_{cluster_name}.h5ad"
                if "cluster" in cluster_name
                else f"rna_cluster_{cluster_name}.h5ad"
            ),
            compression="gzip",
        )

    # save metadata tables
    print("Saving metadata tables.")
    per_cell_metadata.to_csv("barcode_level_metadata.tsv", sep="\t", header=True)
    per_cluster_metadata.to_csv("cluster_level_metadata.tsv", sep="\t", header=True)

    print("Saving off cluster fragment files.")
    atac_counts.obs['CellClusterID'] = atac_counts.obs['CellClusterID'].astype(str)
    snap.ex.export_fragments(atac_counts, groupby='CellClusterID', suffix='.tsv.gz', prefix='atac_fragments_clustered_')

    print("Done.")


if __name__ == "__main__":
    main()
