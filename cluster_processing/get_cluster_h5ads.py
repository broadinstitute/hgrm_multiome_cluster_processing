import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import snapatac2 as snap
import argparse


def main():
    # Read in args
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", dest="expression_h5", type=str, required=True)
    parser.add_argument("-f", dest="fragments_tsv", type=str, required=True)
    parser.add_argument("-c", dest="cluster_labels", type=str, required=True)
    parser.add_argument("-i", dest="input_name", type=str, nargs='?', const='')
    parser.add_argument("--min_num_fragments", type=int,
                        help='optional argument for min_fragments for reading in fragments')
    args = parser.parse_args()

    input_name = args.input_name
    print("Reading in scRNA counts.")
    rna_counts = sc.read_10x_h5(args.expression_h5)

    print("Reading in scATAC fragments. This may take awhile.")
    atac_counts = snap.pp.import_data(
        args.fragments_tsv, chrom_sizes=snap.genome.hg38, sorted_by_barcode=False,
        min_num_fragments=200 if args.min_num_fragments is None else args.min_num_fragments
    )

    print("Reading in cell cluster assignments.")
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

    # Get CPM information, normalized within each cluster group
    for cluster_name in rna_counts.obs.CellClusterID.unique():
        tot = rna_counts[rna_counts.obs.CellClusterID == cluster_name, :].X.A.sum(
            axis=0
        )
        cpm = (tot / tot.sum()) * 1e6
        rna_counts.var[f"{cluster_name}_CPM"] = cpm

    # Run these loops separately to be safe. some clusters may not be in both.
    print("Saving per-cluster atac files.")
    with open('atac_cluster_pairs.tsv', 'w') as f:
        for cluster_name in atac_counts.obs.CellClusterID.unique():
            atac_cluster_only = atac_counts[
                atac_counts.obs["CellClusterID"] == cluster_name
            ]
            fname = f"{input_name}_atac_{cluster_name}.h5ad" if "cluster" in cluster_name else f"{input_name}_atac_cluster_{cluster_name}.h5ad"
            atac_cluster_only.write(fname,
                compression="gzip",
            )
            cluster_str = cluster_name if "cluster" in cluster_name else f"cluster_{cluster_name}"
            f.write(f"{cluster_str}\t{fname}\n")


    print("Saving per-cluster rna files.")
    with open('rna_cluster_pairs.tsv', 'w') as f:
        for cluster_name in rna_counts.obs.CellClusterID.dropna().unique():
            rna_cluster_only = rna_counts[
                rna_counts.obs["CellClusterID"] == cluster_name
            ]
            fname = f"{input_name}_rna_{cluster_name}.h5ad" if "cluster" in cluster_name else f"{input_name}_rna_cluster_{cluster_name}.h5ad"
            rna_cluster_only.write_h5ad(fname,
                compression="gzip",
            )
            cluster_str = cluster_name if "cluster" in cluster_name else f"cluster_{cluster_name}"
            f.write(f"{cluster_str}\t{fname}\n")

    print("Saving off cluster fragment files.")
    atac_counts.obs['CellClusterID'] = atac_counts.obs['CellClusterID'].astype(str)
    snap.ex.export_fragments(atac_counts, groupby='CellClusterID', suffix='.tsv.gz', prefix=f'{input_name}_atac_fragments_clustered_')

    print("Done.")


if __name__ == "__main__":
    main()
