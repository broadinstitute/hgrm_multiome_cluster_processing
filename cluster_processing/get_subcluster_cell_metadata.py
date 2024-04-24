import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import argparse
import snapatac2 as snap

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "subcluster_name", type=str, required=True)
    parser.add_argument("-c", "cluster_labels", type=str, required=True)
    parser.add_argument("-r", "--rna_h5ad", dest="clustered_rna_h5ads", nargs='+', default=[],
                        help="RNA h5ads for this cluster across all inputs")
    parser.add_argument("-a", "--atac_h5ad", dest="clustered_atc_h5ads", nargs='+', default=[],
                        help="ATAC h5ads for this cluster across all inputs")
    args = parser.parse_args()

    cluster_name = args.subcluster_name
    print(f"Reading in all {cluster_name} RNA files.")
    if args.clustered_rna_h5ad is not None:
        all_rna_counts = [sc.read_h5ad(rna) for rna in args.clustered_rna_h5ads]
        rna_counts = ad.concat(all_rna_counts)
        print(f"Writing out concatenated rna h5ad for cluster {cluster_name} across inputs.")
        rna_counts.write_h5ad(f'rna_{cluster_name}.h5ad')

    print(f"Reading in all {cluster_name} ATAC files.")
    if args.clustered_atac_h5ad is not None:
        all_atac_counts = [sc.read_h5ad(atac) for atac in args.clustered_rna_h5ads]
        atac_counts = ad.concat(all_atac_counts)
        print(f"Writing out concatenated atac h5ad for cluster {cluster_name} across inputs.")
        atac_counts.write_h5ad(f'atac_{cluster_name}.h5ad')


    # read in cell cluster info and get in same format as cluster names
    cell_clusters = pd.read_table(args.cluster_labels, index_col=0)
    if (cell_clusters.iloc[:, 0].dtypes == "int" or cell_clusters.iloc[:, 0].dtypes == "float"):
        cell_clusters["cluster_string"] = "cluster_" + cell_clusters.iloc[:, 0].astype(str)
    else:
        cell_clusters.rename(
            columns={cell_clusters.columns[0]: "cluster_string"}, inplace=True
        )

    #TODO: deal with this, write out empty files perhaps
    # account for rare cases where an entire cluster exists only in one of the data types (RNA or ATAC)
    if args.clustered_rna_h5ad is None:
        # make an RNA object with same dimension, filled with NaN for consistency later
        rna_obs = pd.DataFrame(index=atac_counts.obs_names,
                               columns=['CellClusterID', 'n_genes_by_counts', 'total_counts', 'pct_counts_mito', 'pct_counts_ribo'])
    else:
        # copy objects metadata
        rna_obs = rna_counts.obs.copy()
    if args.clustered_atac_h5ad is None:
        atac_obs = pd.DataFrame(index=rna_counts.obs_names,
                                columns=['n_fragment', 'frac_dup', 'frac_mito', 'CellClusterID'])
    else:
        # copy objects metadata
        atac_obs = atac_counts.obs.copy()

    # add type to col names
    rna_obs.columns = "RNA_" + rna_obs.columns
    atac_obs.columns = "ATAC_" + atac_obs.columns
    # Info by cell
    per_cell_metadata = rna_obs.merge(
        atac_obs, left_index=True, right_index=True, how="outer"
    )
    # we use the original assignments, because some may not be present in both objects
    per_cell_metadata['CellClusterID'] = cell_clusters.cluster_string

    per_cell_metadata.to_csv(f'{cluster_name}_per_cell_metadata.txt', sep='\t')



if __name__ == "__main__":
    main()
