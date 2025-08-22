import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", dest="barcode_metadata", type=str, required=True)
    parser.add_argument("-r", "--rna_h5ad", dest="rna_h5ads", nargs='+', default=[],
                        help="RNA h5ads file strings")
    parser.add_argument("-a", "--atac_h5ad", dest="atac_h5ads", nargs='+', default=[],
                        help="ATAC h5ads file strings")
    parser.add_argument("-f", "--fragments", dest="atac_fragments", nargs='+', default=[],
                        help="ATAC fragment file strings")
    args = parser.parse_args()

    barcode_metadata = pd.read_table(args.barcode_metadata, index_col=0)

    cluster_df = pd.DataFrame({'atac_h5ads':args.atac_h5ads, 'rna_h5ads':args.rna_h5ads})
    cluster_df['cluster_id'] = cluster_df.atac_h5ads.str.split('/').str[-1].str.split('.').str[0].str.split('atac_').str[-1]
    assert (cluster_df.rna_h5ads.str.split('/').str[-1].str.split('.').str[0].str.split('rna_').str[-1] == cluster_df['cluster_id']).all()

    cluster_df.set_index('cluster_id', inplace=True)
    full_frag_file_list = []
    for cluster_id in cluster_df.index:
        cluster_list = []
        for frag in args.atac_fragments:
            if f'{cluster_id}.' in frag:
                cluster_list.append(frag)
        full_frag_file_list.append(cluster_list)

    cluster_df['atac_fragments'] = full_frag_file_list