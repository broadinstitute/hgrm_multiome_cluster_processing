import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="cluster_labels", type=str, required=True)
    parser.add_argument("-r", "--rna_h5ad", dest="rna_filenames", nargs='+', default=[],
                        help="RNA h5ads file strings")
    parser.add_argument("-a", "--atac_h5ad", dest="atac_filenames", nargs='+', default=[],
                        help="ATAC h5ads file strings")
    args = parser.parse_args()

    cluster_names = [ln.strip() for ln in open(args.cluster_labels)]

    with open('rna_cluster_pairs.tsv', 'w') as f:
        for filename in args.rna_filenames:
            for cluster in cluster_names:
                if filename.endswith(f'/{cluster}_rna_fake_file.txt'):
                    print(f"{cluster}\t{filename}", file=f)

    with open('atac_cluster_pairs.tsv', 'w') as f:
        for filename in args.atac_filenames:
            for cluster in cluster_names:
                if filename.endswith(f'/{cluster}_atac_fake_file.txt'):
                    print(f"{cluster}\t{filename}", file=f)

if __name__ == "__main__":
    main()