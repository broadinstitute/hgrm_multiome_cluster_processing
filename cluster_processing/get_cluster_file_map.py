import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rna_h5ad", dest="rna_filenames", nargs='+', default=[],
                        help="RNA h5ads file strings")
    parser.add_argument("-a", "--atac_h5ad", dest="atac_filenames", nargs='+', default=[],
                        help="ATAC h5ads file strings")
    args = parser.parse_args()

    rna_clusters = [(fname.split('/')[-2], fname) for fname in args.rna_filenames]
    atac_clusters = [(fname.split('/')[-2], fname) for fname in args.rna_filenames]

    with open('rna_cluster_pairs.tsv', 'w') as f:
        for i in range(len(rna_clusters)):
            f.write(f"{rna_clusters[i]}\t{args.rna_filenames[i]}\n")

    with open('atac_cluster_pairs.tsv', 'w') as f:
        for i in range(len(atac_clusters)):
            f.write(f"{atac_clusters[i]}\t{args.atac_filenames[i]}\n")

if __name__ == "__main__":
    main()