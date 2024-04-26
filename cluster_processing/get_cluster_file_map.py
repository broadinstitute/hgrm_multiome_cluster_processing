import argparse
import pandas as pd
import numpy as np
import json
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="cluster_labels", type=str, required=True)
    parser.add_argument("-r", "--rna_h5ad", dest="rna_filenames", nargs='+', default=[],
                        help="RNA h5ads file strings")
    parser.add_argument("-a", "--atac_h5ad", dest="atac_filenames", nargs='+', default=[],
                        help="ATAC h5ads file strings")
    args = parser.parse_args()

    cluster_names = [ln.strip() for ln in open(args.cluster_labels)]

    rna_files = defaultdict(list)
    for filename in args.rna_filenames:
        for cluster in cluster_names:
            if filename.endswith(f'/{cluster}_rna_fake_file.txt'):
                rna_files[cluster].append(filename)
    json.dump(rna_files, open("rna_files.json", "w"))


    atac_files = defaultdict(list)
    for filename in args.atac_filenames:
        for cluster in cluster_names:
            if filename.endswith(f'/{cluster}_atac_fake_file.txt'):
                rna_files[cluster].append(filename)
    json.dump(atac_files, open("atac_files.json", "w"))


if __name__ == "__main__":
    main()