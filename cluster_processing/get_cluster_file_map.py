import argparse
import pandas as pd
import numpy as np
import json
from pathlib import Path
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="cluster_labels", type=str, required=True)
    parser.add_argument("-r", "--rna_h5ad", dest="rna_filenames", nargs='+', default=[],
                        help="RNA h5ads file strings")
    parser.add_argument("-a", "--atac_fragments", dest="atac_filenames", nargs='+', default=[],
                        help="ATAC fragment file strings")
    args = parser.parse_args()

    cluster_names = [ln.strip() for ln in open(args.cluster_labels)]

    rna_files = defaultdict(list)
    for filename in args.rna_filenames:
        for cluster in cluster_names:
            if filename.endswith(f'/{cluster}_rna.h5ad'):
                rna_files[cluster].append(filename)

    for cluster in cluster_names:
        if cluster not in rna_files:
            rna_files[cluster] = []

    json.dump(rna_files, open("rna_files.json", "w"))


    atac_files = defaultdict(list)
    for filename in args.atac_filenames:
        for cluster in cluster_names:
            if filename.endswith(f'atac_fragments_clustered_{cluster}.tsv.gz'):
                atac_files[cluster].append(filename)

    for cluster in cluster_names:
        if cluster not in atac_files:
            atac_files[cluster] = []

    json.dump(atac_files, open("atac_files.json", "w"))


if __name__ == "__main__":
    main()