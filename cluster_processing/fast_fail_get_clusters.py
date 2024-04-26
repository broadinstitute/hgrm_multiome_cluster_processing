import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="cluster_labels", type=str, required=True)
    args = parser.parse_args()

    cell_clusters = pd.read_table(args.cluster_labels, index_col=0)
    cell_clusters_unique = pd.DataFrame(cell_clusters.iloc[:, 0].unique().astype(str), columns=['cluster'])
    cell_clusters_unique['rna_fake_file'] = 'rna_' + cell_clusters_unique['cluster'] + '_fake_file.txt'
    cell_clusters_unique['atac_fake_file'] = 'atac_' + cell_clusters_unique['cluster'] + '_fake_file.txt'

    for file in cell_clusters_unique.rna_fake_file:
        Path(f'{file}').touch()

    for file in cell_clusters_unique.atac_fake_file:
        Path(f'{file}').touch()

    cell_clusters_unique[['cluster', 'rna_fake_file']].to_csv('rna_cluster_pairs.tsv', sep='\t', header=False)
    cell_clusters_unique[['cluster', 'atac_fake_file']].to_csv('atac_cluster_pairs.tsv', sep='\t', header=False)

    Path("atac_fragments_clustered_fake.tsv.gz").touch()


if __name__ == "__main__":
    main()