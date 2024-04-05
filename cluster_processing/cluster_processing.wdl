version 1.0

workflow multiome_cluster_processing {
    input {
        File expression_h5
        File atac_fragments_tsv
        File cluster_labels
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.2'
        String git_branch = "main"
    }

    call get_cluster_data {
        input:
            expression_h5 = expression_h5,
            atac_fragments_tsv = atac_fragments_tsv,
            cluster_labels = cluster_labels,
            docker_image = docker_image,
            git_branch = git_branch
    }

    output {
        Array[File] h5ads = get_cluster_data.h5ads
        File barcode_level_metadata = get_cluster_data.barcode_level_metadata
        File cluster_level_metadata = get_cluster_data.cluster_level_metadata
    }
}

task get_cluster_data {
    input {
        File expression_h5
        File atac_fragments_tsv
        File cluster_labels
        String docker_image
        String git_branch
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        micromamba run -n tools2 python3 /app/cluster_processing/multiome_cluster_metadata_and_matrices.py ${expression_h5} ${atac_fragments_tsv} ${cluster_labels}
    }

    output {
        Array[File]+ h5ads = glob("*.h5ad")
        File barcode_level_metadata = "barcode_level_metadata.tsv"
        File cluster_level_metadata = "cluster_level_metadata.tsv"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "16GB"
        preemptible: 1
    }
}