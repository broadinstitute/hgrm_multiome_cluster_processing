version 1.0

workflow multiome_cluster_processing {
    input {
        File expression_h5
        File atac_fragments_tsv
        File cluster_labels
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        String git_branch = "main"
    }

    Int atac_size = floor(size(atac_fragments_tsv, "GB"))
    Int rna_size = floor(size(expression_h5, "GB"))

    call get_cluster_data {
        input:
            expression_h5 = expression_h5,
            atac_fragments_tsv = atac_fragments_tsv,
            cluster_labels = cluster_labels,
            docker_image = docker_image,
            git_branch = git_branch,
            atac_size = atac_size,
            rna_size = rna_size
    }

    output {
        Array[File] h5ads = get_cluster_data.h5ads
        File barcode_level_metadata = get_cluster_data.barcode_level_metadata
        File cluster_level_metadata = get_cluster_data.cluster_level_metadata
        Array[File] bigwigs = get_cluster_data.bigwigs

    }
}

task get_cluster_data {
    input {
        File expression_h5
        File atac_fragments_tsv
        File cluster_labels
        String docker_image
        String git_branch
        Int atac_size
        Int rna_size
    }

    Int disk_size = 100 + (3 * (atac_size + rna_size))

    command {
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        /tmp/monitor_script.sh &
        micromamba run -n tools2 python3 /app/cluster_processing/multiome_cluster_metadata_and_matrices.py ${expression_h5} ${atac_fragments_tsv} ${cluster_labels}
    }

    output {
        Array[File]+ h5ads = glob("*.h5ad")
        File barcode_level_metadata = "barcode_level_metadata.tsv"
        File cluster_level_metadata = "cluster_level_metadata.tsv"
        Array[File]+ bigwigs = glob("*.bw")
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "32GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}