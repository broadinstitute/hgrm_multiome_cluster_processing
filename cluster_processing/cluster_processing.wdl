version 1.0

workflow multiome_cluster_processing {
    input {
        Array[File]+ expression_h5s
        Array[File]+ atac_fragments_tsvs
        Array[File]+ input_names
        File cluster_labels
        Int? min_num_fragments
        String docker_image = 'us.gcr.io/landerlab-atacseq-200218/hgrm_multiome_cluster_processing:0.6'
        String git_branch = "main"
    }

    scatter (files in zip(zip(expression_h5s, atac_fragments_tsvs), input_names)) {

        Int atac_size = floor(size(files.left.right, "GB"))
        Int rna_size = floor(size(files.left.left, "GB"))

        call get_cluster_data {
            input:
                expression_h5 = files.left.left,
                atac_fragments_tsv = files.left.right,
                input_name = files.right,
                cluster_labels = cluster_labels,
                docker_image = docker_image,
                git_branch = git_branch,
                atac_size = atac_size,
                rna_size = rna_size,
                min_num_fragments = min_num_fragments
        }
    }

    output {
        Array[File] rna_h5ads = get_cluster_data.rna_h5ads
        Array[File] atc_h5ads = get_cluster_data.atac_h5ads
        File barcode_level_metadata = get_cluster_data.barcode_level_metadata
        File cluster_level_metadata = get_cluster_data.cluster_level_metadata
        Array[File] fragment_files = get_cluster_data.fragment_files
    }
}

task get_cluster_data {
    input {
        File expression_h5
        File atac_fragments_tsv
        File cluster_labels
        String input_name
        String docker_image
        String git_branch
        Int atac_size
        Int rna_size
        Int? min_num_fragments
    }

    Int disk_size = 100 + (3 * (atac_size + rna_size))

    command {
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        /tmp/monitor_script.sh &
        micromamba run -n tools2 python3 /app/cluster_processing/get_cluster_h5ads.py -r ${expression_h5} -a ${atac_fragments_tsv} -c ${cluster_labels} \
            -n ${input_name} ~{"--min_num_fragments=" + min_num_fragments}
    }

    output {
        Array[File]+ rna_h5ads = glob("rna_*.h5ad")
        Array[File]+ atac_h5ads = glob("atac_*.h5ad")
        File barcode_level_metadata = "barcode_level_metadata.tsv"
        File cluster_level_metadata = "cluster_level_metadata.tsv"
        Array[File]+ fragment_files = glob("atac_fragments_clustered_*.tsv.gz")
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "128GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}