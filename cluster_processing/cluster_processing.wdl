version development

workflow multiome_cluster_processing {
    input {
        Array[File]+ expression_h5s
        Array[File]+ atac_fragments_tsvs
        Array[String]+ input_names
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

    Map[String, Array[File]] rna_files = collect_by_key(flatten(get_cluster_data.rna_files))
    Map[String, Array[File]] atac_files = collect_by_key(flatten(get_cluster_data.rna_files))
    Array[String] cluster_names = keys(rna_files)

    scatter (cluster_name in cluster_names) {
        call concatenate_cluster {
            input:
                cluster_name = cluster_name,
                rna_files = rna_files[cluster_name],
                atac_files = atac_files[cluster_name],
                cluster_labels = cluster_labels,
                git_branch = git_branch,
                docker_image = docker_image
        }
    }

    call extract_metadata_tables {
        input:
            cluster_cell_metadatas = concatenate_cluster.cluster_cell_metadata,
            git_branch = git_branch,
            docker_image = docker_image
    }

    output {
        Array[File] rna_h5ads = concatenate_cluster.concatenated_rna_h5ad
        Array[File] atac_h5ads = concatenate_cluster.concatenated_atac_h5ad
        File barcode_level_metadata = extract_metadata_tables.barcode_level_metadata
        File cluster_level_metadata = extract_metadata_tables.cluster_level_metadata
        Array[File] all_fragment_files = flatten(get_cluster_data.fragment_files)
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
        micromamba run -n tools2 python3 /app/cluster_processing/get_cluster_h5ads.py -e ${expression_h5} -f ${atac_fragments_tsv} -c ${cluster_labels} \
            -i ${input_name} ~{"--min_num_fragments=" + min_num_fragments}
    }

    output {
        Array[Pair[String, File]] rna_files = as_pairs(read_map("rna_cluster_pairs.tsv"))
        Array[Pair[String, File]] atac_files = as_pairs(read_map("atac_cluster_pairs.tsv"))
        Array[File] fragment_files = glob("*atac_fragments_clustered_*.tsv.gz")
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "128GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}

task concatenate_cluster {
    input {
        String cluster_name
        Array[File]+ rna_files
        Array[File]+ atac_files
        String cluster_labels
        String git_branch
        String docker_image
    }

    Int atac_size = floor(size(atac_files, "GB"))
    Int rna_size = floor(size(rna_files, "GB"))
    Int disk_size = 100 + (3 * (atac_size + rna_size))

    command {
        echo "writing concatenated file and metadata"
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        /tmp/monitor_script.sh &
        micromamba run -n tools2 python3 /app/cluster_processing/get_subcluster_cell_metadata.py -s ${cluster_name} -c ${cluster_labels} \
            -r ${sep=' ' rna_files} \
            -a ${sep=' ' atac_files}
    }

    output {
        File concatenated_rna_h5ad = "rna_~{cluster_name}.h5ad"
        File concatenated_atac_h5ad = "atac_~{cluster_name}.h5ad"
        File cluster_cell_metadata = "~{cluster_name}_per_cell_metadata.txt"
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "128GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}

task extract_metadata_tables {
    input {
        Array[File] cluster_cell_metadatas
        String git_branch
        String docker_image
    }

    command {
        echo "joining all cell metadata across clusters, and getting cluster level data"
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        /tmp/monitor_script.sh &
        micromamba run -n tools2 python3 /app/cluster_processing/join_cluster_metadata.py -c ${sep=' ' cluster_cell_metadatas}
    }

    output {
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