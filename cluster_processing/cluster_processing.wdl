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

        #Int atac_size_1 = floor(size(files.left.right, "GB"))
        #Int rna_size_1 = floor(size(files.left.left, "GB"))
        Int atac_size_1 = 2
        Int rna_size_1 = 2

        call get_cluster_data {
            input:
                expression_h5 = files.left.left,
                atac_fragments_tsv = files.left.right,
                input_name = files.right,
                cluster_labels = cluster_labels,
                docker_image = docker_image,
                git_branch = git_branch,
                atac_size = atac_size_1,
                rna_size = rna_size_1,
                min_num_fragments = min_num_fragments
        }
    }

    Array[String] all_rna_files = flatten(get_cluster_data.rna_h5ads)
    Array[String] all_atac_files = flatten(get_cluster_data.atac_h5ads)

    call get_cluster_file_map {
        input:
            all_rna_files = all_rna_files,
            all_atac_files = all_atac_files,
            cluster_names = select_first(get_cluster_data.cluster_names),
            git_branch = git_branch,
            docker_image = docker_image
    }

    Array[String] cluster_names = read_lines(select_first(get_cluster_data.cluster_names))

    scatter (cluster_name in cluster_names) {

        Int atac_size_2 = floor(size(get_cluster_file_map.atac_files[cluster_name], "GB"))
        Int rna_size_2 = floor(size(get_cluster_file_map.rna_files[cluster_name], "GB"))

        call concatenate_cluster {
            input:
                cluster_name = cluster_name,
                rna_files = get_cluster_file_map.rna_files[cluster_name],
                atac_files = get_cluster_file_map.atac_files[cluster_name],
                cluster_labels = cluster_labels,
                atac_size = atac_size_2,
                rna_size = rna_size_2,
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

    parameter_meta {
        expression_h5: {
            localization_optional: true
        }
        atac_fragments_tsv: {
            localization_optional: true
        }
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        /tmp/monitor_script.sh &
        micromamba run -n tools2 python3 /app/cluster_processing/fast_fail_get_clusters.py -c ${cluster_labels}
        #micromamba run -n tools2 python3 /app/cluster_processing/get_cluster_h5ads.py -e ${expression_h5} -f ${atac_fragments_tsv} -c ${cluster_labels} \
        #    -i ${input_name} ~{"--min_num_fragments=" + min_num_fragments}
    }

    output {
        Array[File] fragment_files = glob("*atac_fragments_clustered_*.tsv.gz")
        Array[File]+ rna_h5ads = glob("*rna_fake_file.txt")
        Array[File]+ atac_h5ads = glob("*atac_fake_file.txt")
        File cluster_names = "all_unique_clusters.txt"
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "128GB"
        preemptible: 1
        disks: "local-disk ~{disk_size} HDD"
    }
}

task get_cluster_file_map {
    input {
        Array[String] all_rna_files
        Array[String] all_atac_files
        File cluster_names
        String git_branch
        String docker_image
    }

    command {
        set -ex
        (git clone https://github.com/broadinstitute/hgrm_multiome_cluster_processing.git /app ; cd /app ; git checkout ${git_branch})
        /tmp/monitor_script.sh &
        micromamba run -n tools2 python3 /app/cluster_processing/get_cluster_file_map.py -c ${cluster_names} -r ${sep=' ' all_rna_files} \
            -a ${sep=' ' all_atac_files}
    }

    output {
        Map[String, Array[String]] rna_files = read_json("rna_files.json")
        Map[String, Array[String]] atac_files = read_json("atac_files.json")
    }

    runtime {
        docker: docker_image
    }
}

task save_map {
    input {
        Map[String, Array[String]] rna_files
        Map[String, Array[String]] atac_files
        String docker_image
    }

    command <<<
    >>>

    output {
        File rna_map_file_out = write_json(rna_files)
        File atac_map_file_out = write_json(atac_files)
    }

    runtime {
        docker: docker_image
    }
}

task concatenate_cluster {
    input {
        String cluster_name
        Array[File]+ rna_files
        Array[File]+ atac_files
        String cluster_labels
        Int atac_size
        Int rna_size
        String git_branch
        String docker_image
    }

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