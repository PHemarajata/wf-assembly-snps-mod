#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CLUSTER_GENOMES {
    tag "clustering"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.{tsv,txt}"

    input:
    path mash_distances

    output:
    path "clusters.tsv", emit: clusters
    path "cluster_summary.txt", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threshold = params.mash_threshold ?: 0.03
    def max_cluster_size = params.max_cluster_size ?: 100
    """
    # Install required packages
    pip install pandas numpy scipy

    # Run clustering
    python3 ${projectDir}/bin/cluster_mash.py \
        $mash_distances \
        clusters.tsv \
        --threshold $threshold \
        --max-cluster-size $max_cluster_size \
        $args > cluster_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
}