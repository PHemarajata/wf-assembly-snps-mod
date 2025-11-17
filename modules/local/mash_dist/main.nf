#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MASH_DIST {
    tag "pairwise_distances"
    label 'process_high'
    container "quay.io/biocontainers/mash:2.3--he348c14_1"
    
    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.tsv"

    input:
    path sketches

    output:
    path "mash_distances.tsv", emit: distances
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create a combined sketch file
    echo "Combining \$(ls -1 *.msh | wc -l) sketch files..."
    mash paste combined *.msh

    # Calculate pairwise distances with parallelization
    # Using -p flag to utilize multiple threads
    echo "Computing pairwise distances with ${task.cpus} threads..."
    mash dist \
        -p ${task.cpus} \
        $args \
        combined.msh \
        combined.msh > mash_distances.tsv

    echo "Generated distance matrix with \$(wc -l < mash_distances.tsv) entries"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1 | sed 's/^/    /')
    END_VERSIONS
    """
}