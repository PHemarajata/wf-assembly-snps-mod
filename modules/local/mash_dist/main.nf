#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MASH_DIST {
    tag "pairwise_distances"
    label 'process_medium'
    container "quay.io/biocontainers/mash:2.3--he348c14_1"

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
    mash paste combined *.msh

    # Calculate pairwise distances
    mash dist \
        $args \
        combined.msh \
        combined.msh > mash_distances.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1 | sed 's/^/    /')
    END_VERSIONS
    """
}