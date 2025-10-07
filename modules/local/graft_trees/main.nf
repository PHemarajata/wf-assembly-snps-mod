#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GRAFT_TREES {
    tag "tree_grafting"
    label 'process_medium'
    container "quay.io/biocontainers/ete3:3.1.1--py35_0"
    
    // Publish placeholder graft outputs under Final_Results so downstream
    // consumers (and the external python script) find them at
    // <outdir>/Final_Results/...
    publishDir "${params.outdir}/Final_Results", mode: params.publish_dir_mode, pattern: "*.*"

    input:
    path tree_files
    path clusters_file
    path integrated_alignment

    // Grafting is deprecated/disabled in this workflow. Produce small placeholder
    // artifacts so downstream consumers that expect these files still succeed.
    output:
    // Match the filenames used by the external python script invocation
    path "global_grafted.treefile", emit: supertree
    path "grafting_report.txt", emit: report
    path "grafting_log.txt", emit: log
    path "cluster_tree_summary.tsv", emit: tree_summary
    path "supertree_visualization.pdf", optional: true, emit: visualization
    path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

        script:
    """
    # Grafting disabled: produce placeholder artifacts so downstream steps succeed.
    printf '();\n' > global_grafted.treefile
    printf 'GRAFTING DISABLED: The workflow is configured to use the backbone tree as the final tree.\n' > grafting_report.txt
    printf 'GRAFTING LOG\n' > grafting_log.txt
    printf 'cluster_id\ttree_file\tnum_leaves\ttree_length\tnum_internal_nodes\n' > cluster_tree_summary.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>/dev/null | sed 's/Python //' || echo "unknown")
        ete3: not_applicable
        pandas: not_applicable
        numpy: not_applicable
    END_VERSIONS
    """
}