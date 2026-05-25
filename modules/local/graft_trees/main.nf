#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GRAFT_TREES {
    tag "graft_trees"
    label 'process_medium'
    container "quay.io/biocontainers/biopython:1.81"

    publishDir "${params.outdir}/Final_Results", mode: params.publish_dir_mode,
               pattern: "global_grafted.treefile"
    publishDir "${params.outdir}/Final_Results", mode: params.publish_dir_mode,
               pattern: "grafting_{report,log}.txt"

    input:
    path backbone_tree
    path cluster_trees
    path reps_tsv

    output:
    path "global_grafted.treefile", emit: grafted_tree
    path "grafting_report.txt",     emit: report
    path "grafting_log.txt",        emit: log
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def parent_edge_mode = (params.graft_parent_edge_mode ?: 'keep').toString()
    def rename_flag      = params.graft_rename_conflicts ? '--rename-conflicts' : ''
    """
    set -euo pipefail

    # Stage cluster trees into Clusters/ so the script's glob can find them
    # deterministically. Skip empty trees (failed clusters) so they don't
    # pollute the grafting report.
    mkdir -p Clusters
    for f in ${cluster_trees}; do
      [ -e "\$f" ] || continue
      [ -s "\$f" ] || continue
      cp -f "\$f" "Clusters/\$(basename \$f)"
    done

    # Optional representatives TSV. The graft_trees.py script can infer
    # representatives from cluster tip labels if no TSV is provided, but a
    # proper mapping makes the grafting deterministic.
    reps_arg=""
    if [ -s "${reps_tsv}" ] && [ "\$(basename "${reps_tsv}")" != "NO_FILE" ]; then
      reps_arg="--reps ${reps_tsv}"
    fi

    python3 ${projectDir}/graft_trees.py \\
        --backbone "${backbone_tree}" \\
        --clusters 'Clusters/*.final.treefile' \\
        --out-tree global_grafted.treefile \\
        --report   grafting_report.txt \\
        --log      grafting_log.txt \\
        --parent-edge-mode ${parent_edge_mode} \\
        ${rename_flag} \\
        \${reps_arg}

    # Defensive: the script always writes these, but guarantee channel outputs
    # exist even on partial failure so Nextflow can proceed.
    [ -s global_grafted.treefile ] || cp "${backbone_tree}" global_grafted.treefile
    [ -f grafting_report.txt ]     || : > grafting_report.txt
    [ -f grafting_log.txt ]        || : > grafting_log.txt

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python:    \$(python3 --version | sed 's/Python //')
    biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
END_VERSIONS
    """
}
