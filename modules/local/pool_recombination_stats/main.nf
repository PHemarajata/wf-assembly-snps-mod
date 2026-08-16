#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
================================================================================
    POOL_RECOMBINATION_STATS -- per-unit r/m, reference branches excluded
================================================================================

Gubbins reports SNPs inside and outside recombination per BRANCH. r/m is their
ratio pooled over a unit -- but the mapping reference is kept as a taxon, and it
sits outside the population by construction, so its branch is enormous and its
substitutions land in r/m's denominator. On the 2,070-genome run that was 52% of
all outside-recombination SNPs, and pooling naively gave a median r/m of 1.85
against 6.30 once those branches were dropped.

That correction used to live in a script run by hand after the pipeline
finished, which meant the pipeline's own r/m was wrong by a factor of ~3.4 and
stayed wrong unless someone remembered the extra step. It runs here instead.

See bin/pool_recombination_stats.py for which branches are dropped and why.
*/

process POOL_RECOMBINATION_STATS {
    tag "recombination_rm"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"

    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "recombination_rm.tsv"

    input:
    path clusters_tsv
    path per_branch_stats_files
    path final_tree_files

    output:
    path "recombination_rm.tsv", emit: rm
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_taxon = (params.reference_taxon_name ?: 'Reference').toString()
    """
    set -euo pipefail

    # Stage the two artefact classes into named directories so the helper can
    # pair <unit>.per_branch_statistics.csv with <unit>.node_labelled.final_tree.tre
    # by name. Both are already uniquely named per replicon-unit.
    mkdir -p stats trees

    for f in ${per_branch_stats_files}; do
      [ -e "\$f" ] && cp -f "\$f" "stats/\$(basename \$f)" || true
    done
    for f in ${final_tree_files}; do
      [ -e "\$f" ] && cp -f "\$f" "trees/\$(basename \$f)" || true
    done

    python3 ${projectDir}/bin/pool_recombination_stats.py \\
        --stats-dir stats \\
        --trees-dir trees \\
        --clusters-tsv "${clusters_tsv}" \\
        --reference-taxon "${ref_taxon}" \\
        --output recombination_rm.tsv

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    python: \$(python3 --version 2>&1 | sed 's/Python //')
END_VERSIONS
    """
}
