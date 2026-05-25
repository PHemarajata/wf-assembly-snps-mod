#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SUMMARIZE_CLUSTER_PHYLOGENY {
    tag "phylogeny_summary"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"

    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "cluster_phylogeny_summary.csv"
    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "cluster_membership.tsv"
    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "cluster_sizes.tsv"

    input:
    path clusters_tsv
    path diagnostics_files
    path tree_files
    path iqtree_log_files
    path filtered_alignment_files
    path recombination_gff_files

    output:
    path "cluster_phylogeny_summary.csv", emit: summary
    path "cluster_membership.tsv",        emit: membership
    path "cluster_sizes.tsv",             emit: sizes
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    # Stage per-cluster artifacts into named subdirectories so the helper
    # script can resolve files by cluster_id deterministically.
    mkdir -p diagnostics trees filtered recomb

    for f in ${diagnostics_files}; do
      [ -e "\$f" ] && cp -f "\$f" "diagnostics/\$(basename \$f)" || true
    done
    for f in ${tree_files}; do
      [ -e "\$f" ] && cp -f "\$f" "trees/\$(basename \$f)" || true
    done
    for f in ${iqtree_log_files}; do
      [ -e "\$f" ] && cp -f "\$f" "trees/\$(basename \$f)" || true
    done
    for f in ${filtered_alignment_files}; do
      [ -e "\$f" ] && cp -f "\$f" "filtered/\$(basename \$f)" || true
    done
    for f in ${recombination_gff_files}; do
      [ -e "\$f" ] && cp -f "\$f" "recomb/\$(basename \$f)" || true
    done

    # cluster_membership.tsv: same as clusters.tsv but published from this step
    cp -f "${clusters_tsv}" cluster_membership.tsv

    # cluster_sizes.tsv: one row per cluster_id with size
    python3 - <<'PY'
import csv
from collections import Counter
from pathlib import Path

sizes = Counter()
with open("${clusters_tsv}") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for row in reader:
        cid = row.get("cluster_id")
        if cid:
            sizes[cid] += 1

with open("cluster_sizes.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\\t")
    w.writerow(["cluster_id", "n_isolates"])
    for cid, n in sorted(sizes.items()):
        w.writerow([cid, n])
PY

    python3 ${projectDir}/bin/summarize_cluster_phylogeny.py \
        --clusters-tsv "${clusters_tsv}" \
        --diagnostics-dir diagnostics \
        --trees-dir trees \
        --filtered-dir filtered \
        --recomb-dir recomb \
        --output cluster_phylogeny_summary.csv

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //')
END_VERSIONS
    """
}
