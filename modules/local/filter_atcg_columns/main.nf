#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * FILTER_ATCG_COLUMNS
 * Keep only columns where every sequence has A/T/C/G (case-insensitive).
 * Emits (cluster_id, filtered.fa) so we can join on cluster_id later.
 */

process FILTER_ATCG_COLUMNS {
    tag "${cluster_id}"
    label 'process_low'
    container "quay.io/biocontainers/biopython@sha256:10d755c731c82a22d91fc346f338ba47d5fd4f3b357828f5bbc903c9be865614"

    input:
        tuple val(cluster_id), path(alignment)

    output:
        tuple val(cluster_id), path("${cluster_id}.filtered.fa"), emit: filtered_alignment
        path "versions.yml",                                        emit: versions

    script:
    """
    set -euo pipefail

    python3 - <<PYCODE
from Bio import AlignIO

aln = AlignIO.read("${alignment}", "fasta")
L = aln.get_alignment_length()

keep = []
for i in range(L):
    col = aln[:, i]
    if all(base.upper() in "ATCG" for base in col):
        keep.append(i)

with open("${cluster_id}.filtered.fa", "w") as out:
    for rec in aln:
        seq = "".join(str(rec.seq[i]) for i in keep)
        out.write(f">{rec.id}\\n{seq}\\n")
PYCODE

    echo "FILTER_ATCG_COLUMNS:" > versions.yml
    python3 -c 'import Bio; print("  biopython: " + Bio.__version__)' >> versions.yml
    """
}