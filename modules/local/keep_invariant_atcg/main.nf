#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KEEP_INVARIANT_ATCG -- column filter on the whole-genome alignment fed to Gubbins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Gubbins requires a whole-genome alignment ("An alignment of polymorphic sites
    cannot be used as the input to Gubbins") because recombination detection is a
    SPATIAL scanning statistic -- the invariant sites between SNPs carry signal.
    This process therefore prunes columns only for missingness, never for variability.

    WHAT CHANGED (see ALIGNMENT_TREES_NOTES.md):
      1. The per-column Biopython loop was replaced by a vectorised NumPy byte-matrix
         implementation in bin/filter_alignment_columns.py. Measured on 50 taxa x 200 kb:
         6.38 s -> 0.088 s of compute (72x), 0.166 s end-to-end including interpreter
         start (38x). Output is byte-identical at --max_column_missingness 0.0.
      2. The all-or-nothing A/T/C/G rule became a per-column missingness threshold,
         params.max_column_missingness (default 0.10). ** THIS ALTERS SCIENTIFIC OUTPUT. **
         Setting max_column_missingness = 0.0 restores the legacy rule exactly (verified
         byte-identical at 24 and 50 taxa).
      3. The error handler that wrote ">${cluster_id}_dummy / ATCG" on any exception is
         gone. A 4 bp placeholder alignment silently entering Gubbins is data corruption,
         not recovery, so the process now fails.
      4. storeDir was removed. storeDir keys on filename only, so a cached placeholder
         (or an alignment produced under a different missingness threshold) would be
         reused forever and would ignore any change to params. Normal Nextflow work-dir
         caching with `cache true` keys on inputs AND the script text, so changing the
         threshold correctly invalidates. `-resume` still works.
*/

process KEEP_INVARIANT_ATCG {
    tag "cluster_${cluster_id}"
    label 'process_low'
    container "quay.io/biocontainers/biopython@sha256:10d755c731c82a22d91fc346f338ba47d5fd4f3b357828f5bbc903c9be865614"

    publishDir "${params.outdir}/Clusters/cluster_${cluster_id}",
               mode: params.publish_dir_mode,
               pattern: "*.column_filter.stats.tsv"

    input:
    tuple val(cluster_id), path(alignment)

    output:
    tuple val(cluster_id), path("${cluster_id}.core.full.aln"),               emit: core_alignment
    tuple val(cluster_id), path("${cluster_id}.column_filter.stats.tsv"),     emit: stats
    path "versions.yml",                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Per-column missingness threshold. 0.0 == legacy all-or-nothing A/T/C/G rule.
    def max_missing = (params.max_column_missingness == null ? 0.10 : params.max_column_missingness) as double
    def min_cols    = (params.min_kept_columns ?: 1) as int
    """
    set -euo pipefail

    echo "Column-filtering alignment for cluster ${cluster_id}"
    echo "Input alignment: ${alignment}"
    echo "max_column_missingness=${max_missing} (0.0 reproduces the legacy all-ATCG rule)"

    # Fail loudly on a missing/empty input rather than manufacturing a placeholder.
    if [ ! -s "${alignment}" ]; then
        echo "ERROR: input alignment for cluster ${cluster_id} is missing or empty: ${alignment}" >&2
        exit 1
    fi

    filter_alignment_columns.py \\
        --input  "${alignment}" \\
        --output "${cluster_id}.core.full.aln" \\
        --stats  "${cluster_id}.column_filter.stats.tsv" \\
        --max-missingness ${max_missing} \\
        --min-kept-columns ${min_cols} \\
        --label "${cluster_id}"

    # Post-conditions: a real alignment, not a stub. Any violation is a hard error.
    n_seqs=\$(grep -c '^>' "${cluster_id}.core.full.aln")
    n_cols=\$(awk 'NR>1 && /^>/{exit} NR>1{gsub(/[ \\t\\r\\n]/,"");c+=length(\$0)} END{print c+0}' "${cluster_id}.core.full.aln")
    echo "Output: \$n_seqs sequences x \$n_cols columns"
    if [ "\$n_seqs" -lt 3 ]; then
        echo "ERROR: cluster ${cluster_id} filtered alignment has \$n_seqs sequences (<3)" >&2
        exit 1
    fi
    if [ "\$n_cols" -lt ${min_cols} ]; then
        echo "ERROR: cluster ${cluster_id} filtered alignment has \$n_cols columns (< min_kept_columns=${min_cols})" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
