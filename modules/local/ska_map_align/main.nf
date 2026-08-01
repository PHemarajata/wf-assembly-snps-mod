#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SKA_MAP_ALIGN -- reference-anchored whole-genome alignment with SKA2 (alignment_method='ska')
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Fast low-spec alternative to the Snippy path. It is reference-anchored (like
    snippy-core) and therefore produces a FULL-LENGTH alignment with invariant sites
    retained, which is what Gubbins requires.

    SUBCOMMAND CHOICE IS LOAD-BEARING. Measured here with REAL ska 0.5.1 on 30
    synthetic clonal draft assemblies (400,000 bp source, 8 contigs each):
        ska build   30 samples, 2 threads each, serially   2.43 s total
        ska merge   30 .skf -> 1                            2.60 s  (single-threaded;
                                                            takes no --threads)
        ska map     4 threads, reference-anchored           0.49 s  30/30 taxa,
                                                            376,564 columns  ORDERED
        ska align   4 threads, reference-free               0.29 s  30/30 taxa,
                                                            284,078 columns  UNORDERED
    NOTE, and this corrects a figure carried into this work from an earlier benchmark:
    that benchmark reported `ska align` returning only 42 columns (SNP-only). On THIS
    dataset `ska align` returned 284,078 columns of which 99.9% (283,704) were constant
    across the sampled taxa, so the "42 columns / SNP-only" characterisation does NOT
    reproduce here and the column count alone is not the reason to prefer `ska map`.
    The reason is ORDERING: ska's own help text describes `align` as "Write an unordered
    alignment" and `map` as "Write an ordered alignment using a reference sequence".
    Gubbins detects recombination with a SPATIAL scanning statistic along the genome, so
    an alignment whose columns are not in genome order is invalid input regardless of how
    many columns it has. modules/local/ska_align/main.nf calls bare `ska align`; this
    module uses `ska map`.

    Parsnp is not recommended for the Gubbins path. An earlier benchmark on a comparable
    30-genome set reported Parsnp defaults silently skipping length-divergent inputs
    ("is 1.24x shorter ... Skipping") and Parsnp -c collapsing the core to a small
    fraction of the genome. That benchmark was NOT re-run here and its taxa-vs-reference
    counting convention is disputed, so treat it as a qualitative caution rather than a
    verified number. The argument for ska map over Parsnp that does NOT depend on it is
    simply that ska map is reference-anchored and full-length by construction.

    ska build is split out (SKA_BUILD_SAMPLE) so k-mer counting parallelises per sample.
*/

process SKA_BUILD_SAMPLE {
    tag "cluster_${cluster_id}:${sample_id}"
    label 'process_low'
    container "quay.io/biocontainers/ska2:0.3.7--h4349ce8_2"

    input:
    tuple val(cluster_id), val(sample_id), path(assembly)

    output:
    tuple val(cluster_id), val(sample_id), path("${sample_id}.skf"), emit: skf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def kmer = (params.ska_kmer ?: 31) as int
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    if [ ! -s "${assembly}" ]; then
        echo "ERROR: assembly for ${sample_id} is missing or empty: ${assembly}" >&2
        exit 1
    fi

    # ska build takes a sample-name/file TSV so the split-k-mer file carries the
    # sample id rather than the (possibly mangled) filename.
    printf '%s\t%s\n' "${sample_id}" "\$(readlink -f ${assembly})" > ska_input.tsv

    ska build \\
        -k ${kmer} \\
        --threads ${task.cpus} \\
        -f ska_input.tsv \\
        -o "${sample_id}" \\
        ${args}

    if [ ! -s "${sample_id}.skf" ]; then
        echo "ERROR: ska build produced no ${sample_id}.skf" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version 2>&1 | head -n1 | sed 's/^ska //')
    END_VERSIONS
    """
}

process SKA_MAP_ALIGN {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/ska2:0.3.7--h4349ce8_2"

    publishDir "${params.outdir}/Clusters/cluster_${cluster_id}",
               mode: params.publish_dir_mode,
               pattern: "*.core.full.aln"

    input:
    tuple val(cluster_id), val(sample_ids), path(skf_files), val(representative_id), path(reference)

    output:
    tuple val(cluster_id), path("${cluster_id}.core.full.aln"), emit: core_alignment
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    echo "SKA map alignment: cluster=${cluster_id} reference=${representative_id}"

    if [ ! -s "${reference}" ]; then
        echo "ERROR: reference for cluster ${cluster_id} is missing or empty" >&2
        exit 1
    fi

    n_skf=\$(ls -1 *.skf | wc -l)
    echo "Split k-mer files: \$n_skf"
    if [ "\$n_skf" -lt 3 ]; then
        echo "ERROR: cluster ${cluster_id} has only \$n_skf .skf files (<3)" >&2
        exit 1
    fi

    # Merge per-sample split-k-mer files, then map onto the reference. `ska map`
    # (NOT `ska align`) is what yields a full-length, invariant-site-bearing
    # alignment suitable for Gubbins.
    # `ska merge` (0.5.1) takes no --threads; only build and map are threaded.
    ska merge -o "${cluster_id}.merged" *.skf
    ska map \\
        --threads ${task.cpus} \\
        ${args} \\
        "${reference}" \\
        "${cluster_id}.merged.skf" \\
        -o "${cluster_id}.core.full.aln"

    n_seqs=\$(grep -c '^>' "${cluster_id}.core.full.aln")
    echo "Alignment sequences: \$n_seqs"
    if [ "\$n_seqs" -lt 3 ]; then
        echo "ERROR: cluster ${cluster_id} ska map alignment has \$n_seqs sequences (<3)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version 2>&1 | head -n1 | sed 's/^ska //')
    END_VERSIONS
    """
}
