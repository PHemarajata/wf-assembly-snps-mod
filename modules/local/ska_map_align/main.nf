#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SKA_MAP_ALIGN -- reference-anchored whole-genome alignment with SKA2 (alignment_method='ska')
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Fast low-spec alternative to the Snippy path. It is reference-anchored (like
    snippy-core) and therefore produces a FULL-LENGTH alignment with invariant sites
    retained, which is what Gubbins requires.

    *** DO NOT USE THIS PATH FOR RECOMBINATION-AWARE ANALYSIS. ***

    It is 2.7x faster end to end (9m11s vs 24m53s on the 112-genome set, 1.0 vs 3.0
    CPU-hours, alignment CPU 652 s vs 7333 s) and it is WRONG for this workflow,
    because split k-mers cannot see clustered SNPs.

    A split k-mer matches only when both flanking half-k-mers match exactly, so a
    second SNP inside the flank destroys the match and the SNP is never called. The
    loss is therefore a function of SNP SPACING. Measured against the Snippy
    alignment of the same cluster, same samples, same reference -- ratio of SKA
    variable sites to Snippy variable sites, binned by distance to the next SNP:

        gap to next SNP   cluster_12 (21 taxa)   cluster_9 (12 taxa)
             1-10 bp            0.11x                  0.05x
            11-20 bp            0.41x                  0.34x
            21-31 bp            0.75x                  0.72x
            32-50 bp            0.78x                  0.75x
           51-100 bp            0.86x                  0.83x
          101-500 bp            0.96x                  0.93x
         501-5000 bp            1.09x                  1.06x

    Monotonic in spacing and at parity beyond ~100 bp, which is the signature of the
    k-mer flank and not of coverage, missing data or divergence. Both alignments are
    essentially complete (missing 0.16% Snippy vs 0.13% SKA), so this is not a
    data-loss artefact. SNPs within 31 bp of a neighbour: 27.1% of Snippy intervals
    vs 12.2% of SKA's on cluster_12; 24.2% vs 9.1% on cluster_9.

    Dense SNP clusters are exactly what Gubbins keys on -- it detects recombination
    as regions of elevated SNP density. Deleting ~90% of the tightest clusters
    deletes the signal, and the error is amplified on the way through: 23% fewer
    variable sites produced 54% FEWER RECOMBINATION BLOCKS on the 112-genome set
    (2,388 vs 5,148). Per-cluster topologies then disagreed in 3 of the 4 clusters
    large enough to have a non-trivial topology.

    The failure mode is systematic FALSE NEGATIVES for recombination in a highly
    recombinogenic organism, i.e. the one error this workflow exists to avoid. The
    less-recombination result also looks perfectly plausible in the output, which is
    what makes it dangerous.

    SKA may still be reasonable for a recombination-FREE use (quick clustering or
    distance estimates). It is not an alternative to Snippy here.

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
    conda "bioconda::ska2=0.3.7"
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
    conda "bioconda::ska2=0.3.7"
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
        -o "${cluster_id}.ska_raw.aln"

    # ska map emits IUPAC ambiguity codes where a split k-mer matched more than one
    # allele. Gubbins rejects the whole alignment on the first one it sees:
    #
    #   Error with the input FASTA file: <sample> contains disallowed characters,
    #   only ACGTNacgtn- are permitted
    #
    # and exits 1, so the SKA path could never reach recombination detection on real
    # data. It is rare but not negligible at scale -- measured on the 21-taxon
    # cluster_12 of the 112-genome set, 841 of 145,542,579 bases (0.00058%): R 315,
    # Y 290, S 107, M 59, K 40, W 26, V 4. A 3-taxon smoke cluster has too few
    # samples to hit one, which is why this never surfaced before.
    #
    # Ambiguity is folded to N rather than resolved: an ambiguous base is exactly
    # "unknown" for recombination detection, which is what N means to Gubbins.
    # Guessing a majority allele would invent data at the sites least able to
    # support it.
    awk '/^>/ { print; next } { gsub(/[RYSWKMBDHVryswkmbdhv]/, "N"); print }' \\
        "${cluster_id}.ska_raw.aln" > "${cluster_id}.core.full.aln"

    n_amb=\$(grep -v '^>' "${cluster_id}.ska_raw.aln" \\
        | tr -cd 'RYSWKMBDHVryswkmbdhv' | wc -c)
    echo "Folded \$n_amb IUPAC ambiguity code(s) to N for Gubbins compatibility"
    rm -f "${cluster_id}.ska_raw.aln"

    # Nothing outside the Gubbins-permitted set may survive, or Gubbins fails late
    # and the diagnosis is a container error rather than this message.
    n_bad=\$(grep -v '^>' "${cluster_id}.core.full.aln" | tr -d 'ACGTNacgtn\\n-' | wc -c)
    if [ "\$n_bad" -ne 0 ]; then
        echo "ERROR: \$n_bad character(s) outside ACGTNacgtn- remain in \\
${cluster_id}.core.full.aln; Gubbins will reject it" >&2
        grep -v '^>' "${cluster_id}.core.full.aln" | tr -d 'ACGTNacgtn\\n-' \\
            | fold -w1 | sort -u | tr '\\n' ' ' >&2
        exit 1
    fi

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
