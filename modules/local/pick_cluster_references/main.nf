#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * PICK_CLUSTER_REFERENCES
 *
 * Choose a per-cluster MAPPING reference: a complete member where one exists,
 * otherwise the nearest complete genome borrowed from the collection.
 *
 * WHY THIS IS NOT SELECT_CLUSTER_REPRESENTATIVE. The representative used for the
 * BACKBONE and the reference used for MAPPING are different objects with
 * different requirements. The representative is a medoid; the reference must
 * additionally be COMPLETE, because Gubbins cannot use a multi-contig reference
 * and --split_replicons needs one FASTA per replicon. Conflating the two is what
 * previously selected a 135-contig draft as a cluster reference.
 *
 * That distinction is decisive on real data: 2,541 of the 2,802 B. pseudomallei
 * genomes (90.7%) are drafts, so a plain medoid is usually unusable as a
 * reference and --split_replicons would fail on it by design.
 *
 * SELECTION: completeness is a GATE (pass/fail on --max-contigs), centrality is
 * the RANKING (mean then max Mash distance to cluster members). N50 and total
 * length are deliberately not ranking criteria -- once the gate is passed a
 * longer reference only adds positions other members cannot fill.
 *
 * VALIDATED against the manual analysis's own decisions on its 37 analysed
 * units: 36/36 agreement on internal-vs-borrowed, and 29/36 identical reference
 * when given the same 14-reference candidate shortlist.
 *
 * Reuses MASH_TRIANGLE's matrix, so it costs no extra alignment work.
 */

process PICK_CLUSTER_REFERENCES {
    tag "reference_selection"
    label 'process_low'
    conda "conda-forge::python=3.10 conda-forge::numpy conda-forge::pandas"
    container "quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0"

    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.tsv"

    input:
    path clusters_tsv
    path mash_matrix
    path genomes, stageAs: 'genomes/*'
    path blocklist

    output:
    path "cluster_references.tsv",   emit: references
    path "reference_selection.tsv",  emit: report
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_contigs = params.max_replicons ?: 4
    def min_size    = params.min_cluster_size ?: 2
    def pool        = params.reference_pool ? "--reference-pool ${params.reference_pool}" : ''
    // The blocklist is staged as a file so its CONTENTS are part of the task
    // hash: editing it after a bad reference is identified must invalidate the
    // selection, not silently reuse a cached pick made without it.
    def block       = blocklist.name != 'NO_BLOCKLIST' ? "--reference-blocklist ${blocklist}" : ''
    """
    set -euo pipefail

    python3 ${projectDir}/bin/pick_cluster_references.py \\
        --membership ${clusters_tsv} \\
        --fasta-dir genomes \\
        --phylip ${mash_matrix} \\
        --max-contigs ${max_contigs} \\
        --min-size ${min_size} \\
        ${pool} \\
        ${block} \\
        --out cluster_references.tsv \\
        --report reference_selection.tsv

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //')
    max_contigs: ${max_contigs}
    reference_pool: ${params.reference_pool ?: 'all_complete_genomes'}
END_VERSIONS
    """

    stub:
    // The reference path must point at a file that actually exists: the
    // workflow resolves it with checkIfExists, which is exactly the guard that
    // should stay on. Emit a real (empty) FASTA in the task dir rather than a
    // placeholder path.
    // Emit one reference per cluster actually present in the input, and use a
    // REAL staged genome -- the one with fewest contigs -- rather than a
    // synthetic FASTA. SNIPPY_CORE_GATHER, CLUSTER_GENOMES and GUBBINS_CLUSTER
    // define no stub blocks, so `-stub-run` executes their real scripts; a
    // 10 bp placeholder reference makes snippy-core fail and the stub then
    // tests nothing past this point. A stub that covers only 'strain_1' would
    // likewise leave every other cluster unjoined and silently dropped.
    """
    best=""; bestn=999999
    for f in genomes/*; do
      n=\$(grep -c '^>' "\$f" || echo 999999)
      if [ "\$n" -lt "\$bestn" ]; then bestn=\$n; best="\$f"; fi
    done
    cp "\$best" stub_reference.fasta
    echo "stub reference: \$best (\$bestn contigs)"

    printf 'cluster_id\\treference_path\\n'                        > cluster_references.tsv
    printf 'cluster_id\\tn\\treference\\tsource\\n'                > reference_selection.tsv
    for c in \$(awk 'NR>1{print \$1}' ${clusters_tsv} | sort -u); do
      printf '%s\\t%s/stub_reference.fasta\\n' "\$c" "\$PWD"     >> cluster_references.tsv
      printf '%s\\t%s\\tstub_reference\\tinternal\\n' "\$c" "\$bestn" >> reference_selection.tsv
    done
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: 3.10.0
END_VERSIONS
    """
}
