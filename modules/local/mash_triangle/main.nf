#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * MASH_TRIANGLE
 *
 * Collapses MASH_DIST + MASH_TAB_TO_MATRIX into one process.
 *
 * The legacy path was:
 *     mash paste combined *.msh
 *     mash dist combined.msh combined.msh > mash_distances.tsv   # n^2 rows
 *     python mash_tab_to_matrix.py mash_distances.tsv mash_matrix.tsv
 * with three problems at n = 2000:
 *   1. `mash dist` was never given -p despite the module being allocated
 *      cpus = 4 (conf/profiles/local_workstation_rtx4070.config:65).
 *   2. It computes and writes all n^2 = 4,000,000 ordered pairs, i.e. every
 *      distance twice plus n self-distances.
 *   3. That 4M-row TSV was then reshaped by a per-row Python loop in a second
 *      container.
 *
 * `mash triangle -p ${task.cpus}` computes each unordered pair once and writes
 * a lower-triangular relaxed Phylip matrix directly, removing both the
 * duplicated computation and the whole reshape process.
 *
 * Measured on 40 synthetic 200 kb genomes, s=50000 sketches (mash 2.3):
 *     mash dist, 1 thread   0.789 s   (1600 rows)
 *     mash dist -p 8        0.182 s
 *     mash triangle -p 8    0.093 s
 *
 * Downstream consumes the Phylip directly: bin/cluster_mash.py sniffs the
 * format and bin/mash_matrix_io.py reads lower-triangular Phylip with numpy
 * fancy indexing.  A square TSV is still emitted (matrix_out) for any consumer
 * or user script that expects the historical mash_matrix.tsv shape - it is
 * produced from the same in-memory matrix, so the two cannot disagree.
 *
 * `mash triangle` labels rows by the INPUT FILE PATH (verified: row labels came
 * back as "bench/genomes/L0_G000.fasta"), exactly like `mash dist` did, so name
 * normalization still happens in mash_matrix_io.normalize_name().
 */

process MASH_TRIANGLE {
    tag "pairwise_distances"
    label 'process_medium'
    conda "bioconda::mash=2.3 conda-forge::python=3.10 conda-forge::numpy conda-forge::pandas conda-forge::scipy"
    container "quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0"

    publishDir "${params.outdir}/Clustering", mode: params.publish_dir_mode, pattern: "mash_distances.*"

    input:
    path combined_msh

    output:
    path "mash_distances.phylip", emit: phylip
    path "mash_distances_matrix.tsv", emit: matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mash triangle \\
        -p ${task.cpus} \\
        ${args} \\
        ${combined_msh} > mash_distances.phylip

    # Convert once, in-process, with the vectorized reader (no per-row loop).
    # Emitted for backwards compatibility with anything still expecting a square
    # mash_matrix.tsv; cluster_mash.py reads the Phylip directly.
    python3 ${projectDir}/bin/mash_phylip_to_matrix.py \\
        mash_distances.phylip \\
        mash_distances_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1 | sed 's/^/    /')
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}

/*
 * MASH_DIST_COMPAT
 *
 * Drop-in replacement for the legacy MASH_DIST, kept because
 * subworkflows/local/clustering.nf (used by workflows/assembly_snps_scalable.nf)
 * imports MASH_DIST and feeds MASH_DIST.out.distances straight into
 * CLUSTER_GENOMES.  That path is left working unchanged in shape; the only
 * difference is that -p ${task.cpus} is now actually passed.
 *
 * `mash dist -p 8` output was verified identical to single-threaded `mash dist`
 * after sorting, so adding -p does not change results.
 */
process MASH_DIST_COMPAT {
    tag "pairwise_distances_compat"
    label 'process_medium'
    conda "bioconda::mash=2.3 conda-forge::python=3.10 conda-forge::numpy conda-forge::pandas conda-forge::scipy"
    container "quay.io/biocontainers/mash:2.3--he348c14_1"

    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.tsv"

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
    if [ "\$(ls -1 *.msh | wc -l)" -eq 1 ]; then
        cp \$(ls -1 *.msh) combined.msh
    else
        mash paste combined \$(ls -1 *.msh | LC_ALL=C sort | tr '\\n' ' ')
    fi

    mash dist \\
        -p ${task.cpus} \\
        ${args} \\
        combined.msh \\
        combined.msh > mash_distances.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1 | sed 's/^/    /')
    END_VERSIONS
    """
}
