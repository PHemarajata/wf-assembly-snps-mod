#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * CLUSTER_GENOMES
 *
 * Changes relative to the original module:
 *
 * 1. The runtime `pip install pandas numpy scipy` is GONE.  It ran inside every
 *    container invocation, needed outbound network from the compute node, was
 *    completely unpinned (so the numeric stack could differ between two runs of
 *    the same pipeline), and its versions.yml then reported whatever pip
 *    happened to fetch.  The container now ships the stack.
 *
 * 2. Container pinned to quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0.
 *    HOW THIS WAS VERIFIED: the image manifest was pulled from the quay.io
 *    registry API and the conda-meta/ records inside its largest layer were
 *    listed.  The image contains
 *        python 3.8.20, pandas 2.0.3, numpy 1.22.4, scipy 1.7.0,
 *        scikit-learn 1.3.2, and mash 2.3
 *    Note biocontainers/python:3.9--1 (the image the original module used) ships
 *    NONE of these - hence the pip install.  biocontainers/pandas:2.2.1 has
 *    pandas+numpy but NO scipy, which cluster_mash.py requires for
 *    connected_components and the average-linkage split.  pyseer was chosen
 *    because it is the smallest verified image (273 MB) carrying the full
 *    pandas+numpy+scipy set AND the mash binary, so the whole Mash/clustering
 *    front end can share one image and one Docker pull.
 *
 * 3. --split-method is passed through from params.cluster_split_method so that
 *    the phylogeny-aware splitting of oversized clusters can be turned off.
 *
 * 4. --emit-submatrices writes one small per-cluster distance TSV, so
 *    SELECT_CLUSTER_REPRESENTATIVE no longer re-parses the full n x n matrix in
 *    every per-cluster task.
 *
 * The input is whatever MASH_TRIANGLE or MASH_DIST_COMPAT emitted; the format
 * (lower-triangular Phylip / square TSV / raw mash dist edges) is sniffed by
 * bin/mash_matrix_io.py, so this module works on both front ends.
 */

process CLUSTER_GENOMES {
    tag "clustering"
    label 'process_low'
    conda "conda-forge::python=3.10 conda-forge::numpy conda-forge::pandas conda-forge::scipy"
    container "quay.io/biocontainers/pyseer:1.4.2--pyhdfd78af_0"

    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.{tsv,txt}"

    input:
    path mash_distances

    output:
    path "clusters.tsv", emit: clusters
    path "cluster_summary.txt", emit: summary
    path "cluster_matrices/*.matrix.tsv", emit: submatrices, optional: true
    path "threshold_sweep.tsv",  emit: threshold_sweep,  optional: true
    path "chosen_threshold.txt", emit: chosen_threshold, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def max_cluster_size = params.max_cluster_size ?: 100
    def split_method     = params.cluster_split_method ?: 'similarity'
    def consol_tol       = params.cluster_consolidate_tolerance ?: 2.0
    // `--mash_threshold auto` sweeps candidate thresholds and picks one from the
    // collection's own distance distribution. A fixed threshold does not transfer
    // between collections: 0.003, derived on 112 genomes, silently drops 203 of
    // 2795 on a wider set because they land in components of <3 taxa.
    def auto             = (params.mash_threshold ?: '').toString().trim().toLowerCase() == 'auto'
    def threshold        = auto ? '' : (params.mash_threshold ?: 0.03)
    def coherence        = params.auto_threshold_coherence ?: 0.90
    def grid             = params.auto_threshold_grid ? "--auto-grid ${params.auto_threshold_grid}" : ''
    def mode             = auto ? "--auto-threshold --auto-coherence ${coherence} ${grid} --auto-report threshold_sweep.tsv --threshold-out chosen_threshold.txt" : "--threshold ${threshold}"
    """
    mkdir -p cluster_matrices

    python3 ${projectDir}/bin/cluster_mash.py \\
        ${mash_distances} \\
        clusters.tsv \\
        ${mode} \\
        --max-cluster-size ${max_cluster_size} \\
        --split-method ${split_method} \\
        --consolidate-tolerance ${consol_tol} \\
        --emit-submatrices cluster_matrices \\
        ${args} > cluster_summary.txt

    cat cluster_summary.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //')
    pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
    scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    cluster_split_method: ${split_method}
END_VERSIONS
    """
}
