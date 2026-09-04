#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * POPPUNK_CLUSTER
 *
 * Partition the collection with PopPUNK instead of Mash single-linkage, and emit
 * the SAME `cluster_id<TAB>sample_id` contract that CLUSTER_GENOMES and curated
 * mode already produce. Everything downstream is therefore unchanged -- this
 * swaps the clustering front end only.
 *
 * WHY THIS EXISTS. Mash single-linkage cannot partition this collection at all.
 * Measured on the real 2,802-genome B. pseudomallei matrix: one connected
 * component containing every genome at threshold >= 0.007, shattering into
 * 786/1,094 singletons below it, with no value in between that works. What the
 * Mash path then emits is a size-capped chop of one component (Gini 0.059,
 * max/min 2.78) -- an imposed partition, not a found one. PopPUNK's refined fit
 * on the same 2,802 genomes gives 271 strains with a sane size distribution.
 *
 * PARAMETERS are Seng et al. 2024's published in-organism fit (PMID 38972886),
 * the only configuration demonstrated to feed Gubbins successfully in
 * B. pseudomallei: --min-k 15 --max-k 31 --k-step 2 --K 4 --max-a-dist 0.53.
 * Their fit reported network score 0.8961 against a documented bar of >= 0.8;
 * ours is echoed into the log so it can be checked the same way.
 *
 * BOTH FIT STEPS ARE REQUIRED. Measured on this collection:
 *     --fit-model bgmm    -> 208 clusters, largest 1,723, 121 singletons
 *     --fit-model refine  -> 271 clusters, largest   913, 158 singletons
 * Stopping after bgmm leaves a 1,723-genome "cluster" that is not an analysis
 * unit by any definition. The refine pass is not optional.
 *
 * CONTAINER pinned to PopPUNK 2.7.6 -- the exact version whose refined fit
 * produced the validated 271-strain partition in poppunk_bp/. Do not float this
 * to 2.7.7/2.7.8 without re-deriving the partition and checking the strain count
 * and size distribution, for the same reason the Gubbins builder is pinned.
 *
 * TWO WORKAROUNDS carried over from run_poppunk_bp.sh, both load-bearing:
 *
 *  1. PopPUNK builds some output paths as os.path.join(prefix, prefix + ext), so
 *     a multi-component --output such as "out/db" becomes "out/db/out/db.png"
 *     and crashes on a missing directory. Run from inside the work dir and use
 *     FLAT output names.
 *
 *  2. The rfile (<name>\\tab<path>, one line per genome) must be complete. A
 *     truncated single-line rfile against a clusters file with thousands of rows
 *     is what stalled PopPIPE-bp, silently. The guard below refuses to proceed
 *     on a suspiciously short rfile rather than producing a one-genome partition.
 */

process POPPUNK_CLUSTER {
    tag "poppunk"
    label 'process_high'
    conda "bioconda::poppunk=2.7.6"
    container "quay.io/biocontainers/poppunk:2.7.6--py310h4d0eb5b_0"

    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.{tsv,txt}"
    publishDir "${params.outdir}/PopPUNK",   mode: params.publish_dir_mode, pattern: "*_clusters.csv"

    input:
    path genomes, stageAs: 'genomes/*'

    output:
    path "clusters.tsv",           emit: clusters
    path "cluster_summary.txt",    emit: summary
    path "refined_clusters.csv",   emit: poppunk_clusters
    path "poppunk_excluded.tsv",   emit: excluded, optional: true
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_k     = params.poppunk_min_k     ?: 15
    def max_k     = params.poppunk_max_k     ?: 31
    def k_step    = params.poppunk_k_step    ?: 2
    def kval      = params.poppunk_K         ?: 4
    def max_a     = params.poppunk_max_a_dist ?: 0.53
    def min_size  = params.min_cluster_size  ?: 7
    def do_refine = (params.poppunk_refine == null) ? true
                    : params.poppunk_refine.toString().toLowerCase() in ['true','1','yes']
    """
    set -euo pipefail

    # --- rfile: <name>\\t<path>, one line per genome -------------------------
    : > rfile.txt
    for f in genomes/*; do
      n=\$(basename "\$f")
      n="\${n%.*}"
      printf '%s\\t%s\\n' "\$n" "\$(readlink -f "\$f")" >> rfile.txt
    done
    N=\$(wc -l < rfile.txt)
    echo "rfile lines: \$N"

    # Guard: the PopPIPE-bp failure mode was a one-line rfile that ran to
    # completion and produced a meaningless partition. Fail loudly instead.
    if [ "\$N" -lt 10 ]; then
      echo "ERROR: rfile has only \$N lines -- refusing to cluster. This is the" >&2
      echo "       truncated-rfile failure mode, not a small dataset." >&2
      exit 2
    fi

    # --- 1. sketch -----------------------------------------------------------
    echo "=== 1. create database (sketching) ==="
    poppunk --create-db \\
        --r-files rfile.txt \\
        --output db \\
        --min-k ${min_k} --max-k ${max_k} --k-step ${k_step} \\
        --threads ${task.cpus}

    # --- 2. bgmm fit ---------------------------------------------------------
    echo "=== 2. fit model (bgmm, Seng et al. 2024 parameters) ==="
    poppunk --fit-model bgmm \\
        --ref-db db \\
        --output fit \\
        --K ${kval} --max-a-dist ${max_a} \\
        --threads ${task.cpus}

    FINAL_DIR=fit
    """ + (do_refine ? """
    # --- 3. refine -----------------------------------------------------------
    # Not optional: bgmm alone left a 1,723-genome cluster on this collection.
    echo "=== 3. refine model ==="
    poppunk --fit-model refine \\
        --ref-db db \\
        --model-dir fit \\
        --output refined \\
        --threads ${task.cpus}
    FINAL_DIR=refined
    """ : """
    echo "=== 3. refine SKIPPED (poppunk_refine = false) ==="
    echo "WARNING: without the refine pass the largest cluster is typically far" >&2
    echo "         too large to be an analysis unit." >&2
    """) + """

    # Normalise the clusters file name regardless of which stage produced it.
    CL=\$(ls \${FINAL_DIR}/*_clusters.csv 2>/dev/null | grep -v unword | head -1)
    if [ -z "\$CL" ]; then
      echo "ERROR: PopPUNK produced no clusters CSV in \${FINAL_DIR}" >&2
      exit 3
    fi
    cp "\$CL" refined_clusters.csv

    # Network score, against the documented >= 0.8 bar.
    grep -h -i "network score\\|Network summary" -A4 .command.log 2>/dev/null | head -12 || true

    # --- convert to the cluster_id<TAB>sample_id contract --------------------
    python3 ${projectDir}/bin/poppunk_clusters_to_tsv.py \\
        --clusters refined_clusters.csv \\
        --rfile rfile.txt \\
        --min-cluster-size ${min_size} \\
        --out clusters.tsv \\
        --excluded poppunk_excluded.tsv \\
        > cluster_summary.txt

    cat cluster_summary.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    poppunk: \$(poppunk --version 2>&1 | sed 's/^poppunk //')
    poppunk_fit: bgmm${do_refine ? '+refine' : ''}
    min_cluster_size: ${min_size}
END_VERSIONS
    """

    stub:
    // Emit the REAL staged sample names, not placeholders. A stub that invents
    // names produces a clusters.tsv that joins against nothing, so every
    // downstream process is silently skipped and the stub run "passes" having
    // exercised only this module. Splitting the genomes across two clusters
    // also keeps the per-cluster fan-out meaningful.
    """
    printf 'cluster_id\\tsample_id\\n' > clusters.tsv
    printf 'Taxon,Cluster\\n'          > refined_clusters.csv
    i=0
    for f in genomes/*; do
      n=\$(basename "\$f"); n="\${n%.*}"
      c=\$(( i % 2 + 1 ))
      printf 'strain_%s\\t%s\\n' "\$c" "\$n" >> clusters.tsv
      printf '%s,%s\\n' "\$n" "\$c"          >> refined_clusters.csv
      i=\$(( i + 1 ))
    done
    echo "stub: \$(( i )) genomes into 2 clusters" > cluster_summary.txt
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    poppunk: 2.7.6
END_VERSIONS
    """
}
