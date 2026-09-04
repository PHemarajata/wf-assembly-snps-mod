#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * MASH_PASTE
 *
 * Combines the per-batch sketches from MASH_SKETCH_BATCH into a single .msh.
 *
 * Verified with mash 2.3 on 40 synthetic genomes:
 *   - `mash triangle` on one 40-genome sketch built in a single call and on the
 *     same 40 genomes sketched as 4 batches of 10 then pasted produced a
 *     BYTE-IDENTICAL Phylip matrix.
 *   - BUT the row order of the pasted file follows the order of the `mash paste`
 *     arguments: pasting the same 4 batches in reverse order changed the row
 *     order of the matrix (first label went from L0_G000 to L6_G000).
 * Hence the explicit LC_ALL=C sort of the batch files below: without it the
 * matrix row order would depend on Nextflow's staging order and vary between
 * otherwise identical runs.  Distances are unaffected either way; this is purely
 * about reproducible output.
 */

process MASH_PASTE {
    tag "combine_sketches"
    label 'process_low'
    conda "bioconda::mash=2.3"
    container "quay.io/biocontainers/mash:2.3--he348c14_1"

    publishDir "${params.outdir}/Clustering/Sketches", mode: params.publish_dir_mode, pattern: "combined.msh"

    input:
    path batch_sketches

    output:
    path "combined.msh", emit: sketch
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Deterministic paste order -> deterministic matrix row order.
    ls -1 *.msh | LC_ALL=C sort > sketch_list.txt

    if [ "\$(wc -l < sketch_list.txt)" -eq 1 ]; then
        cp "\$(cat sketch_list.txt)" combined.msh
    else
        mash paste combined \$(tr '\\n' ' ' < sketch_list.txt)
    fi

    n=\$(mash info -t combined.msh | tail -n +2 | grep -c . || true)
    echo "Combined sketch contains \$n genomes"
    if [ "\$n" -lt 2 ]; then
        echo "ERROR: combined sketch has \$n entries; clustering needs >= 2." >&2
        exit 1
    fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    mash: \$(mash --version 2>&1 | sed 's/^/    /')
END_VERSIONS
    """
}
