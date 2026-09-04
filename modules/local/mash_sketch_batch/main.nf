#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * MASH_SKETCH_BATCH
 *
 * Replaces the per-genome MASH_SKETCH scatter.  The legacy module ran one task
 * (and therefore one container start) per assembly: at 2000+ genomes that is
 * 2000+ container startups whose combined overhead dominates the actual
 * sketching work, and it never used the CPUs it was allocated.
 *
 * Here a batch of assemblies is sketched in a single `mash sketch -l` call with
 * -p ${task.cpus}.  Batching (rather than one single global task) is chosen so
 * that:
 *   - several batches run concurrently under Nextflow's own scheduler, so a
 *     modest workstation still saturates its cores without needing mash to
 *     scale to 20+ threads inside one process;
 *   - `-resume` stays useful: adding or changing one assembly invalidates one
 *     batch, not the entire sketching stage;
 *   - peak memory per task stays bounded and predictable.
 * The batch sketches are then combined with MASH_PASTE.
 *
 * IMPORTANT mash semantics that this module depends on (verified with mash 2.3):
 *   - `mash sketch -o out -l fof.txt` writes ONE out.msh containing one sketch
 *     per LINE of fof.txt (verified: 40 input files -> 40 sketches).
 *   - `mash info -t` on that file lists all of them; `mash dist`/`mash triangle`
 *     accept it directly.
 *   - `-m <int>` and `-r` IMPLY read-set mode, which POOLS every input into a
 *     single sketch.  With -l and 40 genomes, `-m 1` produced 1 sketch instead
 *     of 40.  The legacy per-genome module passed `-m ${params.mash_min_copies}`
 *     harmlessly because it only ever had one input file; here it would silently
 *     collapse the whole batch.  Therefore -m is only emitted when the user has
 *     explicitly raised it above 1, and doing so is refused for assemblies.
 *   - `-i` sketches individual SEQUENCES, so on multi-contig draft assemblies it
 *     yields one sketch per contig (verified: 3 files x 4 contigs -> 12
 *     sketches).  It must NOT be used here; row labels are normalized
 *     downstream instead.
 */

process MASH_SKETCH_BATCH {
    tag "batch_${batch_index}"
    label 'process_low'
    conda "bioconda::mash=2.3"
    container "quay.io/biocontainers/mash:2.3--he348c14_1"

    input:
    tuple val(batch_index), val(sample_ids), path(assemblies)

    output:
    path "batch_${batch_index}.msh", emit: sketch
    path "batch_${batch_index}.fof.txt", emit: fof
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def sketch_size = params.mash_sketch_size ?: 1000
    def kmer_size   = params.mash_kmer_size ?: 21
    def min_copies  = params.mash_min_copies ?: 1
    // See the -m note above: emitting -m implies -r, which pools the batch.
    def min_copies_arg = (min_copies as int) > 1 ? "-m ${min_copies}" : ''
    """
    if [ "${min_copies}" -gt 1 ]; then
        echo "ERROR: mash_min_copies=${min_copies} implies mash -r (read-set mode)," >&2
        echo "       which pools every file in the batch into a single sketch." >&2
        echo "       Assemblies must be sketched with the default -m 1." >&2
        exit 1
    fi

    # Deterministic file-of-filenames. Sorting matters: mash writes sketches in
    # input order, and MASH_PASTE preserves that order, which in turn fixes the
    # row order of the mash triangle matrix. Sorting here makes the matrix row
    # order reproducible across runs regardless of Nextflow staging order.
    for f in ${assemblies}; do
        echo "\$f"
    done | LC_ALL=C sort > batch_${batch_index}.fof.txt

    n_in=\$(wc -l < batch_${batch_index}.fof.txt)

    mash sketch \\
        -l batch_${batch_index}.fof.txt \\
        -o batch_${batch_index} \\
        -s ${sketch_size} \\
        -k ${kmer_size} \\
        -p ${task.cpus} \\
        ${min_copies_arg} \\
        ${args}

    # Fail loudly rather than silently shipping a pooled or truncated sketch.
    n_out=\$(mash info -t batch_${batch_index}.msh | tail -n +2 | grep -c . || true)
    if [ "\$n_out" -ne "\$n_in" ]; then
        echo "ERROR: sketched \$n_out entries from \$n_in input assemblies." >&2
        echo "       Expected one sketch per input file." >&2
        exit 1
    fi
    echo "Sketched \$n_out assemblies in batch ${batch_index}"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    mash: \$(mash --version 2>&1 | sed 's/^/    /')
END_VERSIONS
    """
}

/*
 * MASH_SKETCH_BATCH_PER_SAMPLE
 *
 * Optional compatibility path, enabled with --mash_per_sample_sketches true.
 *
 * A grep of workflows/ and subworkflows/ shows MASH_SKETCH.out.sketch is
 * consumed in exactly two places, and both do the same thing:
 *   workflows/recombination_aware_snps.nf:201
 *   subworkflows/local/clustering.nf:40
 *       MASH_SKETCH.out.sketch.map{ sample_id, sketch -> sketch }.collect()
 * i.e. the sample_id half of the tuple is discarded and the sketches are
 * collected into MASH_DIST.  Nothing needs a per-sample .msh file, so the
 * batched emit shape is safe.
 *
 * mash 2.3 offers no way to extract one sketch from a combined .msh, so this
 * path simply sketches each assembly separately *in one batched task* rather
 * than one task per genome.  It is still far cheaper than the legacy scatter
 * (no per-genome container start) while producing the same per-sample files
 * some users keep for screening new isolates against an existing panel.
 */
process MASH_SKETCH_BATCH_PER_SAMPLE {
    tag "batch_${batch_index}_per_sample"
    label 'process_low'
    conda "bioconda::mash=2.3"
    container "quay.io/biocontainers/mash:2.3--he348c14_1"

    publishDir "${params.outdir}/Clustering/Sketches", mode: params.publish_dir_mode, pattern: "*.msh"

    input:
    tuple val(batch_index), val(sample_ids), path(assemblies)

    output:
    tuple val(batch_index), path("*.msh"), emit: sketches
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sketch_size = params.mash_sketch_size ?: 1000
    def kmer_size   = params.mash_kmer_size ?: 21
    """
    for f in ${assemblies}; do
        base=\$(basename "\$f")
        base="\${base%.*}"
        mash sketch -o "\${base}" -s ${sketch_size} -k ${kmer_size} -p ${task.cpus} "\$f"
    done

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    mash: \$(mash --version 2>&1 | sed 's/^/    /')
END_VERSIONS
    """
}
