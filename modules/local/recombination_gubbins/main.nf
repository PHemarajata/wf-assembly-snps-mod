process RECOMBINATION_GUBBINS {

    tag { "${meta.snp_package}" }
    label "process_medium"
    // Pinned 2026-08-12 to the SAME BUILD as the validated production analysis:
    // conda `gubbins 3.4.3 py310h5140242_0 bioconda` == this biocontainer tag.
    // The bump from 3.3.5 and the --invariant-site-correction flag below MUST
    // move together: 3.4.2 made that correction OPTIONAL and defaulted it OFF,
    // so 3.4.3 without the flag silently drops a correction that 3.3.5 always
    // applied. Bumping the container alone would change results with no error.
    container "quay.io/biocontainers/gubbins:3.4.3--py310h5140242_0"
    
    publishDir "${params.outdir}/Recombination_Analysis", mode: params.publish_dir_mode, pattern: "*.{txt,tree}"

    input:
    tuple val(meta), path(core_alignment_fasta)
    tuple val(meta_alignment), path(alignment_files)

    output:
    tuple val(meta), path("*.{txt,tree}"), emit: positions_and_tree
    path(".command.{out,err}")
    path("versions.yml")                 , emit: versions

    shell:
    '''
    source $!{projectDir}/bin/bash_functions.sh

    msg "INFO: Performing recombination using Gubbins."

    # DETERMINISM, and a crash. Without --seed Gubbins derives RAxML's parsimony
    # seed from an unseeded randint(0, 10000). That is 0 about 1 time in 10,001,
    # RAxML rejects a non-positive seed, and Gubbins reports it only as "Unable
    # to fit model to data" -- ~16% chance per panel of losing a unit while still
    # exiting 0. See conf/params.config. Setting gubbins_seed = null restores the
    # old random draw and reintroduces both failure modes.
    # THREADS. This module previously passed no --threads at all, so Gubbins used
    # its own default of 1. The path was therefore single-threaded, and with a
    # seed it was accidentally deterministic, while being allocated
    # process_medium CPUs it never used. That is an implicit behaviour nobody
    # chose -- the same class of defect as a default pointing at the wrong run --
    # so it is now explicit and governed by the same parameter as the clustered
    # path. Default: use the allocated CPUs. gubbins_deterministic = true: 1.
    #
    # NOTE FOR THE RECORD: no reported result came from this path. The clustered
    # workflow produced the reported analysis. This is wired for consistency, so
    # that one parameter governs determinism everywhere rather than two paths
    # behaving differently for reasons nobody wrote down.
    # Either parameter turns it on: `deterministic` is the current name and covers
    # IQ-TREE too, `gubbins_deterministic` predates it and still works.
    GUBBINS_DETERMINISTIC="!{params.gubbins_deterministic}"
    DETERMINISTIC="!{params.deterministic}"
    if [ "${GUBBINS_DETERMINISTIC}" = "true" ] || [ "${DETERMINISTIC}" = "true" ]; then
      GUBBINS_THREADS=1
    else
      GUBBINS_THREADS="!{task.cpus}"
    fi

    GUBBINS_SEED="!{params.gubbins_seed}"
    SEED_ARG=""
    if [ -n "${GUBBINS_SEED}" ] && [ "${GUBBINS_SEED}" != "null" ]; then
      # Zero is refused rather than passed through: --seed 0 is the precise value
      # this parameter exists to make impossible, and setting it by hand would
      # reproduce the failure every run instead of 1 time in 10,001.
      if [ "${GUBBINS_SEED}" -le 0 ] 2>/dev/null; then
        msg "ERROR: gubbins_seed must be a positive integer, or null. RAxML rejects a non-positive parsimony seed. Got: ${GUBBINS_SEED}"
        exit 1
      fi
      SEED_ARG="--seed ${GUBBINS_SEED}"
    else
      msg "WARN: gubbins_seed is unset; Gubbins will draw a random RAxML seed and this run is not reproducible."
    fi

    # Build Gubbins command with hybrid tree builders if enabled
    if [[ "!{params.gubbins_use_hybrid}" == "true" ]]; then
        # Use hybrid approach with two tree builders
        run_gubbins.py \
          --starting-tree "!{meta.snp_package}.tree" \
          --prefix "!{meta.snp_package}-Gubbins" \
          --first-tree-builder "!{params.gubbins_first_tree_builder}" \
          --tree-builder "!{params.gubbins_tree_builder}" \
          --invariant-site-correction \
          --filter-percentage "!{params.gubbins_filter_percentage}" \
          --threads "${GUBBINS_THREADS}" \
          ${SEED_ARG} \
          "!{core_alignment_fasta}"
    else
        # Use single tree builder
        run_gubbins.py \
          --starting-tree "!{meta.snp_package}.tree" \
          --prefix "!{meta.snp_package}-Gubbins" \
          --tree-builder "!{params.gubbins_tree_builder}" \
          --invariant-site-correction \
          --filter-percentage "!{params.gubbins_filter_percentage}" \
          --threads "${GUBBINS_THREADS}" \
          ${SEED_ARG} \
          "!{core_alignment_fasta}"
    fi

    # Rename output files
    mv "!{meta.snp_package}-Gubbins.recombination_predictions.gff" \
      "!{meta.snp_package}-Gubbins.recombination_positions.txt"

    mv "!{meta.snp_package}-Gubbins.node_labelled.final_tree.tre" \
      "!{meta.snp_package}-Gubbins.labelled_tree.tree"

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        gubbins: $(run_gubbins.py --version | sed 's/^/    /')
    END_VERSIONS
    '''
}
