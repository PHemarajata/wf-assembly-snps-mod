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
          "!{core_alignment_fasta}"
    else
        # Use single tree builder
        run_gubbins.py \
          --starting-tree "!{meta.snp_package}.tree" \
          --prefix "!{meta.snp_package}-Gubbins" \
          --tree-builder "!{params.gubbins_tree_builder}" \
          --invariant-site-correction \
          --filter-percentage "!{params.gubbins_filter_percentage}" \
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
