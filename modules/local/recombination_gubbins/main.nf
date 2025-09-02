process RECOMBINATION_GUBBINS {

    tag { "${meta.snp_package}" }
    label "process_medium"
    container "quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0"

    input:
    tuple val(meta), path(core_alignment_fasta)
    tuple val(meta_alignment), path(alignment_files)

    output:
    tuple val(meta), path("*.{txt,tree}"), emit: positions_and_tree
    path(".command.{out,err}")
    path("versions.yml")                 , emit: versions

    shell:
    '''
    source bash_functions.sh

    msg "INFO: Performing recombination using Gubbins."

    # Build Gubbins command with hybrid tree builders if enabled
    if [[ "!{params.gubbins_use_hybrid}" == "true" ]]; then
        # Use hybrid approach with two tree builders
        run_gubbins.py \
          --starting-tree "!{meta.snp_package}.tree" \
          --prefix "!{meta.snp_package}-Gubbins" \
          --first-tree-builder "!{params.gubbins_first_tree_builder}" \
          --tree-builder "!{params.gubbins_tree_builder}" \
          "!{core_alignment_fasta}"
    else
        # Use single tree builder
        run_gubbins.py \
          --starting-tree "!{meta.snp_package}.tree" \
          --prefix "!{meta.snp_package}-Gubbins" \
          --tree-builder "!{params.gubbins_tree_builder}" \
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
