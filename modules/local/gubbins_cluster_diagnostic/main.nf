#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Gubbins clustering with strict diagnostics and output checks.
 * Emits .diagnostics.log for troubleshooting.
 */

process GUBBINS_CLUSTER_DIAGNOSTIC {
    tag "cluster_${cluster_id}"
    label 'process_high'
    container "quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0"

    input:
        val cluster_id
        path alignment
        path starting_tree

    output:
        tuple val(cluster_id), path("${cluster_id}.filtered_polymorphic_sites.fasta"), emit: filtered_alignment
        tuple val(cluster_id), path("${cluster_id}.recombination_predictions.gff"), emit: recombination_gff
        tuple val(cluster_id), path("${cluster_id}.node_labelled.final_tree.tre"), emit: final_tree
        path "versions.yml", emit: versions
        path ".diagnostics.log", emit: diagnostics

    script:
    def iterations = params.gubbins_iterations ?: 3
    def tree_builder = params.gubbins_tree_builder ?: 'iqtree'
    def first_tree_builder = params.gubbins_first_tree_builder ?: 'rapidnj'
    def min_snps = params.gubbins_min_snps ?: 2
    def use_hybrid = params.gubbins_use_hybrid ?: true
    def cpus = task.cpus
    """
    set -euo pipefail

    echo "=== Gubbins Diagnostics for cluster ${cluster_id} ===" > .diagnostics.log
    seq_count=\$(grep -c "^>" ${alignment})
    aln_len=\$(awk '/^[^>]/ {sum += length(\$0)} END {print sum}' ${alignment})
    echo "Sequence count: \$seq_count" >> .diagnostics.log
    echo "Alignment length: \$aln_len" >> .diagnostics.log

    # Check if input files exist and are readable
    echo "Input file checks:" >> .diagnostics.log
    echo "Alignment file: ${alignment} - \$(ls -la ${alignment})" >> .diagnostics.log
    echo "Starting tree: ${starting_tree} - \$(ls -la ${starting_tree})" >> .diagnostics.log

    # Logic to omit --first-tree-builder if starting_tree is provided
    if [ -s "${starting_tree}" ]; then
        gubbins_first_tree_builder_arg=""
    else
        gubbins_first_tree_builder_arg="--first-tree-builder ${first_tree_builder}"
    fi

    if [ \$seq_count -lt 3 ]; then
        echo "WARNING: Alignment has only \$seq_count sequences. Skipping Gubbins." >> .diagnostics.log
        touch ${cluster_id}.filtered_polymorphic_sites.fasta
        touch ${cluster_id}.recombination_predictions.gff
        touch ${cluster_id}.node_labelled.final_tree.tre
        echo "Gubbins not run due to insufficient sequences." >> .diagnostics.log
    else
        echo "Running Gubbins with the following parameters:" >> .diagnostics.log
        echo "  iterations: ${iterations}" >> .diagnostics.log
        echo "  tree_builder: ${tree_builder}" >> .diagnostics.log
        echo "  first_tree_builder: ${first_tree_builder}" >> .diagnostics.log
        echo "  min_snps: ${min_snps}" >> .diagnostics.log
        echo "  use_hybrid: ${use_hybrid}" >> .diagnostics.log
        echo "  cpus: ${cpus}" >> .diagnostics.log

        if [ "${use_hybrid}" = "true" ]; then
            echo "Running Gubbins with hybrid mode..." >> .diagnostics.log
            run_gubbins.py \
                --starting-tree ${starting_tree} \
                --prefix ${cluster_id} \
                \${gubbins_first_tree_builder_arg} \
                --tree-builder ${tree_builder} \
                --iterations ${iterations} \
                --min-snps ${min_snps} \
                --threads ${cpus} \
                ${alignment} >> .diagnostics.log 2>&1
            gubbins_exit_code=\$?
        else
            echo "Running Gubbins without hybrid mode..." >> .diagnostics.log
            run_gubbins.py \
                --starting-tree ${starting_tree} \
                --prefix ${cluster_id} \
                \${gubbins_first_tree_builder_arg} \
                --tree-builder ${tree_builder} \
                --iterations ${iterations} \
                --min-snps ${min_snps} \
                --threads ${cpus} \
                ${alignment} >> .diagnostics.log 2>&1
            gubbins_exit_code=\$?
        fi

        echo "Gubbins exit code: \$gubbins_exit_code" >> .diagnostics.log

        if [ \$gubbins_exit_code -ne 0 ]; then
            echo "ERROR: Gubbins failed for cluster ${cluster_id} with exit code \$gubbins_exit_code" >> .diagnostics.log
            echo "Creating empty output files..." >> .diagnostics.log
            touch ${cluster_id}.filtered_polymorphic_sites.fasta
            touch ${cluster_id}.recombination_predictions.gff
            touch ${cluster_id}.node_labelled.final_tree.tre
        fi

        # Check output file sizes and log
        for f in ${cluster_id}.filtered_polymorphic_sites.fasta ${cluster_id}.recombination_predictions.gff ${cluster_id}.node_labelled.final_tree.tre; do
            if [ ! -s "\$f" ]; then
                echo "WARNING: Output file \$f is empty." >> .diagnostics.log
                touch "\$f"
            else
                echo "Output file \$f size: \$(stat -c%s "\$f") bytes" >> .diagnostics.log
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    gubbins: \$(run_gubbins.py --version | sed 's/^/    /')
END_VERSIONS
    """
}