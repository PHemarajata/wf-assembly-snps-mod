process GUBBINS_CLUSTER {
    tag "cluster_${cluster_id}"
    label 'process_high'
    container "quay.io/biocontainers/gubbins:3.3.5--py39pl5321he4a0461_0"

    input:
    tuple val(cluster_id), path(alignment), path(starting_tree)

    output:
    tuple val(cluster_id), path("${cluster_id}.filtered_polymorphic_sites.fasta"), emit: filtered_alignment
    tuple val(cluster_id), path("${cluster_id}.recombination_predictions.gff"), emit: recombination_gff
    tuple val(cluster_id), path("${cluster_id}.node_labelled.final_tree.tre"), emit: final_tree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def iterations = params.gubbins_iterations ?: 3
    def tree_builder = params.gubbins_tree_builder ?: 'iqtree'
    def first_tree_builder = params.gubbins_first_tree_builder ?: 'rapidnj'
    def min_snps = params.gubbins_min_snps ?: 5
    def use_hybrid = params.gubbins_use_hybrid ?: true
    """
    # Check if alignment has at least 3 sequences (minimum for phylogenetic analysis)
    seq_count=\$(grep -c "^>" $alignment)
    
    if [ \$seq_count -lt 3 ]; then
        echo "WARNING: Alignment has only \$seq_count sequences. Gubbins requires at least 3 sequences for phylogenetic analysis."
        echo "Skipping Gubbins analysis for cluster ${cluster_id}"
        
        # Create empty output files to satisfy pipeline expectations
        touch ${cluster_id}.filtered_polymorphic_sites.fasta
        touch ${cluster_id}.recombination_predictions.gff
        touch ${cluster_id}.node_labelled.final_tree.tre
        
        # Create versions file for small clusters
        echo '"${task.process}":' > versions.yml
        echo '    gubbins: '\$(run_gubbins.py --version | sed 's/^/    /') >> versions.yml
        
        exit 0
    fi

    # Build Gubbins command with hybrid tree builders if enabled
    if [ "$use_hybrid" = "true" ]; then
        # Use hybrid approach with two tree builders
        run_gubbins.py \\
            --starting-tree $starting_tree \\
            --prefix ${cluster_id} \\
            --first-tree-builder $first_tree_builder \\
            --tree-builder $tree_builder \\
            --iterations $iterations \\
            --min-snps $min_snps \\
            --threads ${task.cpus} \\
            $args \\
            $alignment || {
            echo "WARNING: Gubbins failed for cluster ${cluster_id}. Creating empty output files."
            touch ${cluster_id}.filtered_polymorphic_sites.fasta
            touch ${cluster_id}.recombination_predictions.gff
            touch ${cluster_id}.node_labelled.final_tree.tre
        }
    else
        # Use single tree builder
        run_gubbins.py \\
            --starting-tree $starting_tree \\
            --prefix ${cluster_id} \\
            --tree-builder $tree_builder \\
            --iterations $iterations \\
            --min-snps $min_snps \\
            --threads ${task.cpus} \\
            $args \\
            $alignment || {
            echo "WARNING: Gubbins failed for cluster ${cluster_id}. Creating empty output files."
            touch ${cluster_id}.filtered_polymorphic_sites.fasta
            touch ${cluster_id}.recombination_predictions.gff
            touch ${cluster_id}.node_labelled.final_tree.tre
        }
    fi

    # Ensure all required output files exist (in case Gubbins partially failed)
    if [ ! -f "${cluster_id}.filtered_polymorphic_sites.fasta" ]; then
        echo "WARNING: Missing filtered_polymorphic_sites.fasta for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.filtered_polymorphic_sites.fasta
    fi
    
    if [ ! -f "${cluster_id}.recombination_predictions.gff" ]; then
        echo "WARNING: Missing recombination_predictions.gff for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.recombination_predictions.gff
    fi
    
    if [ ! -f "${cluster_id}.node_labelled.final_tree.tre" ]; then
        echo "WARNING: Missing node_labelled.final_tree.tre for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.node_labelled.final_tree.tre
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version | sed 's/^/    /')
    END_VERSIONS
    """
}