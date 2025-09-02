process IQTREE_FAST {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    input:
    tuple val(cluster_id), path(alignment)

    output:
    tuple val(cluster_id), path("${cluster_id}.treefile"), emit: tree
    tuple val(cluster_id), path("${cluster_id}.iqtree"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model = params.iqtree_model ?: 'GTR+ASC'
    """
    # Check if alignment has at least 3 sequences (minimum for tree building)
    seq_count=\$(grep -c "^>" $alignment)
    
    if [ \$seq_count -lt 3 ]; then
        echo "WARNING: Alignment has only \$seq_count sequences. IQ-TREE requires at least 3 sequences for tree building."
        echo "Skipping tree construction for cluster ${cluster_id}"
        
        # Create empty output files to satisfy pipeline expectations
        touch ${cluster_id}.treefile
        touch ${cluster_id}.iqtree
        
        # Create versions file for small clusters
        echo '"${task.process}":' > versions.yml
        echo '    iqtree: '\$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /') >> versions.yml
        
        exit 0
    fi

    # Run IQ-TREE with fast mode
    iqtree2 \\
        -s $alignment \\
        -st DNA \\
        -m MFP \\
        --fast \\
        -nt AUTO \\
        --prefix ${cluster_id} \\
        $args || {
        echo "WARNING: IQ-TREE failed for cluster ${cluster_id}. Creating empty output files."
        touch ${cluster_id}.treefile
        touch ${cluster_id}.iqtree
    }

    # Ensure all required output files exist (in case IQ-TREE partially failed)
    if [ ! -f "${cluster_id}.treefile" ]; then
        echo "WARNING: Missing treefile for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.treefile
    fi
    
    if [ ! -f "${cluster_id}.iqtree" ]; then
        echo "WARNING: Missing iqtree log for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.iqtree
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}