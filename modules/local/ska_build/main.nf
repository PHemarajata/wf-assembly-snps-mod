process SKA_BUILD {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/ska2:0.3.7--h4349ce8_2"

    input:
    tuple val(cluster_id), val(sample_ids), path(assemblies)

    output:
    tuple val(cluster_id), path("${cluster_id}.skf"), emit: ska_file
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Create input content using staged file names
    def input_content = [sample_ids, assemblies].transpose().collect{ sample_id, assembly -> 
        "${sample_id}\t${assembly.name}" 
    }.join('\n')
    """
    # Create input file for SKA
    cat > ${cluster_id}_input.tsv <<'EOFSKA'
${input_content}
EOFSKA

    # Verify input file was created correctly
    echo "Input file contents:"
    cat ${cluster_id}_input.tsv
    
    # Verify that all files exist
    echo "Checking file existence:"
    for file in *.fasta *.fa *.fna *.fas *.fsa; do
        if [ -f "\$file" ]; then
            echo "Found: \$file"
        fi
    done 2>/dev/null || echo "No FASTA files found with standard extensions"
    
    # Build SKA file
    ska build \\
        -o ${cluster_id} \\
        -f ${cluster_id}_input.tsv \\
        ${args} || {
        echo "WARNING: SKA build failed for cluster ${cluster_id}. Creating empty SKA file."
        touch ${cluster_id}.skf
    }

    # Ensure SKA file exists
    if [ ! -f "${cluster_id}.skf" ]; then
        echo "WARNING: Missing SKA file for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.skf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}