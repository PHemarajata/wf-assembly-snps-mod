process SKA_ALIGN {
    tag "cluster_${cluster_id}"
    label 'process_medium'
    container "quay.io/biocontainers/ska2:0.3.7--h4349ce8_2"

    input:
    tuple val(cluster_id), path(ska_file)

    output:
    tuple val(cluster_id), path("${cluster_id}.aln.fa"), emit: alignment
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create alignment from SKA file
    ska align \\
        $args \\
        $ska_file > ${cluster_id}.aln.fa || {
        echo "WARNING: SKA align failed for cluster ${cluster_id}. Creating empty alignment file."
        touch ${cluster_id}.aln.fa
    }

    # Ensure alignment file exists
    if [ ! -f "${cluster_id}.aln.fa" ]; then
        echo "WARNING: Missing alignment file for cluster ${cluster_id}. Creating empty file."
        touch ${cluster_id}.aln.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ska: \$(ska --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}