process USHER_PLACE {
    tag "$meta.id"
    label 'process_medium'
    container "quay.io/biocontainers/usher:0.6.2--h9ee0642_0"

    input:
    tuple val(meta), path(mat)
    path new_vcf

    output:
    tuple val(meta), path("*.pb"), emit: updated_mat
    tuple val(meta), path("*.txt"), emit: placement_stats, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_updated"
    """
    # Place new samples on existing MAT
    usher \\
        -i $mat \\
        -v $new_vcf \\
        -o ${prefix}.pb \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usher: \$(usher --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}