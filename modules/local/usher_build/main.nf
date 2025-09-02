process USHER_BUILD {
    tag "$meta.id"
    label 'process_high'
    container "quay.io/biocontainers/usher:0.6.2--h9ee0642_0"

    input:
    tuple val(meta), path(vcf)
    path reference

    output:
    tuple val(meta), path("*.pb"), emit: mat
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Build mutation-annotated tree (MAT) from VCF
    usher \\
        -v $vcf \\
        -t /dev/null \\
        -d $reference \\
        -o ${prefix}.pb \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usher: \$(usher --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}