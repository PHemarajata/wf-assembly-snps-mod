process SNP_SITES {
    tag "$meta.id"
    label 'process_low'
    container "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.phylip"), emit: phylip, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snp-sites \\
        $args \\
        -v \\
        -o ${prefix}.vcf \\
        $alignment

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-sites: \$(snp-sites -V 2>&1 | sed 's/snp-sites //')
    END_VERSIONS
    """
}