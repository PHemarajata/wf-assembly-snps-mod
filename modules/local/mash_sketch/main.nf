process MASH_SKETCH {
    tag "$sample_id"
    label 'process_low'
    container "quay.io/biocontainers/mash:2.3--he348c14_1"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.msh"), emit: sketch
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sketch_size = params.mash_sketch_size ?: 1000
    def kmer_size = params.mash_kmer_size ?: 21
    def min_copies = params.mash_min_copies ?: 1
    """
    mash sketch \\
        -o ${sample_id} \\
        -s ${sketch_size} \\
        -k ${kmer_size} \\
        -m ${min_copies} \\
        $args \\
        $assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1 | sed 's/^/    /')
    END_VERSIONS
    """
}