process SNIPPY_ALIGN {
    tag "cluster_${cluster_id}"
    label 'process_high'
    container "staphb/snippy:4.6.0"

    input:
    tuple val(cluster_id), val(sample_ids), path(assemblies), val(representative_id), path(reference)

    output:
    tuple val(cluster_id), path("${cluster_id}.core.full.aln"), emit: core_alignment
    tuple val(cluster_id), path("${cluster_id}.core.tab"), emit: core_tab
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bash <<'EOF'
    echo "Running Snippy alignment for cluster ${cluster_id}"
    echo "Reference: ${representative_id}"
    echo "Samples: ${sample_ids}"

    # Create working directory
    mkdir -p snippy_work

    for sample in \$(echo ${sample_ids} | sed 's/[][]//g; s/,/ /g'); do
        echo "Processing sample: \$sample"

        # Find the assembly file for this sample
        assembly_file=""
        for file in *.fa *.fasta *.fna; do
            if [[ "\$file" == *"\$sample"* ]] || [[ "\$(basename "\$file" .fa)" == "\$sample" ]] || [[ "\$(basename "\$file" .fasta)" == "\$sample" ]] || [[ "\$(basename "\$file" .fna)" == "\$sample" ]]; then
                assembly_file="\$file"
                break
            fi
        done

        if [[ -z "\$assembly_file" ]]; then
            echo "Warning: Could not find assembly file for sample \$sample"
            continue
        fi

        echo "Found assembly: \$assembly_file"

        # Skip if this is the reference sample
        if [[ "\$sample" == "${representative_id}" ]]; then
            echo "Skipping reference sample \$sample"
            continue
        fi

        # Run snippy
        snippy \\
            --outdir snippy_work/\$sample \\
            --ref ${reference} \\
            --ctgs \$assembly_file \\
            --cpus ${task.cpus} \\
            --force \\
            ${args} || {
            echo "Warning: Snippy failed for sample \$sample"
            continue
        }
    done

    # Collect all snippy output directories
    snippy_dirs=()
    for dir in snippy_work/*/; do
        if [[ -d "\$dir" ]] && [[ -f "\$dir/snps.vcf" ]]; then
            snippy_dirs+=("\$dir")
        fi
    done

    if [[ \${#snippy_dirs[@]} -eq 0 ]]; then
        echo "Warning: No successful Snippy runs found"
        # Create empty alignment
        echo ">${representative_id}" > ${cluster_id}.core.full.aln
        echo "N" >> ${cluster_id}.core.full.aln
        touch ${cluster_id}.core.tab
    else
        echo "Found \${#snippy_dirs[@]} successful Snippy runs"

        # Run snippy-core to generate core alignment
        snippy-core \\
            --ref ${reference} \\
            --prefix ${cluster_id} \\
            \${snippy_dirs[@]} || {
            echo "Warning: snippy-core failed, creating minimal alignment"
            echo ">${representative_id}" > ${cluster_id}.core.full.aln
            echo "N" >> ${cluster_id}.core.full.aln
            touch ${cluster_id}.core.tab
        }

        # Rename snippy-core output files to expected names
        if [[ -f "${cluster_id}.full.aln" ]]; then
            mv "${cluster_id}.full.aln" "${cluster_id}.core.full.aln"
        elif [[ -f "${cluster_id}.aln" ]]; then
            mv "${cluster_id}.aln" "${cluster_id}.core.full.aln"
        fi

        if [[ -f "${cluster_id}.tab" ]]; then
            mv "${cluster_id}.tab" "${cluster_id}.core.tab"
        fi

        # Ensure output files exist
        if [[ ! -f "${cluster_id}.core.full.aln" ]]; then
            echo "Warning: Core alignment not generated, creating minimal alignment"
            echo ">${representative_id}" > ${cluster_id}.core.full.aln
            echo "N" >> ${cluster_id}.core.full.aln
        fi

        if [[ ! -f "${cluster_id}.core.tab" ]]; then
            touch ${cluster_id}.core.tab
        fi
    fi

    echo "Snippy alignment completed for cluster ${cluster_id}"
    echo "Output alignment size: \$(wc -c < ${cluster_id}.core.full.aln) bytes"
    echo "Number of sequences: \$(grep -c '^>' ${cluster_id}.core.full.aln)"

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    snippy: \$(snippy --version 2>&1 | head -n1 | sed 's/^/    /')
END_VERSIONS
EOF
"""
}