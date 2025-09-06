#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process IQTREE_ASC {
    tag "cluster_${cluster_id}"
    label 'process_high'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"

    input:
    tuple val(cluster_id), path(filtered_snps), val(representative_id)

    output:
    tuple val(cluster_id), path("${cluster_id}.final.treefile"), val(representative_id), emit: final_tree
    tuple val(cluster_id), path("${cluster_id}.final.iqtree"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model = params.iqtree_asc_model ?: 'GTR+ASC'
    def asc_correction = params.iqtree_asc_correction ?: 'lewis'
    """
    echo "Building final ML tree for cluster ${cluster_id} with ASC correction"
    echo "Representative: ${representative_id}"
    
    # Check if alignment has at least 3 sequences and variable sites
    seq_count=\$(grep -c "^>" $filtered_snps)
    
    if [ \$seq_count -lt 3 ]; then
        echo "WARNING: Alignment has only \$seq_count sequences. IQ-TREE requires at least 3 sequences."
        echo "Creating empty output files for cluster ${cluster_id}"
        
        touch ${cluster_id}.final.treefile
        touch ${cluster_id}.final.iqtree
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
        END_VERSIONS
        
        exit 0
    fi
    
    # Check if alignment has variable sites
    if [ ! -s "$filtered_snps" ]; then
        echo "WARNING: Empty or no variable sites in alignment for cluster ${cluster_id}"
        echo "Creating empty output files"
        
        touch ${cluster_id}.final.treefile
        touch ${cluster_id}.final.iqtree
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
        END_VERSIONS
        
        exit 0
    fi
    
    # Count variable sites
    var_sites=\$(awk '/^[^>]/ {for(i=1;i<=length(\$0);i++) chars[i]=chars[i]substr(\$0,i,1)} END {for(i in chars) {split(chars[i],arr,""); asort(arr); if(length(arr)>1) count++} print count+0}' $filtered_snps)
    
    echo "Found \$var_sites variable sites in alignment"
    
    if [ "\$var_sites" -eq 0 ]; then
        echo "WARNING: No variable sites found in filtered alignment for cluster ${cluster_id}"
        echo "Creating minimal tree"
        
        # Create a star tree with all samples
        sample_names=\$(grep "^>" $filtered_snps | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
        echo "(\$sample_names);" > ${cluster_id}.final.treefile
        touch ${cluster_id}.final.iqtree
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
        END_VERSIONS
        
        exit 0
    fi
    
    echo "Running IQ-TREE with ASC correction (${asc_correction}) and model ${model}"
    
    # Run IQ-TREE with ASC correction for SNP-only alignments
    iqtree2 \\
        -s $filtered_snps \\
        -st DNA \\
        -m ${model} \\
        --asc-corr ${asc_correction} \\
        -bb 1000 \\
        -alrt 1000 \\
        -nt AUTO \\
        --prefix ${cluster_id}.final \\
        $args || {
        echo "WARNING: IQ-TREE failed for cluster ${cluster_id}. Creating minimal tree."
        
        # Create a star tree as fallback
        sample_names=\$(grep "^>" $filtered_snps | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
        echo "(\$sample_names);" > ${cluster_id}.final.treefile
        touch ${cluster_id}.final.iqtree
    }
    
    # Ensure output files exist
    if [ ! -f "${cluster_id}.final.treefile" ]; then
        echo "WARNING: Missing final treefile for cluster ${cluster_id}. Creating minimal tree."
        sample_names=\$(grep "^>" $filtered_snps | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
        echo "(\$sample_names);" > ${cluster_id}.final.treefile
    fi
    
    if [ ! -f "${cluster_id}.final.iqtree" ]; then
        touch ${cluster_id}.final.iqtree
    fi
    
    echo "Final ML tree construction completed for cluster ${cluster_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}