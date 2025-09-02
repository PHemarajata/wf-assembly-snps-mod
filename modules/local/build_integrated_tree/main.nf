process BUILD_INTEGRATED_TREE {
    tag "integrated_phylogeny"
    label 'process_high'
    container "quay.io/biocontainers/iqtree:2.2.6--h21ec9f0_0"
    
    publishDir "${params.outdir}/Integrated_Results", mode: params.publish_dir_mode, pattern: "*.{treefile,iqtree,log}"

    input:
    path integrated_alignment
    path sample_mapping

    output:
    path "integrated_core_snps.treefile", emit: tree
    path "integrated_core_snps.iqtree", emit: log
    path "integrated_phylogeny_report.txt", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model = params.iqtree_model ?: 'GTR+ASC'
    """
    echo "Checking integrated alignment file: ${integrated_alignment}"
    
    # Check if alignment file exists and has content
    if [ ! -f "${integrated_alignment}" ]; then
        echo "ERROR: Alignment file ${integrated_alignment} does not exist"
        seq_count=0
    elif [ ! -s "${integrated_alignment}" ]; then
        echo "WARNING: Alignment file ${integrated_alignment} is empty"
        seq_count=0
    else
        # Count sequences safely
        seq_count=\$(grep -c "^>" "${integrated_alignment}" 2>/dev/null || echo "0")
        echo "Found \$seq_count sequences in alignment"
    fi
    
    # Check if we have enough sequences for phylogenetic analysis
    if [ "\$seq_count" -lt 3 ]; then
        echo "WARNING: Integrated alignment has only \$seq_count sequences. Cannot build phylogenetic tree."
        echo "Minimum 3 sequences required for phylogenetic analysis."
        
        # Create empty output files
        touch integrated_core_snps.treefile
        touch integrated_core_snps.iqtree
        
        echo "Insufficient sequences for phylogenetic analysis (\$seq_count sequences found)" > integrated_phylogeny_report.txt
        echo "Minimum 3 sequences required for tree construction" >> integrated_phylogeny_report.txt
        
    else
        echo "Building phylogenetic tree from integrated core SNPs alignment (\$seq_count sequences)"
        
        # Check alignment content
        echo "Alignment file size: \$(wc -c < "${integrated_alignment}") bytes"
        echo "First few lines of alignment:"
        head -10 "${integrated_alignment}" || echo "Could not read alignment file"
        
        # Build tree with IQ-TREE
        iqtree2 \\
            -s ${integrated_alignment} \\
            -st DNA \\
            -m MFP \\
            -bb 1000 \\
            -alrt 1000 \\
            -nt AUTO \\
            --prefix integrated_core_snps \\
            || {
            echo "WARNING: IQ-TREE failed. Creating empty output files."
            touch integrated_core_snps.treefile
            touch integrated_core_snps.iqtree
        }
        
        # Create phylogeny report without pandas (since pip not available in this container)
        python3 << 'EOF'
import os

# Create basic phylogeny report without pandas
try:
    # Read sample mapping manually
    mapping_data = []
    try:
        with open("${sample_mapping}", 'r') as f:
            lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\\t')
            if len(parts) >= 2:
                mapping_data.append({'sample_id': parts[0], 'cluster_id': parts[1]})
    except Exception as e:
        print(f"Warning: Could not read sample mapping: {e}")
    
    with open("integrated_phylogeny_report.txt", 'w') as f:
        f.write("INTEGRATED PHYLOGENETIC ANALYSIS REPORT\\n")
        f.write("=" * 50 + "\\n\\n")
        
        f.write(f"Total samples in phylogeny: {len(mapping_data)}\\n")
        clusters = set([item['cluster_id'] for item in mapping_data])
        f.write(f"Clusters represented: {len(clusters)}\\n\\n")
        
        f.write("Phylogenetic method: IQ-TREE with model selection\\n")
        f.write("Bootstrap support: 1000 replicates\\n")
        f.write("SH-aLRT support: 1000 replicates\\n")
        
        # Check if tree was successfully built
        if os.path.exists("integrated_core_snps.treefile") and os.path.getsize("integrated_core_snps.treefile") > 0:
            f.write("\\nTree construction: SUCCESSFUL\\n")
        else:
            f.write("\\nTree construction: FAILED\\n")
            f.write("Possible reasons:\\n")
            f.write("- Empty or invalid alignment file\\n")
            f.write("- Insufficient variable sites\\n")
            f.write("- All sequences identical\\n")

except Exception as e:
    with open("integrated_phylogeny_report.txt", 'w') as f:
        f.write(f"Error creating phylogeny report: {e}\\n")
EOF
    fi

    # Ensure all output files exist
    for file in integrated_core_snps.treefile integrated_core_snps.iqtree integrated_phylogeny_report.txt; do
        if [ ! -f "\$file" ]; then
            echo "Creating missing file: \$file"
            touch "\$file"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 --version 2>&1 | head -n1 | sed 's/^/    /')
    END_VERSIONS
    """
}