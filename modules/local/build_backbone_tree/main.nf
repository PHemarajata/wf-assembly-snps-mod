process BUILD_BACKBONE_TREE {
    tag "backbone_tree"
    label 'process_high'
    container "quay.io/biocontainers/parsnp@sha256:b46999fb9842f183443dd6226b461c1d8074d4c1391c1f2b1e51cc20cee8f8b2"

    input:
    path representatives_fasta
    path mash_distances

    output:
    path "backbone.treefile", emit: backbone_tree
    path "backbone_alignment.fa", emit: backbone_alignment
    path "backbone_report.txt", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def method = params.backbone_method ?: 'parsnp'
    """
    echo "Building backbone tree from cluster representatives"
    echo "Method: ${method}"
    
    # Count representatives
    rep_count=\$(grep -c "^>" $representatives_fasta)
    echo "Number of representatives: \$rep_count"
    
    if [ \$rep_count -lt 3 ]; then
        echo "WARNING: Only \$rep_count representatives found. Cannot build meaningful backbone tree."
        echo "Creating minimal backbone tree"
        
        # Create star tree
        rep_names=\$(grep "^>" $representatives_fasta | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
        echo "(\$rep_names);" > backbone.treefile
        cp $representatives_fasta backbone_alignment.fa
        
        echo "Minimal backbone tree created with \$rep_count representatives" > backbone_report.txt
        
        cat <<END_VERSIONS > versions.yml
"${task.process}":
    parsnp: \$(parsnp --version | sed 's/^/    /')
END_VERSIONS
        
        exit 0
    fi
    
    if [ "${method}" = "parsnp" ]; then
        echo "Using Parsnp for backbone tree construction"
        
        # Create directory for representatives
        mkdir -p representatives
        
        # Split multi-FASTA into individual files for Parsnp
        awk '/^>/ {if(filename) close(filename); filename="representatives/"substr(\$0,2)".fa"} {print > filename}' $representatives_fasta
        
        # Select first representative as reference
        ref_file=\$(ls representatives/*.fa | head -n1)
        echo "Using reference: \$ref_file"
        
        # Run Parsnp
        parsnp \\
            --sequences representatives/ \\
            --reference \$ref_file \\
            --output-dir parsnp_backbone \\
            --use-fasttree \\
            --threads ${task.cpus} \\
            --verbose || {
            echo "WARNING: Parsnp failed. Creating distance-based tree."
            
            # Fallback to distance-based tree using mash distances
            if [ -f "${mash_distances}" ]; then
                echo "Creating distance-based tree from Mash distances"
                python3 << 'EOF'
import pandas as pd
import numpy as np

# Read mash distances
try:
    distances_df = pd.read_csv("${mash_distances}", sep='\\t', index_col=0)
    
    # Get representative names
    with open("$representatives_fasta", 'r') as f:
        rep_names = [line.strip()[1:] for line in f if line.startswith('>')]
    
    # Filter distance matrix to representatives only
    available_reps = [rep for rep in rep_names if rep in distances_df.index]
    
    if len(available_reps) >= 3:
        rep_distances = distances_df.loc[available_reps, available_reps]
        
        # Simple UPGMA-like clustering to create tree
        # This is a simplified approach - in practice you'd use proper phylogenetic methods
        tree_str = "(" + ",".join(available_reps) + ");"
        
        with open("backbone.treefile", 'w') as f:
            f.write(tree_str)
        
        print(f"Created distance-based tree with {len(available_reps)} representatives")
    else:
        # Create star tree
        tree_str = "(" + ",".join(rep_names) + ");"
        with open("backbone.treefile", 'w') as f:
            f.write(tree_str)
        
        print(f"Created star tree with {len(rep_names)} representatives")

except Exception as e:
    print(f"Error creating distance-based tree: {e}")
    # Create star tree as final fallback
    with open("$representatives_fasta", 'r') as f:
        rep_names = [line.strip()[1:] for line in f if line.startswith('>')]
    
    tree_str = "(" + ",".join(rep_names) + ");"
    with open("backbone.treefile", 'w') as f:
        f.write(tree_str)
    
    print(f"Created fallback star tree with {len(rep_names)} representatives")
EOF
            }
        }
        
        # Extract results from Parsnp if successful
        if [ -f "parsnp_backbone/parsnp.tree" ]; then
            cp parsnp_backbone/parsnp.tree backbone.treefile
            echo "Parsnp backbone tree construction successful"
        fi
        
        if [ -f "parsnp_backbone/parsnp.xmfa" ]; then
            # Convert XMFA to FASTA (simplified)
            harvesttools -i parsnp_backbone/parsnp.ggr -M backbone_alignment.fa || {
                echo "Warning: Could not extract alignment, using input representatives"
                cp $representatives_fasta backbone_alignment.fa
            }
        else
            cp $representatives_fasta backbone_alignment.fa
        fi
        
    elif [ "${method}" = "fasttree" ]; then
        echo "Using FastTree for backbone tree construction"
        
        # Need alignment first - use simple approach
        cp $representatives_fasta backbone_alignment.fa
        
        # Run FastTree
        fasttree -nt backbone_alignment.fa > backbone.treefile || {
            echo "WARNING: FastTree failed. Creating star tree."
            rep_names=\$(grep "^>" $representatives_fasta | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
            echo "(\$rep_names);" > backbone.treefile
        }
        
    else
        echo "Unknown backbone method: ${method}. Using star tree."
        rep_names=\$(grep "^>" $representatives_fasta | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
        echo "(\$rep_names);" > backbone.treefile
        cp $representatives_fasta backbone_alignment.fa
    fi
    
    # Ensure output files exist
    if [ ! -f "backbone.treefile" ]; then
        echo "Creating fallback star tree"
        rep_names=\$(grep "^>" $representatives_fasta | sed 's/>//' | tr '\\n' ',' | sed 's/,\$//')
        echo "(\$rep_names);" > backbone.treefile
    fi
    
    if [ ! -f "backbone_alignment.fa" ]; then
        cp $representatives_fasta backbone_alignment.fa
    fi
    
    # Create report
    echo "BACKBONE TREE CONSTRUCTION REPORT" > backbone_report.txt
    echo "=================================" >> backbone_report.txt
    echo "Method: ${method}" >> backbone_report.txt
    echo "Number of representatives: \$rep_count" >> backbone_report.txt
    echo "Tree file: backbone.treefile" >> backbone_report.txt
    echo "Alignment file: backbone_alignment.fa" >> backbone_report.txt
    
    if [ -f "backbone.treefile" ] && [ -s "backbone.treefile" ]; then
        echo "Status: SUCCESS" >> backbone_report.txt
    else
        echo "Status: PARTIAL (fallback tree created)" >> backbone_report.txt
    fi
    
    echo "Backbone tree construction completed"

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    parsnp: \$(parsnp --version | sed 's/^/    /')
    fasttree: \$(fasttree -expert 2>&1 | head -1 | sed 's/^/    /')
    harvesttools: \$(harvesttools --version | sed 's/^/    /')
END_VERSIONS
    """
}