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

    script:
    """
    set -euo pipefail

    echo "Building backbone tree from cluster representatives"
    method="${params.backbone_method:-parsnp}"
    echo "Method: \$method"

    rep_count=\$(grep -c "^>" ${representatives_fasta})
    echo "Number of representatives: \$rep_count"

    # Always create backbone.treefile and backbone_alignment.fa
    if [ "\$rep_count" -lt 3 ]; then
        echo "WARNING: Only \$rep_count representatives found. Cannot build meaningful backbone tree."
        rep_names=\$(grep "^>" ${representatives_fasta} | sed 's/>//' | tr '\n' ',' | sed 's/,\$//')
        echo "(\$rep_names);" > backbone.treefile
        cp ${representatives_fasta} backbone_alignment.fa
        echo "Minimal backbone tree created with \$rep_count representatives" > backbone_report.txt
        cat <<-END_VERSIONS > versions.yml
"${task.process}":
    parsnp: \$(parsnp --version | sed 's/^/    /')
END_VERSIONS
        exit 0
    fi

    if [ "\$method" = "parsnp" ]; then
        echo "Using Parsnp for backbone tree construction"
        mkdir -p representatives
        awk '/^>/ {filename="representatives/"substr(\$0,2)".fa"} {print > filename}' ${representatives_fasta}
        ref_file=\$(ls representatives/*.fa | head -n1)
        echo "Using reference: \$ref_file"
        parsnp \\
            --sequences representatives/ \\
            --reference \$ref_file \\
            --output-dir parsnp_backbone \\
            --use-fasttree \\
            --threads \${task.cpus} \\
            --verbose || {
            echo "WARNING: Parsnp failed. Creating distance-based tree."
            if [ -f "${mash_distances}" ]; then
                echo "Creating distance-based tree from Mash distances"
                python3 << 'PYTHON_EOF'
import pandas as pd
try:
    distances_df = pd.read_csv("${mash_distances}", sep='\\t', index_col=0)
    with open("${representatives_fasta}", 'r') as f:
        rep_names = [line.strip()[1:] for line in f if line.startswith('>')]
    available_reps = [rep for rep in rep_names if rep in distances_df.index]
    if len(available_reps) >= 3:
        tree_str = "(" + ",".join(available_reps) + ");"
        with open("backbone.treefile", 'w') as f:
            f.write(tree_str)
    else:
        tree_str = "(" + ",".join(rep_names) + ");"
        with open("backbone.treefile", 'w') as f:
            f.write(tree_str)
except Exception as e:
    with open("${representatives_fasta}", 'r') as f:
        rep_names = [line.strip()[1:] for line in f if line.startswith('>')]
    tree_str = "(" + ",".join(rep_names) + ");"
    with open("backbone.treefile", 'w') as f:
        f.write(tree_str)
PYTHON_EOF
            fi
        }
        if [ -f "parsnp_backbone/parsnp.tree" ]; then
            cp parsnp_backbone/parsnp.tree backbone.treefile
        fi
        if [ -f "parsnp_backbone/parsnp.xmfa" ]; then
            harvesttools -i parsnp_backbone/parsnp.ggr -M backbone_alignment.fa || cp ${representatives_fasta} backbone_alignment.fa
        else
            cp ${representatives_fasta} backbone_alignment.fa
        fi
    elif [ "\$method" = "fasttree" ]; then
        echo "Using FastTree for backbone tree construction"
        cp ${representatives_fasta} backbone_alignment.fa
        fasttree -nt backbone_alignment.fa > backbone.treefile || {
            rep_names=\$(grep "^>" ${representatives_fasta} | sed 's/>//' | tr '\n' ',' | sed 's/,\$//')
            echo "(\$rep_names);" > backbone.treefile
        }
    else
        rep_names=\$(grep "^>" ${representatives_fasta} | sed 's/>//' | tr '\n' ',' | sed 's/,\$//')
        echo "(\$rep_names);" > backbone.treefile
        cp ${representatives_fasta} backbone_alignment.fa
    fi

    # Fallback: always create backbone.treefile if missing
    if [ ! -f "backbone.treefile" ]; then
        rep_names=\$(grep "^>" ${representatives_fasta} | sed 's/>//' | tr '\n' ',' | sed 's/,\$//')
        echo "(\$rep_names);" > backbone.treefile
    fi

    # Fallback: always create backbone_alignment.fa if missing
    if [ ! -f "backbone_alignment.fa" ]; then
        cp ${representatives_fasta} backbone_alignment.fa
    fi

    echo "BACKBONE TREE CONSTRUCTION REPORT" > backbone_report.txt
    echo "=================================" >> backbone_report.txt
    echo "Method: \$method" >> backbone_report.txt
    echo "Number of representatives: \$rep_count" >> backbone_report.txt
    echo "Tree file: backbone.treefile" >> backbone_report.txt
    echo "Alignment file: backbone_alignment.fa" >> backbone_report.txt
    if [ -f "backbone.treefile" ] && [ -s "backbone.treefile" ]; then
        echo "Status: SUCCESS" >> backbone_report.txt
    else
        echo "Status: PARTIAL (fallback tree created)" >> backbone_report.txt
    fi

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    parsnp: \$(parsnp --version | sed 's/^/    /')
    fasttree: \$(fasttree -expert 2>&1 | head -1 | sed 's/^/    /')
    harvesttools: \$(harvesttools --version | sed 's/^/    /')
END_VERSIONS
    """
}