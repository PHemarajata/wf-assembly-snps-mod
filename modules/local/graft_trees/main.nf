#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GRAFT_TREES {
    tag "tree_grafting"
    label 'process_medium'
    container "quay.io/biocontainers/ete3:3.1.1--py35_0"
    
    publishDir "${params.outdir}/Integrated_Results", mode: params.publish_dir_mode, pattern: "*.{tre,txt,pdf}"

    input:
    path tree_files
    path clusters_file
    path integrated_alignment

    output:
    path "grafted_supertree.tre", emit: supertree
    path "tree_grafting_report.txt", emit: report
    path "cluster_tree_summary.tsv", emit: tree_summary
    path "supertree_visualization.pdf", optional: true, emit: visualization
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install compatible versions of numpy and pandas
    pip install --upgrade numpy>=1.15.4
    pip install pandas

    python3 << 'EOF'
import os
try:
    import pandas as pd
    pandas_available = True
except ImportError:
    pandas_available = False
    print("Warning: pandas not available, using basic functionality")

from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import glob
from collections import defaultdict
import sys

def graft_cluster_trees():
    # Graft cluster trees into a supertree using tree grafting algorithms
    
    # Read cluster assignments
    if pandas_available:
        clusters_df = pd.read_csv("${clusters_file}", sep='\\t')
        sample_to_cluster = dict(zip(clusters_df['sample_id'], clusters_df['cluster_id']))
    else:
        # Basic file reading without pandas
        sample_to_cluster = {}
        with open("${clusters_file}", 'r') as f:
            lines = f.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split('\\t')
            if len(parts) >= 2:
                sample_to_cluster[parts[1]] = parts[0]  # sample_id -> cluster_id
    
    cluster_to_samples = defaultdict(list)
    for sample, cluster in sample_to_cluster.items():
        cluster_to_samples[cluster].append(sample)
    
    # Load all cluster trees
    cluster_trees = {}
    tree_files = glob.glob("*.treefile") + glob.glob("*.tre") + glob.glob("*.tree")
    
    tree_summary_data = []
    
    for tree_file in tree_files:
        try:
            # Extract cluster ID from filename
            cluster_id = None
            for possible_cluster in cluster_to_samples.keys():
                if str(possible_cluster) in tree_file:
                    cluster_id = possible_cluster
                    break
            
            if cluster_id is None:
                print(f"Warning: Could not determine cluster ID for {tree_file}")
                continue
            
            # Load tree
            tree = Tree(tree_file)
            cluster_trees[cluster_id] = tree
            
            # Collect tree statistics
            tree_summary_data.append({
                'cluster_id': cluster_id,
                'tree_file': tree_file,
                'num_leaves': len(tree.get_leaves()),
                'tree_length': tree.get_distance(tree.get_farthest_leaf()[0]),
                'num_internal_nodes': len([n for n in tree.traverse() if not n.is_leaf()])
            })
            
            print(f"Loaded tree for cluster {cluster_id}: {len(tree.get_leaves())} leaves")
            
        except Exception as e:
            print(f"Warning: Could not load tree from {tree_file}: {e}")
    
    if not cluster_trees:
        print("Error: No valid trees found for grafting")
        # Create empty outputs
        with open("grafted_supertree.tre", 'w') as f:
            f.write("();\\n")
        with open("tree_grafting_report.txt", 'w') as f:
            f.write("No trees available for grafting\\n")
        
        if pandas_available:
            pd.DataFrame().to_csv("cluster_tree_summary.tsv", sep='\\t', index=False)
        else:
            with open("cluster_tree_summary.tsv", 'w') as f:
                f.write("cluster_id\\ttree_file\\tnum_leaves\\ttree_length\\tnum_internal_nodes\\n")
        return
    
    # Create supertree using simple grafting approach
    print(f"Grafting {len(cluster_trees)} cluster trees...")
    
    # Method 1: Create a star tree with cluster representatives, then graft subtrees
    supertree = Tree()
    supertree.name = "root"
    
    # Add each cluster as a subtree
    for cluster_id, cluster_tree in cluster_trees.items():
        # Create a copy of the cluster tree
        cluster_subtree = cluster_tree.copy()
        
        # Add cluster label to internal nodes
        for node in cluster_subtree.traverse():
            if hasattr(node, 'name') and node.name:
                node.name = f"{node.name}_cluster_{cluster_id}"
            else:
                node.name = f"cluster_{cluster_id}_node"
        
        # Attach to supertree
        supertree.add_child(cluster_subtree)
    
    # Write supertree
    supertree.write(format=1, outfile="grafted_supertree.tre")
    
    # Create tree summary
    if pandas_available:
        tree_summary_df = pd.DataFrame(tree_summary_data)
        tree_summary_df.to_csv("cluster_tree_summary.tsv", sep='\\t', index=False)
    else:
        with open("cluster_tree_summary.tsv", 'w') as f:
            f.write("cluster_id\\ttree_file\\tnum_leaves\\ttree_length\\tnum_internal_nodes\\n")
            for data in tree_summary_data:
                f.write(f"{data['cluster_id']}\\t{data['tree_file']}\\t{data['num_leaves']}\\t{data['tree_length']:.6f}\\t{data['num_internal_nodes']}\\n")
    
    # Create detailed report
    with open("tree_grafting_report.txt", 'w') as f:
        f.write("TREE GRAFTING ANALYSIS REPORT\\n")
        f.write("=" * 50 + "\\n\\n")
        
        f.write(f"Total cluster trees processed: {len(cluster_trees)}\\n")
        f.write(f"Total leaves in supertree: {len(supertree.get_leaves())}\\n")
        f.write(f"Supertree depth: {supertree.get_farthest_leaf()[1]}\\n\\n")
        
        f.write("Cluster tree statistics:\\n")
        f.write("-" * 30 + "\\n")
        
        for data in tree_summary_data:
            f.write(f"Cluster {data['cluster_id']}:\\n")
            f.write(f"  - Leaves: {data['num_leaves']}\\n")
            f.write(f"  - Internal nodes: {data['num_internal_nodes']}\\n")
            f.write(f"  - Tree length: {data['tree_length']:.6f}\\n")
            f.write(f"  - Source file: {data['tree_file']}\\n\\n")
        
        f.write("Grafting method: Star topology with cluster subtrees\\n")
        f.write("Note: Each cluster forms a monophyletic group in the supertree\\n")
    
    # Create visualization (optional, may fail if display not available)
    try:
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        ts.mode = "c"  # circular
        ts.arc_start = -180
        ts.arc_span = 180
        
        # Color nodes by cluster
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        cluster_colors = {}
        
        for i, cluster_id in enumerate(cluster_trees.keys()):
            cluster_colors[cluster_id] = colors[i % len(colors)]
        
        for node in supertree.traverse():
            if "cluster_" in str(node.name):
                for cluster_id, color in cluster_colors.items():
                    if f"cluster_{cluster_id}" in str(node.name):
                        nstyle = NodeStyle()
                        nstyle["fgcolor"] = color
                        nstyle["size"] = 10
                        node.set_style(nstyle)
                        break
        
        supertree.render("supertree_visualization.pdf", tree_style=ts)
        print("Created supertree visualization")
        
    except Exception as e:
        print(f"Could not create visualization: {e}")
    
    print(f"Tree grafting complete: {len(cluster_trees)} clusters grafted into supertree")

# Run tree grafting
graft_cluster_trees()
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        ete3: \$(python -c "import ete3; print(ete3.__version__)")
        pandas: \$(python -c "try: import pandas; print(pandas.__version__); except: print('not available')")
        numpy: \$(python -c "try: import numpy; print(numpy.__version__); except: print('not available')")
    END_VERSIONS
    """
}