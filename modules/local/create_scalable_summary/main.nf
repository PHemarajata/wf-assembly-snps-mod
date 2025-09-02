process CREATE_SCALABLE_SUMMARY {
    tag "scalable_summary"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Summaries", mode: params.publish_dir_mode, pattern: "*.{tsv,txt}"

    input:
    path clusters_file
    path cluster_summary
    path alignments
    path trees

    output:
    path "Scalable_Analysis_Summary.tsv", emit: summary
    path "Scalable_Analysis_Report.txt", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install required packages
    pip install pandas

    # Create summary
    python3 << 'EOF'
import pandas as pd
import os
import glob

# Read cluster assignments
clusters_df = pd.read_csv("${clusters_file}", sep='\\t')
total_samples = len(clusters_df)
total_clusters = clusters_df['cluster_id'].nunique()

# Count files in each category
alignment_files = glob.glob("*.aln.fa")
tree_files = glob.glob("*.treefile")

# Create summary dataframe
summary_data = {
    'Metric': [
        'Total Samples',
        'Total Clusters',
        'Clusters with Alignments',
        'Clusters with Trees',
        'Average Samples per Cluster',
        'Largest Cluster Size',
        'Smallest Cluster Size'
    ],
    'Value': [
        total_samples,
        total_clusters,
        len(alignment_files),
        len(tree_files),
        round(total_samples / total_clusters, 2) if total_clusters > 0 else 0,
        clusters_df.groupby('cluster_id').size().max() if total_clusters > 0 else 0,
        clusters_df.groupby('cluster_id').size().min() if total_clusters > 0 else 0
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('Scalable_Analysis_Summary.tsv', sep='\\t', index=False)

# Create detailed report
with open('Scalable_Analysis_Report.txt', 'w') as f:
    f.write("SCALABLE ASSEMBLY SNPs ANALYSIS REPORT\\n")
    f.write("=" * 50 + "\\n\\n")
    
    f.write(f"Total samples processed: {total_samples}\\n")
    f.write(f"Number of clusters formed: {total_clusters}\\n")
    f.write(f"Clusters with alignments: {len(alignment_files)}\\n")
    f.write(f"Clusters with phylogenetic trees: {len(tree_files)}\\n\\n")
    
    if total_clusters > 0:
        cluster_sizes = clusters_df.groupby('cluster_id').size()
        f.write("CLUSTER SIZE DISTRIBUTION:\\n")
        f.write("-" * 30 + "\\n")
        for cluster_id, size in cluster_sizes.items():
            f.write(f"Cluster {cluster_id}: {size} samples\\n")
        
        f.write(f"\\nAverage cluster size: {cluster_sizes.mean():.2f}\\n")
        f.write(f"Largest cluster: {cluster_sizes.max()} samples\\n")
        f.write(f"Smallest cluster: {cluster_sizes.min()} samples\\n")
    
    f.write("\\nOUTPUT DIRECTORIES:\\n")
    f.write("-" * 20 + "\\n")
    f.write("- Clusters/: Individual cluster results\\n")
    f.write("- Summaries/: Clustering summaries and reports\\n")
    f.write("- Mash_Sketches/: Mash sketch files\\n")

print("Summary files created successfully!")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}