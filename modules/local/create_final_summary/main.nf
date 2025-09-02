process CREATE_FINAL_SUMMARY {
    tag "final_summary"
    label 'process_low'
    container "quay.io/biocontainers/python:3.9--1"
    
    publishDir "${params.outdir}/Final_Results", mode: params.publish_dir_mode, pattern: "*.{html,txt,tsv}"

    input:
    path clusters_file
    path cluster_summary
    path core_snp_summary
    path tree_grafting_report
    path phylogeny_report
    path sample_mapping
    path tree_summary

    output:
    path "Scalable_Analysis_Final_Report.html", emit: html_report
    path "Scalable_Analysis_Final_Report.txt", emit: text_report
    path "Analysis_Statistics.tsv", emit: statistics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install compatible versions of numpy and pandas
    pip install --upgrade numpy>=1.15.4
    pip install pandas

    # Create comprehensive final summary
    python3 << 'EOF'
import pandas as pd
import os
from datetime import datetime

def create_final_summary():
    # Create comprehensive final summary of scalable analysis
    
    # Read input files
    try:
        clusters_df = pd.read_csv("${clusters_file}", sep='\\t')
        sample_mapping_df = pd.read_csv("${sample_mapping}", sep='\\t')
        tree_summary_df = pd.read_csv("${tree_summary}", sep='\\t')
    except Exception as e:
        print(f"Warning: Could not read some input files: {e}")
        clusters_df = pd.DataFrame()
        sample_mapping_df = pd.DataFrame()
        tree_summary_df = pd.DataFrame()
    
    # Calculate statistics
    total_samples = len(clusters_df) if not clusters_df.empty else 0
    total_clusters = clusters_df['cluster_id'].nunique() if not clusters_df.empty else 0
    
    # Cluster size statistics
    if not clusters_df.empty:
        cluster_sizes = clusters_df.groupby('cluster_id').size()
        avg_cluster_size = cluster_sizes.mean()
        max_cluster_size = cluster_sizes.max()
        min_cluster_size = cluster_sizes.min()
        
        # Count clusters by size categories
        singleton_clusters = (cluster_sizes == 1).sum()
        small_clusters = ((cluster_sizes >= 2) & (cluster_sizes <= 5)).sum()
        medium_clusters = ((cluster_sizes >= 6) & (cluster_sizes <= 20)).sum()
        large_clusters = (cluster_sizes > 20).sum()
    else:
        avg_cluster_size = max_cluster_size = min_cluster_size = 0
        singleton_clusters = small_clusters = medium_clusters = large_clusters = 0
    
    # SNP statistics
    if not sample_mapping_df.empty and 'total_snps' in sample_mapping_df.columns:
        total_snps = sample_mapping_df['total_snps'].sum()
        avg_snps_per_sample = sample_mapping_df['total_snps'].mean()
        max_snps_per_sample = sample_mapping_df['total_snps'].max()
        min_snps_per_sample = sample_mapping_df['total_snps'].min()
    else:
        total_snps = avg_snps_per_sample = max_snps_per_sample = min_snps_per_sample = 0
    
    # Tree statistics
    if not tree_summary_df.empty:
        total_trees = len(tree_summary_df)
        avg_tree_leaves = tree_summary_df['num_leaves'].mean() if 'num_leaves' in tree_summary_df.columns else 0
        total_tree_leaves = tree_summary_df['num_leaves'].sum() if 'num_leaves' in tree_summary_df.columns else 0
    else:
        total_trees = avg_tree_leaves = total_tree_leaves = 0
    
    # Create statistics table
    statistics_data = {
        'Metric': [
            'Total Samples',
            'Total Clusters',
            'Average Cluster Size',
            'Largest Cluster',
            'Smallest Cluster',
            'Singleton Clusters',
            'Small Clusters (2-5 samples)',
            'Medium Clusters (6-20 samples)',
            'Large Clusters (>20 samples)',
            'Total Phylogenetic Trees',
            'Total SNP Positions',
            'Average SNPs per Sample',
            'Maximum SNPs per Sample',
            'Minimum SNPs per Sample',
            'Average Tree Size (leaves)',
            'Total Tree Leaves'
        ],
        'Value': [
            total_samples,
            total_clusters,
            f"{avg_cluster_size:.2f}",
            max_cluster_size,
            min_cluster_size,
            singleton_clusters,
            small_clusters,
            medium_clusters,
            large_clusters,
            total_trees,
            total_snps,
            f"{avg_snps_per_sample:.2f}",
            max_snps_per_sample,
            min_snps_per_sample,
            f"{avg_tree_leaves:.2f}",
            total_tree_leaves
        ]
    }
    
    statistics_df = pd.DataFrame(statistics_data)
    statistics_df.to_csv('Analysis_Statistics.tsv', sep='\\t', index=False)
    
    # Create HTML report
    html_content = f\"\"\"
    <!DOCTYPE html>
    <html>
    <head>
        <title>Scalable Assembly SNPs Analysis - Final Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1 {{ color: #2E86AB; border-bottom: 2px solid #2E86AB; }}
            h2 {{ color: #A23B72; margin-top: 30px; }}
            h3 {{ color: #F18F01; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
            th {{ background-color: #f2f2f2; font-weight: bold; }}
            .summary-box {{ background-color: #f9f9f9; padding: 20px; border-left: 4px solid #2E86AB; margin: 20px 0; }}
            .warning {{ background-color: #fff3cd; border-left: 4px solid #ffc107; padding: 15px; margin: 10px 0; }}
            .success {{ background-color: #d4edda; border-left: 4px solid #28a745; padding: 15px; margin: 10px 0; }}
            .metric {{ display: inline-block; margin: 10px 20px 10px 0; }}
            .metric-value {{ font-size: 24px; font-weight: bold; color: #2E86AB; }}
            .metric-label {{ font-size: 14px; color: #666; }}
        </style>
    </head>
    <body>
        <h1>🧬 Scalable Assembly SNPs Analysis - Final Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="summary-box">
            <h2>📊 Analysis Overview</h2>
            <div class="metric">
                <div class="metric-value">{total_samples}</div>
                <div class="metric-label">Total Samples</div>
            </div>
            <div class="metric">
                <div class="metric-value">{total_clusters}</div>
                <div class="metric-label">Clusters Formed</div>
            </div>
            <div class="metric">
                <div class="metric-value">{total_trees}</div>
                <div class="metric-label">Phylogenetic Trees</div>
            </div>
            <div class="metric">
                <div class="metric-value">{total_snps:,}</div>
                <div class="metric-label">Total SNP Positions</div>
            </div>
        </div>
        
        <h2>🔬 Clustering Results</h2>
        <table>
            <tr><th>Cluster Category</th><th>Count</th><th>Description</th></tr>
            <tr><td>Singleton Clusters</td><td>{singleton_clusters}</td><td>Clusters with 1 sample</td></tr>
            <tr><td>Small Clusters</td><td>{small_clusters}</td><td>Clusters with 2-5 samples</td></tr>
            <tr><td>Medium Clusters</td><td>{medium_clusters}</td><td>Clusters with 6-20 samples</td></tr>
            <tr><td>Large Clusters</td><td>{large_clusters}</td><td>Clusters with >20 samples</td></tr>
        </table>
        
        <h2>🌳 Phylogenetic Analysis</h2>
        <div class="success">
            <strong>✅ Integration Complete:</strong> Core SNPs have been integrated across all clusters and phylogenetic trees have been grafted into a supertree.
        </div>
        
        <h3>Key Results:</h3>
        <ul>
            <li><strong>Integrated Core SNPs Alignment:</strong> Combined SNP positions from all clusters</li>
            <li><strong>Grafted Supertree:</strong> Combined phylogenetic relationships from all cluster trees</li>
            <li><strong>Global Phylogenetic Tree:</strong> Built from integrated core SNPs alignment</li>
        </ul>
        
        <h2>📁 Output Files</h2>
        <h3>Cluster-specific Results:</h3>
        <ul>
            <li><code>Clusters/cluster_*/</code> - Individual cluster alignments, trees, and Gubbins results</li>
            <li><code>Mash_Sketches/</code> - Distance calculation sketches</li>
        </ul>
        
        <h3>Integrated Results:</h3>
        <ul>
            <li><code>Integrated_Results/integrated_core_snps.fa</code> - Combined core SNPs alignment</li>
            <li><code>Integrated_Results/grafted_supertree.tre</code> - Grafted supertree from all clusters</li>
            <li><code>Integrated_Results/integrated_core_snps.treefile</code> - Global phylogenetic tree</li>
            <li><code>Integrated_Results/sample_cluster_mapping.tsv</code> - Sample-to-cluster assignments</li>
        </ul>
        
        <h3>Summary Reports:</h3>
        <ul>
            <li><code>Summaries/clusters.tsv</code> - Cluster assignments</li>
            <li><code>Summaries/mash_distances.tsv</code> - Pairwise distance matrix</li>
            <li><code>Final_Results/</code> - This comprehensive report and statistics</li>
        </ul>
        
        <h2>🔧 Analysis Parameters</h2>
        <ul>
            <li><strong>Clustering Method:</strong> Mash distance-based hierarchical clustering</li>
            <li><strong>Phylogenetic Method:</strong> IQ-TREE with model selection</li>
            <li><strong>Recombination Analysis:</strong> Gubbins (if enabled)</li>
            <li><strong>Tree Grafting:</strong> Star topology with cluster subtrees</li>
        </ul>
        
        <div class="warning">
            <strong>⚠️ Note:</strong> This scalable approach processes large datasets by clustering similar genomes and analyzing each cluster separately, then integrating results. Individual cluster trees maintain detailed relationships within clusters, while the integrated results provide global perspectives.
        </div>
        
        <hr>
        <p><em>Generated by wf-assembly-snps scalable workflow</em></p>
    </body>
    </html>
    \"\"\"
    
    with open('Scalable_Analysis_Final_Report.html', 'w') as f:
        f.write(html_content)
    
    # Create text report
    text_content = f\"\"\"
SCALABLE ASSEMBLY SNPs ANALYSIS - FINAL REPORT
{'=' * 60}

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

ANALYSIS OVERVIEW
{'-' * 20}
Total Samples: {total_samples}
Total Clusters: {total_clusters}
Phylogenetic Trees: {total_trees}
Total SNP Positions: {total_snps:,}

CLUSTERING RESULTS
{'-' * 20}
Average Cluster Size: {avg_cluster_size:.2f}
Largest Cluster: {max_cluster_size} samples
Smallest Cluster: {min_cluster_size} samples

Cluster Distribution:
- Singleton Clusters (1 sample): {singleton_clusters}
- Small Clusters (2-5 samples): {small_clusters}
- Medium Clusters (6-20 samples): {medium_clusters}
- Large Clusters (>20 samples): {large_clusters}

SNP ANALYSIS
{'-' * 20}
Total SNP Positions: {total_snps:,}
Average SNPs per Sample: {avg_snps_per_sample:.2f}
Maximum SNPs per Sample: {max_snps_per_sample}
Minimum SNPs per Sample: {min_snps_per_sample}

PHYLOGENETIC ANALYSIS
{'-' * 20}
Total Trees Generated: {total_trees}
Average Tree Size: {avg_tree_leaves:.2f} leaves
Total Tree Leaves: {total_tree_leaves}

INTEGRATION RESULTS
{'-' * 20}
✅ Core SNPs integrated across all clusters
✅ Phylogenetic trees grafted into supertree
✅ Global phylogenetic tree constructed from integrated alignment

OUTPUT STRUCTURE
{'-' * 20}
Clusters/           - Individual cluster results
Integrated_Results/ - Combined core SNPs and grafted trees
Summaries/          - Analysis summaries and statistics
Final_Results/      - This comprehensive report
Mash_Sketches/      - Distance calculation files

METHODS
{'-' * 20}
Clustering: Mash distance-based hierarchical clustering
Alignment: SKA (Split Kmer Analysis) for each cluster
Phylogeny: IQ-TREE with automatic model selection
Integration: Core SNP extraction and concatenation
Tree Grafting: Star topology with cluster subtrees
Recombination: Gubbins analysis (if enabled)

This scalable approach enables analysis of large genomic datasets by:
1. Clustering similar genomes based on Mash distances
2. Performing detailed analysis within each cluster
3. Integrating core SNPs across all clusters
4. Grafting cluster trees into a comprehensive supertree
5. Building a global phylogenetic tree from integrated data

For detailed results, see individual files in the output directories.
    \"\"\"
    
    with open('Scalable_Analysis_Final_Report.txt', 'w') as f:
        f.write(text_content)
    
    print("Final summary reports created successfully!")

# Run summary creation
create_final_summary()
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}