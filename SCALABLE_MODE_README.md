# Scalable Mode for wf-assembly-snps

## Overview

The scalable mode has been implemented to handle large datasets by:

1. **Clustering genomes** using Mash distances
2. **Building phylogenetic trees** for each cluster using SKA and IQ-TREE
3. **Running Gubbins** on each cluster for recombination analysis
4. **Combining results** into a comprehensive summary

## Key Features Fixed

### 1. Clustering Pipeline
- ✅ Mash sketching of all input genomes
- ✅ Pairwise distance calculation
- ✅ Hierarchical clustering with configurable thresholds
- ✅ Cluster size management and singleton handling

### 2. Per-Cluster Analysis
- ✅ SKA-based alignment generation for each cluster
- ✅ Fast ML tree building with IQ-TREE
- ✅ Gubbins recombination analysis (when enabled)
- ✅ Proper result publishing to organized directories

### 3. Result Organization
- ✅ `Clusters/`: Individual cluster results (alignments, trees, Gubbins output)
- ✅ `Summaries/`: Clustering summaries and analysis reports
- ✅ `Mash_Sketches/`: Mash sketch files for distance calculations

### 4. Configuration Parameters

```bash
# Enable scalable mode
--scalable_mode

# Clustering parameters
--mash_threshold 0.03          # Distance threshold for clustering (lower = more clusters)
--max_cluster_size 100         # Maximum samples per cluster
--merge_singletons             # Merge singleton clusters for analysis

# Gubbins parameters
--run_gubbins                  # Enable recombination analysis
--gubbins_iterations 3         # Number of Gubbins iterations
--gubbins_tree_builder iqtree  # Tree builder for Gubbins
--gubbins_use_hybrid           # Use hybrid tree building approach
```

## Usage Example

```bash
nextflow run main.nf \\
  -profile docker \\
  --input /path/to/assemblies \\
  --outdir results \\
  --scalable_mode \\
  --run_gubbins \\
  --mash_threshold 0.05 \\
  --merge_singletons
```

## Output Structure

```
results/
├── Clusters/
│   ├── cluster_1/
│   │   ├── cluster_1.aln.fa                    # SKA alignment
│   │   ├── cluster_1.treefile                  # IQ-TREE phylogeny
│   │   ├── cluster_1.filtered_polymorphic_sites.fasta  # Gubbins filtered alignment
│   │   ├── cluster_1.recombination_predictions.gff     # Gubbins recombination
│   │   └── cluster_1.node_labelled.final_tree.tre      # Gubbins final tree
│   └── cluster_2/
│       └── ...
├── Summaries/
│   ├── clusters.tsv                           # Cluster assignments
│   ├── cluster_summary.txt                    # Clustering statistics
│   ├── mash_distances.tsv                     # Pairwise distances
│   ├── Scalable_Analysis_Summary.tsv          # Overall summary
│   └── Scalable_Analysis_Report.txt           # Detailed report
└── Mash_Sketches/
    ├── sample1.msh
    ├── sample2.msh
    └── ...
```

## Troubleshooting

### No clusters formed
- Increase `--mash_threshold` (e.g., 0.1) to allow more diverse samples to cluster
- Use `--merge_singletons` to combine all singletons into one large cluster

### Small clusters skipped
- Phylogenetic analysis requires ≥3 samples per cluster
- Use `--merge_singletons` to combine small clusters

### Gubbins failures
- Some clusters may fail Gubbins analysis due to insufficient variation
- Empty output files are created to maintain pipeline flow
- Check individual cluster logs for details

## Performance Notes

- Clustering scales well with large datasets
- Per-cluster analysis is parallelized
- Memory usage is distributed across clusters
- Consider adjusting `--max_cluster_size` for very large clusters