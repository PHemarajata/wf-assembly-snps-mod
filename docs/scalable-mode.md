# Scalable Mode for Large Datasets

This document describes the scalable mode implementation that allows the workflow to handle hundreds to thousands of bacterial genome assemblies efficiently.

## Overview

The original workflow uses Parsnp for global core genome alignment followed by Gubbins for recombination detection. While effective for smaller datasets, this approach doesn't scale well beyond ~200 genomes due to:

1. **Memory bottlenecks**: Global multiple sequence alignment requires substantial RAM
2. **Computational complexity**: Gubbins iterative ML tree building scales poorly
3. **No incremental capability**: Adding new samples requires complete re-analysis

The scalable mode implements a **divide-and-conquer approach** with the following components:

## Architecture

### 1. Pre-clustering (Divide)
- **Mash sketching**: Fast k-mer based distance estimation
- **Distance-based clustering**: Single-linkage clustering to group similar genomes
- **Cluster size control**: Automatic splitting of large clusters (default: max 100 genomes per cluster)

### 2. Per-cluster Analysis (Conquer)
- **SKA alignment**: Fast split-k-mer based SNP alignment (reference-free)
- **IQ-TREE2**: Rapid maximum likelihood tree construction
- **Optional Gubbins**: Recombination detection only within clusters (with optimized settings)

### 3. Global Integration
- **Backbone tree**: Global phylogeny from cluster representatives or concatenated SNPs
- **UShER integration**: Mutation-annotated trees for ultra-fast phylogenetic placement

### 4. Incremental Updates
- **UShER placement**: Add new genomes in seconds without re-running entire analysis
- **Cluster-based updates**: Only re-analyze affected clusters when adding samples

## Usage

### Basic Scalable Mode

```bash
nextflow run main.nf \
  --input assemblies/ \
  --outdir results \
  --scalable_mode true \
  --profile local_workstation
```

### Advanced Configuration

```bash
nextflow run main.nf \
  --input assemblies/ \
  --outdir results \
  --scalable_mode true \
  --workflow_mode cluster \
  --mash_threshold 0.025 \
  --max_cluster_size 80 \
  --run_gubbins true \
  --gubbins_iterations 3 \
  --build_usher_mat true \
  --profile dgx_station
```

### Incremental Updates (Future Implementation)

```bash
# Initial analysis
nextflow run main.nf \
  --input initial_assemblies/ \
  --scalable_mode true \
  --build_usher_mat true \
  --outdir results

# Add new samples
nextflow run main.nf \
  --input new_assemblies/ \
  --workflow_mode place \
  --existing_mat results/global.pb \
  --outdir results_updated
```

## Parameters

### Clustering Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `scalable_mode` | `false` | Enable scalable clustering-based workflow |
| `workflow_mode` | `"cluster"` | Workflow mode: `cluster`, `place`, or `global` |
| `mash_threshold` | `0.03` | Distance threshold for clustering (lower = more clusters) |
| `max_cluster_size` | `100` | Maximum genomes per cluster |

### Gubbins Optimization

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_gubbins` | `true` | Run Gubbins within clusters |
| `gubbins_iterations` | `3` | Maximum Gubbins iterations (reduced from default 5) |
| `gubbins_use_hybrid` | `true` | Enable hybrid tree building approach |
| `gubbins_first_tree_builder` | `"rapidnj"` | Fast tree builder for initial tree (rapidnj, fasttree) |
| `gubbins_tree_builder` | `"iqtree"` | Accurate tree builder for refinement (iqtree, raxml) |
| `gubbins_min_snps` | `5` | Minimum SNPs to call recombination |

### Tree Building

| Parameter | Default | Description |
|-----------|---------|-------------|
| `iqtree_model` | `"GTR+ASC"` | IQ-TREE substitution model |

### UShER Integration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `build_usher_mat` | `false` | Build mutation-annotated tree for incremental updates |
| `existing_mat` | `null` | Path to existing MAT file for placement mode |

## Resource Profiles

### Local Workstation (12 cores, 64GB RAM)

```bash
--profile local_workstation
```

Optimized for:
- Limited concurrent Gubbins processes (maxForks = 1)
- Conservative memory allocation
- Suitable for 200-500 genomes

### DGX Station A100 (128 cores, 512GB RAM)

```bash
--profile dgx_station
```

Optimized for:
- High parallelization (multiple concurrent Gubbins processes)
- Large memory allocation
- Suitable for 1000+ genomes

## Performance Expectations

### Scalability Comparison

| Dataset Size | Original Workflow | Scalable Mode |
|--------------|-------------------|---------------|
| 50 genomes | ~2 hours | ~1 hour |
| 200 genomes | ~8 hours | ~2 hours |
| 500 genomes | >24 hours* | ~4 hours |
| 1000 genomes | Not feasible* | ~8 hours |

*Estimates based on memory limitations and computational complexity

### Resource Usage

**Local Workstation (12 cores, 64GB RAM):**
- Recommended: Up to 500 genomes
- Peak memory: ~32GB (vs >64GB for original)
- Runtime: ~4-6 hours for 500 genomes

**DGX Station (128 cores, 512GB RAM):**
- Recommended: Up to 2000+ genomes
- Peak memory: ~256GB
- Runtime: ~6-12 hours for 1000+ genomes

## Algorithm Details

### Clustering Strategy

1. **Mash sketching**: Create k-mer sketches for each assembly
2. **Distance calculation**: Compute pairwise Mash distances
3. **Graph clustering**: Build distance graph and find connected components
4. **Size control**: Split clusters exceeding `max_cluster_size`

### Per-cluster Processing

1. **SKA build**: Create split-k-mer alignment database
2. **SKA align**: Generate SNP alignment (much faster than Parsnp)
3. **IQ-TREE2**: Build ML tree with `--fast` mode
4. **Gubbins** (optional): Detect recombination with optimized settings

### Optimization Strategies

**For B. pseudomallei (high recombination):**
- Use cluster-local Gubbins masks (avoid over-aggressive global masking)
- Reduce Gubbins iterations (3 vs default 5)
- Use hybrid tree building approach (RapidNJ → IQ-TREE refinement)
- Set minimum SNP threshold to avoid noise

## Troubleshooting

### Memory Issues
- Reduce `max_cluster_size` (try 50-80)
- Increase `mash_threshold` (try 0.04-0.05) to create fewer, larger clusters
- Use `local_workstation` profile for conservative memory allocation

### Too Many Small Clusters
- Decrease `mash_threshold` (try 0.02-0.025)
- Check input data quality (contamination, misassembly)

### Slow Performance
- Disable Gubbins for initial analysis (`--run_gubbins false`)
- Use DGX profile for high-performance systems
- Consider increasing cluster sizes for fewer parallel processes

## Future Enhancements

1. **Complete UShER integration**: Full implementation of incremental placement mode
2. **Adaptive clustering**: Dynamic threshold selection based on dataset characteristics  
3. **Hybrid recombination masking**: Combine cluster-local and global masks
4. **Cloud deployment**: Optimized profiles for AWS, GCP, Azure
5. **Real-time monitoring**: Progress tracking and resource usage visualization

## References

- **Mash**: Ondov et al. (2016) "Mash: fast genome and metagenome distance estimation using MinHash"
- **SKA**: Harris et al. (2018) "SKA: Split Kmer Analysis Toolkit for the rapid and standardized investigation of bacterial genomic epidemiology"
- **IQ-TREE2**: Minh et al. (2020) "IQ-TREE 2: New models and efficient methods for phylogenetic inference"
- **UShER**: Turakhia et al. (2021) "Ultrafast Sample placement on Existing tRees (UShER)"
- **Gubbins**: Croucher et al. (2015) "Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences"