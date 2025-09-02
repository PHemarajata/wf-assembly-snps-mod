# Recombination-Aware Assembly SNPs Workflow

## Overview

This workflow implements a recombination-aware approach for phylogenetic analysis of bacterial genomes, specifically designed for *B. pseudomallei* and similar organisms. The workflow performs per-cluster recombination detection followed by tree grafting to create a comprehensive phylogenetic tree.

## Workflow Steps

### 1. Cluster with Mash
- **Process**: `MASH_SKETCH` → `MASH_DIST` → `CLUSTER_GENOMES`
- **Purpose**: Group similar assemblies using Mash distance-based clustering
- **Output**: `cluster_id` → list of assemblies

### 2. Per-cluster Whole/Core Alignment
- **Process**: `SELECT_CLUSTER_REPRESENTATIVE` → `SNIPPY_ALIGN` (preferred) or `PARSNP_ALIGN` → `KEEP_INVARIANT_ATCG`
- **Purpose**: Create whole genome alignments preserving invariant A/T/C/G sites
- **Key Features**:
  - Uses medoid/high-quality assembly as per-cluster reference
  - Keeps invariant sites (not SNP-only) for Gubbins
  - Snippy preferred for accuracy, Parsnp as fallback

### 3. Gubbins on the WGA
- **Process**: `IQTREE_FAST` (starting tree) → `GUBBINS_CLUSTER`
- **Purpose**: Detect and mask recombinant regions
- **Output per cluster**:
  - `cluster_X.recombination_predictions.gff`
  - `cluster_X.node_labelled.final_tree.tre`
  - `cluster_X.filtered_polymorphic_sites.fasta` (recombination-masked SNPs)

### 4. Per-cluster Final ML Tree
- **Process**: `IQTREE_ASC`
- **Purpose**: Build well-resolved local topology using ASC correction
- **Features**:
  - Uses GTR+ASC model with Lewis correction for SNP-only data
  - Bootstrap support (1000 replicates)
  - SH-aLRT support (1000 replicates)

### 5. Select Representatives per Cluster
- **Process**: `COLLECT_REPRESENTATIVES`
- **Purpose**: Gather medoid/high-quality assemblies used as references
- **Output**: Representatives FASTA and mapping file

### 6. Backbone Tree on Representatives
- **Process**: `BUILD_BACKBONE_TREE`
- **Purpose**: Create backbone phylogeny from cluster representatives
- **Methods**: Parsnp core alignment (preferred) or distance-based FastTree

### 7. Graft Subtrees onto Backbone
- **Process**: `GRAFT_SUBTREES`
- **Purpose**: Combine cluster trees into comprehensive phylogeny
- **Method**: `gotree graft` - attach cluster subtrees at representative positions
- **Output**: `global_grafted.treefile`

## Usage

### Basic Command
```bash
nextflow run main.nf \
  -profile docker \
  --recombination_aware_mode \
  --input assemblies/ \
  --outdir results_recombination_aware
```

### With Custom Parameters
```bash
nextflow run main.nf \
  -profile docker \
  --recombination_aware_mode \
  --input assemblies/ \
  --outdir results_recombination_aware \
  --alignment_method snippy \
  --mash_threshold 0.03 \
  --gubbins_iterations 5 \
  --backbone_method parsnp
```

### Test Configuration
```bash
nextflow run main.nf \
  -c test_recombination_aware.config \
  --input test_input/
```

## Parameters

### Core Parameters
- `--recombination_aware_mode`: Enable the recombination-aware workflow
- `--input`: Path to directory containing assembly files
- `--outdir`: Output directory for results

### Alignment Parameters
- `--alignment_method`: Method for per-cluster alignment (`snippy` or `parsnp`)
- `--mash_threshold`: Distance threshold for clustering (default: 0.03)
- `--max_cluster_size`: Maximum samples per cluster (default: 100)

### Gubbins Parameters
- `--gubbins_iterations`: Number of Gubbins iterations (default: 3)
- `--gubbins_min_snps`: Minimum SNPs for recombination detection (default: 3)
- `--gubbins_tree_builder`: Tree builder for Gubbins (`iqtree` or `fasttree`)

### IQ-TREE ASC Parameters
- `--iqtree_asc_model`: Model for ASC correction (default: `GTR+ASC`)
- `--iqtree_asc_correction`: ASC correction method (default: `lewis`)

### Tree Construction Parameters
- `--backbone_method`: Method for backbone tree (`parsnp` or `fasttree`)
- `--rescale_grafted_branches`: Rescale branch lengths after grafting

## Output Structure

```
results_recombination_aware/
├── Clusters/                          # Per-cluster results
│   ├── cluster_1/
│   │   ├── cluster_1.core.full.aln   # Whole genome alignment
│   │   ├── cluster_1.recombination_predictions.gff
│   │   ├── cluster_1.filtered_polymorphic_sites.fasta
│   │   └── cluster_1.final.treefile  # Per-cluster ML tree
│   └── cluster_N/
├── Representatives/
│   ├── representatives.fa            # All cluster representatives
│   └── cluster_representatives.tsv   # Representative mapping
├── Backbone/
│   ├── backbone.treefile             # Backbone tree from representatives
│   └── backbone_alignment.fa         # Backbone alignment
├── Final_Results/
│   ├── global_grafted.treefile       # Final grafted tree
│   ├── grafting_report.txt           # Grafting summary
│   └── tree_grafting_report.txt      # Detailed grafting log
└── Summaries/
    ├── clusters.tsv                  # Cluster assignments
    ├── mash_distances.tsv            # Distance matrix
    └── Summary.QC_File_Checks.tsv    # QC summary
```

## Key Features

### Guardrails (as specified)
- ✅ **Do not feed SNP-only alignments to Gubbins** - Uses whole genome alignments with invariant sites
- ✅ **Keep invariant A/T/C/G columns** - Preserved in Step 2 alignment
- ✅ **Small cluster handling** - Uses `--min-snps 1-3` for small clusters
- ✅ **Consistent leaf names** - Representatives maintained between cluster and backbone trees

### Scalability
- **Clustering approach**: Handles large datasets by processing similar genomes together
- **Per-cluster analysis**: Detailed recombination detection within clusters
- **Tree grafting**: Combines local topologies into global phylogeny
- **Resource efficiency**: Parallel processing of independent clusters

### Recombination Awareness
- **Gubbins integration**: Detects and masks recombinant regions per cluster
- **Whole genome context**: Preserves genomic context for accurate recombination detection
- **ASC correction**: Proper handling of SNP-only data in final trees
- **Local resolution**: High-quality trees within clusters before grafting

## Comparison with Other Modes

| Feature | Standard Mode | Scalable Mode | Recombination-Aware Mode |
|---------|---------------|---------------|--------------------------|
| Recombination Detection | ❌ | ❌ | ✅ (Gubbins per cluster) |
| Tree Grafting | ❌ | ✅ (Simple) | ✅ (Sophisticated) |
| Whole Genome Alignment | ❌ | ❌ | ✅ (Per cluster) |
| ASC Correction | ❌ | ❌ | ✅ (For SNP trees) |
| Scalability | Limited | High | High |
| Accuracy | Medium | Medium | High |

## Requirements

### Software Dependencies
- Nextflow ≥22.04.3
- Docker or Singularity
- Mash (for clustering)
- Snippy (preferred alignment) or Parsnp (fallback)
- Gubbins (recombination detection)
- IQ-TREE2 (phylogenetic inference)
- gotree (tree grafting)

### Resource Recommendations
- **Minimum**: 8 CPUs, 16GB RAM
- **Recommended**: 16+ CPUs, 32+ GB RAM
- **Storage**: ~10GB per 100 assemblies

## Troubleshooting

### Common Issues
1. **Small clusters**: Adjust `--mash_threshold` or use `--merge_singletons`
2. **Snippy failures**: Fallback to `--alignment_method parsnp`
3. **Gubbins errors**: Check `--gubbins_min_snps` for small clusters
4. **Grafting failures**: Verify representative names match between trees

### Performance Optimization
- Use `--max_cluster_size` to limit cluster sizes
- Adjust `--gubbins_iterations` based on dataset complexity
- Use appropriate resource profiles for your system

## Citation

If you use this recombination-aware workflow, please cite:
- The original wf-assembly-snps pipeline
- Gubbins for recombination detection
- IQ-TREE for phylogenetic inference
- gotree for tree manipulation
- Snippy for genome alignment (if used)