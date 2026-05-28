# wf-assembly-snps-mod

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.3-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-blue.svg)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/singularity-blue.svg)](https://sylabs.io/docs/)
[![Conda](https://img.shields.io/badge/conda-green.svg)](https://docs.conda.io/en/latest/)

A modified and enhanced Nextflow workflow for bacterial genome assembly-based SNP identification and phylogenetic analysis.

> **Note**: This is a modified version of [bacterial-genomics/wf-assembly-snps](https://github.com/bacterial-genomics/wf-assembly-snps) with additional features and optimizations.

## 🚀 Quick Start

```bash
# Typical usage for Burkholderia pseudomallei SNP analysis on APHL Analysis Laptop
nextflow run main.nf -profile local_workstation_rtx4070,docker \
  --input /path/to/assemblies  \
  --outdir /path/to/results \
  --recombination_aware_mode true \
  --integrate_results true \
  --mash_threshold 0.028  \
  --max_cluster_size 50  \   --merge_singletons true \
  --mash_sketch_size 50000  \
  --recombination gubbins \
  --snp_package parsnp \
  --run_gubbins

# Basic usage with input directory
nextflow run PHemarajata/wf-assembly-snps-mod \
  -profile docker \
  --input /path/to/assemblies \
  --outdir results

# Scalable mode for large datasets (>200 genomes)
nextflow run PHemarajata/wf-assembly-snps-mod \
  -profile docker \
  --input /path/to/assemblies \
  --outdir results \
  --scalable_mode true

# Recombination-aware analysis
nextflow run PHemarajata/wf-assembly-snps-mod \
  -profile docker \
  --input /path/to/assemblies \
  --outdir results \
  --recombination_aware_mode true
```

## 📋 Overview

This workflow identifies single nucleotide polymorphisms (SNPs) from bacterial genome assemblies and constructs phylogenetic trees. It offers three distinct analysis modes:

### 🔬 **Standard Mode**
- Core genome alignment using **Parsnp**
- SNP distance matrix calculation
- Maximum likelihood phylogeny construction
- Suitable for datasets up to ~200 genomes

### 📈 **Scalable Mode** 
- **Divide-and-conquer** approach for large datasets (hundreds to thousands of genomes)
- Pre-clustering using **Mash** k-mer distances
- Per-cluster analysis with **SKA** alignment and **IQ-TREE2**
- Global integration with backbone phylogeny
- Support for **UShER** incremental updates

### 🧬 **Recombination-Aware Mode**
- Recombination detection using **Gubbins**
- Masked SNP analysis excluding recombinant regions
- Enhanced phylogenetic accuracy for highly recombinogenic species

## 🎯 Key Features

- **Multiple analysis modes** optimized for different dataset sizes and biological questions
- **Flexible input formats**: Directory of FASTA files or CSV samplesheet
- **Comprehensive output**: SNP matrices, phylogenetic trees, quality reports
- **HPC ready** with built-in profiles for various compute environments
- **Containerized** with Docker, Singularity, and Conda support
- **Reproducible** with detailed provenance tracking

## 📊 Input Requirements

### Accepted File Formats
- **FASTA files** with extensions: `.fasta`, `.fas`, `.fna`, `.fsa`, `.fa`
- **Optional compression**: gzip (`.gz`)
- **File naming**: No spaces, unique filenames required
- **Minimum size**: 45kb (configurable with `--min_input_filesize`)

### Input Methods

#### Option 1: Directory Input
```bash
--input /path/to/assemblies/
```

#### Option 2: Samplesheet Input
```csv
sample,file
SAMPLE_1,/path/to/SAMPLE_1.fasta
SAMPLE_2,/path/to/SAMPLE_2.fasta
SAMPLE_3,/path/to/SAMPLE_3.fasta
```

## ⚙️ Configuration Modes

### Standard Mode (Default)
```bash
nextflow run PHemarajata/wf-assembly-snps-mod \
  --input assemblies/ \
  --outdir results
```

### Scalable Mode
For datasets with >200 genomes:
```bash
nextflow run PHemarajata/wf-assembly-snps-mod \
  --input assemblies/ \
  --outdir results \
  --scalable_mode true \
  --mash_threshold 0.025 \
  --max_cluster_size 100
```

### Recombination-Aware Mode
For accurate phylogeny in highly recombinogenic species:
```bash
nextflow run PHemarajata/wf-assembly-snps-mod \
  --input assemblies/ \
  --outdir results \
  --recombination_aware_mode true \
  --recombination gubbins
```

## 🛠️ Parameters

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `null` | Input directory or samplesheet (required) |
| `--outdir` | `null` | Output directory (required) |
| `--ref` | `null` | Reference genome (optional) |
| `--snp_package` | `parsnp` | SNP calling tool |
| `--min_input_filesize` | `45k` | Minimum input file size |

### Workflow Mode Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--scalable_mode` | `false` | Enable scalable clustering workflow |
| `--recombination_aware_mode` | `true`* | Enable recombination detection |
| `--workflow_mode` | `cluster` | Workflow mode: cluster/place/global |

### Clustering Parameters (Scalable Mode)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mash_threshold` | `0.028` | Distance threshold for clustering |
| `--max_cluster_size` | `50` | Maximum genomes per cluster |
| `--merge_singletons` | `true` | Merge singleton clusters |
| `--mash_sketch_size` | `50000` | Mash sketch size for large datasets |
| `--mash_kmer_size` | `21` | K-mer size for Mash |
| `--mash_min_copies` | `1` | Minimum k-mer copies |

### Recombination Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--recombination` | `gubbins` | Recombination detection tool |
| `--run_gubbins` | `true` | Enable Gubbins analysis |
| `--gubbins_iterations` | `3` | Maximum Gubbins iterations |
| `--gubbins_use_hybrid` | `true` | Use hybrid tree building |
| `--gubbins_first_tree_builder` | `rapidnj` | Fast initial tree builder |
| `--gubbins_tree_builder` | `iqtree` | Refined tree builder |
| `--gubbins_min_snps` | `2` | Minimum SNPs for analysis |

### Parsnp Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--curated_input` | `false` | Use curated input mode |
| `--tree_method` | `fasttree` | Tree method: fasttree/raxml |
| `--max_partition_size` | `15000` | Maximum partition size |

### IQ-TREE Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--iqtree_model` | `GTR+ASC` | Evolutionary model |
| `--iqtree_asc_model` | `GTR+ASC` | Ascertainment bias correction |

### Integration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--integrate_results` | `true` | Integrate cluster results |
| `--alignment_method` | `snippy` | Alignment method (snippy/parsnp) |
| `--backbone_method` | `parsnp` | Backbone tree method |

### UShER Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--build_usher_mat` | `false` | Build UShER mutation tree |
| `--existing_mat` | `null` | Existing UShER tree file |

### Output Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--create_excel_outputs` | `false` | Create Excel format outputs |
| `--excel_sheet_name` | `Sheet1` | Excel sheet name |
| `--publish_dir_mode` | `copy` | How to publish outputs |

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_memory` | `128.GB` | Maximum memory per process |
| `--max_cpus` | `16` | Maximum CPUs per process |
| `--max_time` | `240.h` | Maximum runtime per process |

> **Note**: Default values marked with * may be overridden in specific configuration profiles

## 🖥️ Compute Profiles

Pre-configured profiles for different computing environments:

- **`docker`** - Docker containers (default)
- **`singularity`** - Singularity containers
- **`conda`** - Conda environments
- **`local_workstation`** - Local workstation (12 cores, 64GB RAM)
- **`dgx_station`** - DGX Station A100 (128 cores, 512GB RAM)
- **`aspen_hpc`** - Aspen HPC cluster
- **`rosalind_hpc`** - Rosalind HPC cluster

## 📈 Scalable Mode Details

The scalable mode implements a sophisticated divide-and-conquer approach:

### 1. **Pre-clustering Phase**
- Fast k-mer distance estimation with **Mash**
- Single-linkage clustering to group similar genomes
- Automatic cluster size optimization

### 2. **Per-cluster Analysis**
- Reference-free SNP alignment with **SKA**
- Maximum likelihood phylogeny with **IQ-TREE2**
- Optional recombination detection with **Gubbins**

### 3. **Global Integration**
- Backbone phylogeny construction
- Cross-cluster SNP distance matrices
- Optional **UShER** mutation-annotated trees

## 📁 Output Structure

```
results/
├── alignments/          # Core genome alignments
├── snp_distances/       # SNP distance matrices
├── phylogeny/          # Phylogenetic trees (Newick format)
├── reports/            # Quality control reports
├── gubbins/            # Recombination analysis (if enabled)
├── clusters/           # Per-cluster results (scalable mode)
│   ├── cluster_1/      # Individual cluster results
│   └── cluster_N/      
├── backbone.treefile   # Backbone tree (scalable mode)
├── cluster_representatives.tsv  # Cluster representative mappings
├── final_grafted.treefile      # Complete grafted tree (if successful)
├── grafting_report.txt         # Tree grafting summary
├── grafting_log.txt           # Detailed grafting log
└── pipeline_info/      # Execution reports and logs
```

## 🔧 Installation

### Prerequisites
- **Nextflow** ≥22.04.3
- **Docker**, **Singularity**, or **Conda**

### Quick Installation
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Test the workflow
nextflow run PHemarajata/wf-assembly-snps-mod -profile test,docker
```

## 📖 Documentation

Comprehensive documentation is available in the [`docs/`](docs/) directory:

- **[Usage Guide](docs/usage.md)** - Detailed usage instructions and examples
- **[Output Description](docs/output.md)** - Complete output file descriptions
- **[Scalable Mode](docs/scalable-mode.md)** - In-depth scalable mode documentation
- **[HPC Configuration](docs/HPC-UGE-scheduler.md)** - HPC cluster setup

## 🚀 Convenience Wrapper

For easier usage, a wrapper script is provided:

```bash
# Make the script executable
chmod +x run_workflow.sh

# Standard analysis
./run_workflow.sh --input assemblies/ --outdir results

# Scalable mode for large datasets
./run_workflow.sh --input assemblies/ --mode scalable --profile local_workstation

# Recombination-aware analysis
./run_workflow.sh --input assemblies/ --mode recombination

# Pass additional parameters
./run_workflow.sh --input assemblies/ --mode scalable -- --mash_threshold 0.025
```

## 🧮 Downsampling Contextual Genomes

When you have far more contextual (public) assemblies than you can analyze — especially for a homogeneous species like *B. pseudomallei*, where dense, redundant sampling causes single‑linkage Mash clustering to chain unrelated clades into one mega‑cluster (which then gets sliced into arbitrary `max_cluster_size` chunks) — use the standalone `bin/downsample_contextual.py` script to build a balanced, de‑redundant input **before** running the pipeline.

It applies a **hybrid strategy**: (1) collapse near‑identical genomes *within each country* using Mash distance (one representative per redundancy group — this breaks the chaining and auto‑scales how many genomes are kept per region to the diversity actually present), then (2) enforce a per‑country floor and cap so rare regions stay represented and an over‑sampled region cannot dominate. Study isolates (e.g. CDC `IP-`/`IE-`) are **always kept in full**. The output is a `sample,file` samplesheet that `--input` accepts directly, so the expensive high‑resolution Mash run (`--mash_sketch_size 100000`) happens only on the reduced set.

### Prerequisites

```bash
# Only the `mash` binary (>=2.x) and Python 3 standard library are required.
mash --version
```

### Workflow

```bash
# Make the script executable
chmod +x bin/downsample_contextual.py

# Step 1 (optional): sweep dereplication thresholds to see how aggressively the
# population collapses, then pick --derep-threshold from the report.
./bin/downsample_contextual.py \
  --contextual-dir contextual_fasta \
  --cdc-dir cdc_fasta \
  --metadata contextual_fasta/megamix_bestshot_cleaned_dropGCF_on_Fdups.tsv \
  --derep-sketch-size 50000 --threads 16 \
  --sweep 0.0001,0.0005,0.001,0.005,0.01 \
  --outdir downsample_out
#   -> downsample_out/derep_sweep.tsv  (and a cached downsample_out/mash_distances.tsv)

# Step 2: select genomes. Reuse the cached distances from Step 1 with --mash-dist.
./bin/downsample_contextual.py \
  --contextual-dir contextual_fasta \
  --cdc-dir cdc_fasta \
  --metadata contextual_fasta/megamix_bestshot_cleaned_dropGCF_on_Fdups.tsv \
  --mash-dist downsample_out/mash_distances.tsv \
  --derep-threshold 0.001 \
  --target-total 1000 \
  --min-per-country 3 --max-per-country 200 \
  --seed 42 \
  --outdir downsample_out

# Step 3: feed the samplesheet to the pipeline at high Mash resolution.
nextflow run main.nf \
  -profile bp,docker \
  --input downsample_out/samplesheet.csv \
  --recombination_aware_mode --integrate_results \
  --mash_sketch_size 100000 \
  --mash_threshold 0.001 --max_cluster_size 50 \
  -resume
```

### Key Options

| Option | Default | Description |
|--------|---------|-------------|
| `--target-total` | 1000 | Target total genomes (study isolates + downsampled contextual). |
| `--derep-threshold` | 0.0005 | Mash distance at/below which genomes are treated as redundant. |
| `--derep-sketch-size` | 50000 | Sketch size for the one‑time selection Mash run. |
| `--min-per-country` / `--max-per-country` | 3 / 200 | Per‑country floor and cap for geographic balancing. |
| `--mash-dist` | — | Reuse a precomputed `mash triangle`/`mash dist` file (skips running Mash). |
| `--sweep` | — | Comma‑separated thresholds; write a group‑count report and exit. |
| `--unmatched` | drop | `drop` excludes genomes with no metadata row; `unknown` keeps them in an `Unknown` bucket. |
| `--seed` | 42 | Seed for reproducible representative selection. |

### Outputs (in `--outdir`)

- **`samplesheet.csv`** — `sample,file` for the pipeline (study isolates + selected contextual).
- **`selection_report.tsv`** — per‑genome fate: country, subregion, date, group, role (`cdc`/`representative`/`redundant`), kept, reason.
- **`country_summary.tsv`** — per‑country `n_input`, `n_after_derep`, `n_final`.
- **`unmatched_no_metadata.tsv`** — genomes with no metadata row and how they were dispositioned.
- **`mash_distances.tsv`** — cached Mash distances (reuse via `--mash-dist`).

### Metadata matching notes

Genomes are matched to the metadata TSV by **accession with the version ignored** (e.g. on‑disk `GCA_000182195.1.fasta` matches a `GCA_000182195_2` metadata row), with `FASTA_name`/`sample_id` as fallbacks; this is robust to the dot‑vs‑underscore version style differences between filenames and the TSV. Duplicate filename copies of the same genome (e.g. `GCA_x.2.fasta` and `GCA_x_2.fasta`) are collapsed automatically.

## 🌳 Standalone Tree Grafting

When the main workflow encounters issues in the final tree grafting step, you can use the standalone `graft_trees.py` script to complete the phylogenetic analysis separately.

### Background
In scalable mode, the workflow generates:
1. **Backbone tree** - Global phylogeny from cluster representatives
2. **Cluster trees** - Individual phylogenies for each genome cluster
3. **Final step** - Grafting cluster subtrees onto the backbone tree

The tree grafting step can sometimes fail due to memory constraints, label conflicts, or tree structure issues. The standalone script provides a robust solution with detailed logging and error handling.

### Prerequisites

```bash
# Install required Python package
pip install biopython
```

### Basic Usage

```bash
# Make the script executable
chmod +x graft_trees.py

# Basic tree grafting
./graft_trees.py \
  --backbone results/backbone.treefile \
  --clusters 'results/Clusters/**/cluster_*.final.treefile' \
  --reps results/cluster_representatives.tsv \
  --out-tree results/final_grafted.treefile \
  --report results/grafting_report.txt \
  --log results/grafting_log.txt
```

### Advanced Options

```bash
# With conflict resolution and detailed logging
./graft_trees.py \
  --backbone results/backbone.treefile \
  --clusters 'results/Clusters/**/cluster_*.final.treefile' \
  --reps results/cluster_representatives.tsv \
  --out-tree results/final_grafted.treefile \
  --report results/grafting_report.txt \
  --log results/grafting_log.txt \
  --rename-conflicts \
  --parent-edge-mode keep

# Dry run to check what would be done
./graft_trees.py \
  --backbone results/backbone.treefile \
  --clusters 'results/Clusters/**/cluster_*.final.treefile' \
  --reps results/cluster_representatives.tsv \
  --out-tree results/final_grafted.treefile \
  --report results/grafting_report.txt \
  --log results/grafting_log.txt \
  --dry-run
```

### Parameters

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--backbone` | ✅ | Backbone Newick tree file |
| `--clusters` | ✅ | Glob pattern for cluster tree files (repeatable) |
| `--reps` | ❌ | TSV file with cluster→representative mappings |
| `--out-tree` | ✅ | Output combined tree file |
| `--report` | ✅ | Summary report file |
| `--log` | ✅ | Detailed log file |
| `--rename-conflicts` | ❌ | Rename conflicting tip labels |
| `--parent-edge-mode` | ❌ | Branch length handling: `keep` or `zero` |
| `--dry-run` | ❌ | Plan only, don't write output tree |

### Expected Input Files

From a scalable workflow run, you'll typically find:

```
results/
├── backbone.treefile              # Global backbone phylogeny
├── cluster_representatives.tsv    # Cluster→representative mapping
├── Clusters/
│   ├── cluster_1/
│   │   └── cluster_1.final.treefile
│   ├── cluster_2/
│   │   └── cluster_2.final.treefile
│   └── ...
```

### Troubleshooting

**Common Issues:**

1. **Missing representative file**: If `cluster_representatives.tsv` is missing, the script will infer representatives automatically
2. **Label conflicts**: Use `--rename-conflicts` to automatically rename conflicting tip labels
3. **Memory issues**: The standalone script is more memory-efficient than the Nextflow process
4. **Tree structure problems**: Check the detailed log file for specific grafting failures

**Check your results:**

```bash
# Verify the final tree structure
python -c "
from Bio import Phylo
tree = Phylo.read('results/final_grafted.treefile', 'newick')
print(f'Final tree has {len(tree.get_terminals())} tips')
print(f'Tree depth: {tree.depth()}')
"

# View the grafting report
cat results/grafting_report.txt
```

### Example Usage Script

A complete example script is provided to demonstrate the typical tree grafting workflow:

```bash
# Run the example (after completing a scalable workflow)
chmod +x examples/run_tree_grafting_example.sh
./examples/run_tree_grafting_example.sh
```

## 🧪 Testing

```bash
# Quick test with sample data
nextflow run PHemarajata/wf-assembly-snps-mod -profile test,docker

# Test scalable mode
nextflow run PHemarajata/wf-assembly-snps-mod -profile test,docker --scalable_mode true

# Use the convenient wrapper script
./run_workflow.sh --input test_data/ --mode scalable --profile docker

# Test the tree grafting script (requires Python + Biopython)
pip install biopython
./graft_trees.py --help
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## 📄 License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Original [bacterial-genomics/wf-assembly-snps](https://github.com/bacterial-genomics/wf-assembly-snps) workflow
- [nf-core](https://nf-co.re/) community for workflow development best practices
- All the amazing bioinformatics tool developers whose software powers this workflow

## 📧 Support

For questions or support:
- Open an [issue](https://github.com/PHemarajata/wf-assembly-snps-mod/issues)
- Check the [documentation](docs/)
- Review the [usage examples](docs/usage.md)

---

**Citation**: If you use this workflow in your research, please cite the original tools and consider citing this repository.
