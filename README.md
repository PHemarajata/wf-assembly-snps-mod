# wf-assembly-snps-mod

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.04.3-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-blue.svg)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/singularity-blue.svg)](https://sylabs.io/docs/)
[![Conda](https://img.shields.io/badge/conda-green.svg)](https://docs.conda.io/en/latest/)

A modified and enhanced Nextflow workflow for bacterial genome assembly-based SNP identification and phylogenetic analysis.

> **Note**: This is a modified version of [bacterial-genomics/wf-assembly-snps](https://github.com/bacterial-genomics/wf-assembly-snps) with additional features and optimizations.

## 🚀 Quick Start

```bash
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

## 🧪 Testing

```bash
# Quick test with sample data
nextflow run PHemarajata/wf-assembly-snps-mod -profile test,docker

# Test scalable mode
nextflow run PHemarajata/wf-assembly-snps-mod -profile test,docker --scalable_mode true

# Use the convenient wrapper script
./run_workflow.sh --input test_data/ --mode scalable --profile docker
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