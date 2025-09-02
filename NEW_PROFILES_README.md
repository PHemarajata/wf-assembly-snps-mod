# New Nextflow Profiles

Two new profiles have been added to the wf-assembly-snps pipeline:

## 1. Local Workstation RTX 4070 Profile (`local_workstation_rtx4070`)

**Hardware Specifications:**
- 22 cores total (20 cores allocated, 2 cores reserved for system overhead)
- 64GB RAM total (56GB allocated, 8GB reserved for system overhead)
- NVIDIA GeForce RTX™ 4070 Laptop GPU

**Usage:**
```bash
nextflow run main.nf -profile docker,local_workstation_rtx4070 --input data/ --outdir results/
```

**Key Features:**
- Optimized for local workstation with moderate parallelization
- Conservative resource allocation to maintain system responsiveness
- GPU acceleration support prepared for future implementations
- Process-specific optimizations for clustering workflows

## 2. DGX Station A100 Updated Profile (`dgx_station_a100_updated`)

**Hardware Specifications:**
- 128 cores total (120 cores allocated, 8 cores reserved for system overhead)
- 512GB RAM total (480GB allocated, 32GB reserved for system overhead)
- NVIDIA A100 GPUs

**Usage:**
```bash
nextflow run main.nf -profile docker,dgx_station_a100_updated --input data/ --outdir results/
```

**Key Features:**
- High-performance configuration for large-scale analyses
- Aggressive parallelization with multiple concurrent processes
- Optimized for handling large genomic datasets
- GPU acceleration support for compatible processes
- Enhanced memory allocation for memory-intensive operations

## Resource Allocation Strategy

Both profiles implement conservative overhead management:

### Local Workstation RTX 4070:
- **CPU Overhead:** 2 cores (9% overhead)
- **Memory Overhead:** 8GB (12.5% overhead)
- **Max concurrent processes:** Limited to prevent system overload

### DGX Station A100 Updated:
- **CPU Overhead:** 8 cores (6.25% overhead)
- **Memory Overhead:** 32GB (6.25% overhead)
- **Max concurrent processes:** Higher limits for better throughput

## Process-Specific Optimizations

Both profiles include optimized settings for key processes:
- `MASH_SKETCH` and `MASH_DIST`: Distance calculation processes
- `CLUSTER_GENOMES`: Genome clustering
- `SKA_BUILD` and `SKA_ALIGN`: Split k-mer analysis
- `IQTREE_FAST`: Phylogenetic tree construction
- `GUBBINS_CLUSTER`: Recombination detection
- `USHER_BUILD` and `USHER_PLACE`: Phylogenetic placement

## Error Handling

Both profiles include:
- Automatic retry on common exit codes (71,104,134,137,139,140,143,255)
- Maximum 2 retries per process
- Graceful failure handling with unlimited error tolerance

## Future GPU Support

Both profiles are prepared for GPU acceleration with placeholder configurations for:
- RTX 4070: `nvidia-rtx-4070` accelerator type
- A100: `nvidia-a100` accelerator type

These can be activated when GPU-accelerated versions of bioinformatics tools become available.