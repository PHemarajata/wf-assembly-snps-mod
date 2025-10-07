# Changelog

All notable changes to this modified version of wf-assembly-snps are documented in this file.

## [1.0.3-mod] - 2024-10-07

### Added

#### 🚀 New Workflow Modes
- **Scalable Mode**: Divide-and-conquer approach for large datasets (hundreds to thousands of genomes)
  - Mash-based pre-clustering for distance estimation
  - Per-cluster analysis with SKA alignment and IQ-TREE2
  - Global integration with backbone phylogeny
  - Support for UShER incremental updates
- **Recombination-Aware Mode**: Enhanced analysis with Gubbins integration
  - Automatic recombination detection and masking
  - Improved phylogenetic accuracy for highly recombinogenic species

#### 📊 Performance Optimizations
- Clustering-based analysis for improved scalability
- Optimized Gubbins parameters (reduced iterations: 3 vs default 5)
- Hybrid tree building approach for better performance
- Resume optimization with checksum-based caching

#### 🛠️ Enhanced Configuration
- Comprehensive parameter system with clear defaults
- Multiple compute environment profiles
- Flexible alignment methods (Snippy, Parsnp)
- Configurable clustering thresholds and cluster sizes

#### 🌳 Standalone Tree Grafting Tool
- **`graft_trees.py`**: Standalone Python script for tree grafting when the main workflow fails
  - Robust leaf-expansion grafting algorithm using Biopython
  - Detailed logging and error reporting
  - Conflict resolution for tip label collisions
  - Dry-run capability for testing
  - Memory-efficient alternative to Nextflow-based grafting

#### 📁 New Parameters
- `--scalable_mode`: Enable scalable clustering workflow
- `--recombination_aware_mode`: Enable recombination detection
- `--mash_threshold`: Distance threshold for clustering (default: 0.028)
- `--max_cluster_size`: Maximum genomes per cluster (default: 50)
- `--merge_singletons`: Merge singleton clusters (default: true)
- `--mash_sketch_size`: Sketch size for large datasets (default: 50,000)
- `--gubbins_iterations`: Maximum Gubbins iterations (default: 3)
- `--integrate_results`: Integrate cluster results (default: true)
- `--build_usher_mat`: Build UShER mutation tree (default: false)

#### 🖥️ HPC Integration
- Pre-configured profiles for institutional HPC systems
- UGE (Univa Grid Engine) execution scripts
- Optimized resource allocation for different compute environments

#### 📖 Documentation
- Comprehensive README with parameter tables and examples
- Detailed scalable mode documentation
- Usage guides for different analysis modes
- Complete parameter reference with defaults

### Changed

#### 🔄 Workflow Architecture
- Modular workflow design with three distinct analysis modes
- Enhanced process organization and dependency management
- Improved error handling and logging

#### ⚙️ Default Parameters
- `recombination_aware_mode`: `true` (was `false`)
- `integrate_results`: `true` (new parameter)
- `mash_threshold`: `0.028` (optimized for bacterial genomes)
- `max_cluster_size`: `50` (balanced for performance)
- `merge_singletons`: `true` (improved clustering)
- `run_gubbins`: `true` (enabled by default)

#### 📊 Analysis Pipeline
- Enhanced SNP distance calculation with cluster integration
- Improved phylogenetic tree construction workflow
- Better handling of large datasets through clustering

#### 🏗️ Infrastructure
- Updated manifest information for modified version
- Enhanced container and environment management
- Improved reproducibility features

### Technical Details

#### Dependencies
- Nextflow ≥22.04.3
- Support for Docker, Singularity, and Conda
- Enhanced tool integration (Mash, SKA, IQ-TREE2, Gubbins, UShER)

#### Compatibility
- Maintains backward compatibility with original workflow parameters
- Supports all original input formats and methods
- Compatible with existing compute environment configurations

### Notes

This modified version builds upon the excellent foundation of the original 
[bacterial-genomics/wf-assembly-snps](https://github.com/bacterial-genomics/wf-assembly-snps) 
workflow, adding significant enhancements for scalability and recombination-aware analysis 
while maintaining full compatibility with existing usage patterns.